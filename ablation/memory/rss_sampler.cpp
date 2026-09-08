#include <algorithm>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {

volatile std::sig_atomic_t interruptedSignal = 0;

extern "C" void rememberSignal(int signal) { interruptedSignal = signal; }

struct Options {
  std::uint64_t intervalUs = 0;
  std::uint64_t timeoutSeconds = 0;
  std::string outputPath;
  int commandIndex = -1;
};

struct MemorySample {
  std::uint64_t rssKb = 0;
  std::uint64_t hwmKb = 0;
};

void usage(std::ostream &out, const char *program) {
  out << "Usage: " << program
      << " --interval-us N --timeout-seconds N --output FILE -- COMMAND [ARGS...]\n";
}

std::uint64_t parsePositive(const char *text, const char *option) {
  if (text == nullptr || *text == '\0') {
    throw std::invalid_argument(std::string("missing value for ") + option);
  }
  std::uint64_t value = 0;
  const char *end = text + std::strlen(text);
  const auto parsed = std::from_chars(text, end, value);
  if (parsed.ec != std::errc{} || parsed.ptr != end || value == 0 ||
      value > static_cast<std::uint64_t>(
                  std::numeric_limits<std::int64_t>::max())) {
    throw std::invalid_argument(std::string(option) +
                                " must be a positive integer");
  }
  return value;
}

Options parseOptions(int argc, char **argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "-h" || argument == "--help") {
      usage(std::cout, argv[0]);
      std::exit(0);
    }
    if (argument == "--") {
      options.commandIndex = i + 1;
      break;
    }
    if (argument != "--interval-us" && argument != "--timeout-seconds" &&
        argument != "--output") {
      throw std::invalid_argument("unknown option: " + argument);
    }
    if (++i >= argc) {
      throw std::invalid_argument("missing value after " + argument);
    }
    if (argument == "--interval-us") {
      options.intervalUs = parsePositive(argv[i], "--interval-us");
    } else if (argument == "--timeout-seconds") {
      options.timeoutSeconds = parsePositive(argv[i], "--timeout-seconds");
    } else {
      options.outputPath = argv[i];
    }
  }
  if (options.intervalUs == 0 || options.timeoutSeconds == 0 ||
      options.outputPath.empty() || options.commandIndex < 0 ||
      options.commandIndex >= argc) {
    throw std::invalid_argument("interval, timeout, output, and command are required");
  }
  return options;
}

std::optional<MemorySample> readMemory(pid_t pid) {
  std::ifstream status("/proc/" + std::to_string(pid) + "/status");
  if (!status) {
    return std::nullopt;
  }
  MemorySample sample;
  bool foundRss = false;
  bool foundHwm = false;
  std::string key;
  while (status >> key) {
    if (key == "VmRSS:" || key == "VmHWM:") {
      std::uint64_t value = 0;
      std::string unit;
      if (!(status >> value >> unit) || unit != "kB") {
        return std::nullopt;
      }
      if (key == "VmRSS:") {
        sample.rssKb = value;
        foundRss = true;
      } else {
        sample.hwmKb = value;
        foundHwm = true;
      }
    } else {
      std::string remainder;
      std::getline(status, remainder);
    }
  }
  if (!foundRss) {
    return std::nullopt;
  }
  if (!foundHwm) {
    sample.hwmKb = sample.rssKb;
  }
  return sample;
}

void signalProcessGroup(pid_t pid, int signal) {
  if (kill(-pid, signal) == 0) {
    return;
  }
  if (kill(pid, signal) != 0 && errno != ESRCH) {
    std::cerr << "warning: failed to signal child: " << std::strerror(errno)
              << '\n';
  }
}

void installSignalHandlers() {
  struct sigaction action {};
  action.sa_handler = rememberSignal;
  sigemptyset(&action.sa_mask);
  for (const int signal : {SIGINT, SIGTERM, SIGHUP}) {
    if (sigaction(signal, &action, nullptr) != 0) {
      throw std::runtime_error("failed to install signal handler");
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    const Options options = parseOptions(argc, argv);
    const int traceFd = open(options.outputPath.c_str(),
                             O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
    if (traceFd < 0) {
      throw std::runtime_error("cannot open trace file: " + options.outputPath);
    }
    FILE *trace = fdopen(traceFd, "w");
    if (trace == nullptr) {
      close(traceFd);
      throw std::runtime_error("cannot create trace stream");
    }
    std::fprintf(trace, "# finalCode RSS trace v1\n");
    std::fprintf(trace, "# interval_us=%llu\n",
                 static_cast<unsigned long long>(options.intervalUs));
    std::fprintf(trace, "# columns=elapsed_us rss_kb vm_hwm_kb\n");
    std::fflush(trace);

    int execPipe[2];
    if (pipe2(execPipe, O_CLOEXEC) != 0) {
      std::fclose(trace);
      throw std::runtime_error("cannot create exec synchronization pipe");
    }

    const pid_t child = fork();
    if (child < 0) {
      close(execPipe[0]);
      close(execPipe[1]);
      std::fclose(trace);
      throw std::runtime_error("fork failed");
    }
    if (child == 0) {
      close(execPipe[0]);
      setpgid(0, 0);
      std::vector<char *> command;
      for (int i = options.commandIndex; i < argc; ++i) {
        command.push_back(argv[i]);
      }
      command.push_back(nullptr);
      execvp(command[0], command.data());
      const int execError = errno;
      const ssize_t ignored = write(execPipe[1], &execError, sizeof(execError));
      static_cast<void>(ignored);
      std::cerr << "exec failed: " << std::strerror(execError) << '\n';
      _exit(127);
    }

    close(execPipe[1]);
    if (setpgid(child, child) != 0 && errno != EACCES && errno != ESRCH) {
      signalProcessGroup(child, SIGKILL);
      close(execPipe[0]);
      std::fclose(trace);
      throw std::runtime_error("failed to create child process group");
    }

    int execError = 0;
    ssize_t execResult;
    do {
      execResult = read(execPipe[0], &execError, sizeof(execError));
    } while (execResult < 0 && errno == EINTR);
    close(execPipe[0]);
    if (execResult > 0) {
      int status = 0;
      waitpid(child, &status, 0);
      std::fclose(trace);
      throw std::runtime_error(std::string("child exec failed: ") +
                               std::strerror(execError));
    }
    if (execResult < 0) {
      signalProcessGroup(child, SIGKILL);
      waitpid(child, nullptr, 0);
      std::fclose(trace);
      throw std::runtime_error("failed to synchronize with child exec");
    }

    installSignalHandlers();
    using Clock = std::chrono::steady_clock;
    const auto interval = std::chrono::microseconds(options.intervalUs);
    const auto start = Clock::now();
    const auto timeoutDeadline =
        start + std::chrono::seconds(options.timeoutSeconds);
    auto nextSample = start;
    auto killDeadline = Clock::time_point::max();
    bool timedOut = false;
    bool terminationSent = false;
    bool forceKillSent = false;
    int forwardedSignal = 0;
    int childStatus = 0;
    struct rusage usage {};
    std::uint64_t samples = 0;
    std::uint64_t peakSampledRssKb = 0;
    std::uint64_t peakObservedHwmKb = 0;

    while (true) {
      auto now = Clock::now();
      if (now >= nextSample) {
        if (const auto sample = readMemory(child)) {
          const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
                                   now - start)
                                   .count();
          std::fprintf(trace, "%lld\t%llu\t%llu\n",
                       static_cast<long long>(elapsed),
                       static_cast<unsigned long long>(sample->rssKb),
                       static_cast<unsigned long long>(sample->hwmKb));
          ++samples;
          peakSampledRssKb = std::max(peakSampledRssKb, sample->rssKb);
          peakObservedHwmKb = std::max(peakObservedHwmKb, sample->hwmKb);
          if (samples % 256U == 0U) {
            std::fflush(trace);
          }
        }
        const auto afterSample = Clock::now();
        nextSample += interval;
        if (nextSample <= afterSample) {
          nextSample = afterSample + interval;
        }
      }

      const pid_t waited = wait4(child, &childStatus, WNOHANG, &usage);
      if (waited == child) {
        break;
      }
      if (waited < 0 && errno != EINTR) {
        signalProcessGroup(child, SIGKILL);
        waitpid(child, nullptr, 0);
        std::fclose(trace);
        throw std::runtime_error("wait4 failed");
      }

      now = Clock::now();
      if (interruptedSignal != 0 && !terminationSent) {
        forwardedSignal = interruptedSignal;
        signalProcessGroup(child, forwardedSignal);
        terminationSent = true;
        killDeadline = now + std::chrono::seconds(10);
      }
      if (!timedOut && now >= timeoutDeadline) {
        timedOut = true;
        signalProcessGroup(child, SIGTERM);
        terminationSent = true;
        killDeadline = now + std::chrono::seconds(10);
      }
      if (terminationSent && !forceKillSent && now >= killDeadline) {
        signalProcessGroup(child, SIGKILL);
        forceKillSent = true;
      }

      auto wake = nextSample;
      if (!timedOut) {
        wake = std::min(wake, timeoutDeadline);
      }
      if (terminationSent && !forceKillSent) {
        wake = std::min(wake, killDeadline);
      }
      if (wake > Clock::now()) {
        std::this_thread::sleep_until(wake);
      }
    }

    const auto elapsedUs =
        std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - start)
            .count();
    const std::uint64_t wait4PeakRssKb =
        usage.ru_maxrss > 0 ? static_cast<std::uint64_t>(usage.ru_maxrss) : 0U;
    const int childExitCode = WIFEXITED(childStatus) ? WEXITSTATUS(childStatus) : -1;
    const int childSignal = WIFSIGNALED(childStatus) ? WTERMSIG(childStatus) : 0;

    std::fprintf(trace, "# samples=%llu\n",
                 static_cast<unsigned long long>(samples));
    std::fprintf(trace, "# peak_sampled_rss_kb=%llu\n",
                 static_cast<unsigned long long>(peakSampledRssKb));
    std::fprintf(trace, "# peak_observed_hwm_kb=%llu\n",
                 static_cast<unsigned long long>(peakObservedHwmKb));
    std::fprintf(trace, "# wait4_peak_rss_kb=%llu\n",
                 static_cast<unsigned long long>(wait4PeakRssKb));
    std::fprintf(trace, "# elapsed_us=%lld\n", static_cast<long long>(elapsedUs));
    std::fprintf(trace, "# child_exit_code=%d\n", childExitCode);
    std::fprintf(trace, "# child_signal=%d\n", childSignal);
    std::fprintf(trace, "# timed_out=%d\n", timedOut ? 1 : 0);
    std::fclose(trace);

    std::cout << "sampler.trace=" << options.outputPath << '\n'
              << "sampler.interval_us=" << options.intervalUs << '\n'
              << "sampler.samples=" << samples << '\n'
              << "sampler.peak_sampled_rss_kb=" << peakSampledRssKb << '\n'
              << "sampler.peak_observed_hwm_kb=" << peakObservedHwmKb << '\n'
              << "sampler.wait4_peak_rss_kb=" << wait4PeakRssKb << '\n'
              << "sampler.elapsed_us=" << elapsedUs << '\n'
              << "sampler.timed_out=" << (timedOut ? 1 : 0) << '\n'
              << "sampler.child_exit_code=" << childExitCode << '\n'
              << "sampler.child_signal=" << childSignal << '\n';

    if (timedOut) {
      return 124;
    }
    if (forwardedSignal != 0) {
      return 128 + forwardedSignal;
    }
    if (childSignal != 0) {
      return 128 + childSignal;
    }
    return childExitCode;
  } catch (const std::exception &error) {
    std::cerr << "rss_sampler: " << error.what() << '\n';
    usage(std::cerr, argv[0]);
    return 2;
  }
}
