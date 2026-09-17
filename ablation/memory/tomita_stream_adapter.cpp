#include <array>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

constexpr std::size_t kBufferSize = 1U << 20;

class FastInput {
public:
    explicit FastInput(const char* path) : input_(std::fopen(path, "rb")) {
        if (input_ == nullptr) {
            throw std::runtime_error(std::string("cannot open input: ") + path);
        }
    }

    ~FastInput() {
        if (input_ != nullptr) {
            std::fclose(input_);
        }
    }

    bool readUInt(std::uint32_t& value) {
        int c = nextChar();
        while (c != EOF && (c == ' ' || c == '\t' || c == '\r' || c == '\n')) {
            c = nextChar();
        }
        if (c == EOF) {
            return false;
        }
        if (c < '0' || c > '9') {
            throw std::runtime_error("input contains a non-numeric token");
        }

        std::uint64_t parsed = 0;
        do {
            parsed = parsed * 10 + static_cast<unsigned>(c - '0');
            if (parsed > std::numeric_limits<std::uint32_t>::max()) {
                throw std::runtime_error("vertex identifier exceeds uint32_t");
            }
            c = nextChar();
        } while (c >= '0' && c <= '9');

        if (c != EOF && c != ' ' && c != '\t' && c != '\r' && c != '\n') {
            throw std::runtime_error("input contains an invalid delimiter");
        }
        value = static_cast<std::uint32_t>(parsed);
        return true;
    }

private:
    int nextChar() {
        if (position_ == available_) {
            available_ = std::fread(buffer_.data(), 1, buffer_.size(), input_);
            position_ = 0;
            if (available_ == 0) {
                if (std::ferror(input_) != 0) {
                    throw std::runtime_error("failed while reading input");
                }
                return EOF;
            }
        }
        return static_cast<unsigned char>(buffer_[position_++]);
    }

    std::FILE* input_ = nullptr;
    std::array<char, kBufferSize> buffer_ {};
    std::size_t position_ = 0;
    std::size_t available_ = 0;
};

class FastOutput {
public:
    ~FastOutput() {
        flush();
    }

    void writeUInt(std::uint64_t value) {
        char digits[32];
        std::size_t length = 0;
        do {
            digits[length++] = static_cast<char>('0' + value % 10);
            value /= 10;
        } while (value != 0);
        while (length != 0) {
            put(digits[--length]);
        }
    }

    void put(char value) {
        if (position_ == buffer_.size()) {
            flush();
        }
        buffer_[position_++] = value;
    }

    void writeArc(std::uint32_t source, std::uint32_t target) {
        writeUInt(source);
        put(',');
        writeUInt(target);
        put('\n');
    }

    void flush() {
        if (position_ == 0) {
            return;
        }
        if (std::fwrite(buffer_.data(), 1, position_, stdout) != position_) {
            std::exit(EXIT_FAILURE);
        }
        position_ = 0;
    }

private:
    std::array<char, kBufferSize> buffer_ {};
    std::size_t position_ = 0;
};

std::uint64_t parseArgument(const char* text, const char* name) {
    errno = 0;
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') {
        throw std::runtime_error(std::string("invalid ") + name + ": " + text);
    }
    return value;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "Usage: %s EDGE_LIST NUM_VERTICES NUM_EDGES\n", argv[0]);
        return 2;
    }

    try {
        const std::uint64_t n = parseArgument(argv[2], "NUM_VERTICES");
        const std::uint64_t m = parseArgument(argv[3], "NUM_EDGES");
        if (n > std::numeric_limits<std::uint32_t>::max() ||
            m > std::numeric_limits<std::uint64_t>::max() / 2) {
            throw std::runtime_error("graph dimensions exceed adapter limits");
        }

        FastInput input(argv[1]);
        FastOutput output;
        output.writeUInt(n);
        output.put('\n');
        output.writeUInt(2 * m);
        output.put('\n');

        for (std::uint64_t edge = 0; edge < m; ++edge) {
            std::uint32_t u = 0;
            std::uint32_t v = 0;
            if (!input.readUInt(u) || !input.readUInt(v)) {
                throw std::runtime_error("edge list ended before its declared edge count");
            }
            if (u >= n || v >= n || u >= v) {
                throw std::runtime_error("edge list is not normalized as 0<=u<v<n");
            }
            output.writeArc(u, v);
            output.writeArc(v, u);
        }

        std::uint32_t extra = 0;
        if (input.readUInt(extra)) {
            throw std::runtime_error("edge list contains data beyond its declared edge count");
        }
        output.flush();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "rmce_stream_adapter: %s\n", error.what());
        return 1;
    }
    return 0;
}
