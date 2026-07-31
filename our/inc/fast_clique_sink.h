#pragma once

#include "common.h"

#include <functional>

// Optional validation/output hook used by the FastList family. Production
// enumeration leaves the function empty and retains the count-only hot path.
// When installed, every reported clique is sorted before it reaches the sink.
using FastCliqueSink = std::function<void(const std::vector<ui> &)>;

