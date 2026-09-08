#pragma once

#include "common.h"

#include <functional>

// Output hook used by the 3-plex terminal. Pure installs it for every call so
// each terminal clique reaches the retained maximal-clique list.
using FastCliqueSink = std::function<void(const std::vector<ui> &)>;
