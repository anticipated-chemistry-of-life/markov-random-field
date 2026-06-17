#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

static constexpr uint32_t MAX_NUMBER_OF_MOLECULES = (1 << 24) - 1;
static constexpr size_t NUMBER_OF_TREES           = 2;
using DimensionIndex                              = std::array<size_t, NUMBER_OF_TREES>;
