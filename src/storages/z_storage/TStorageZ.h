//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEZ_H
#define TSTORAGEZ_H

#include "coretools/Main/TError.h"
#include <cstdint>

class TStorageZ {
private:
	uint32_t _value                      = 0;
	/// bit equivalent : 1000 0000 0000 0000 0000 0000 0000 0000
	static constexpr uint32_t STATE_MASK = 0x80000000;
	/// bit equivalent : 0111 1111 1111 1111 1111 1111 1111 1111
	static constexpr uint32_t INDEX_MASK = 0x7FFFFFFF;

public:
	static constexpr uint32_t MAX_INDEX = INDEX_MASK;
	TStorageZ()                         = default;
	~TStorageZ()                        = default;
	explicit TStorageZ(const uint32_t linear_index_in_Z_space) {
		set_linear_index_in_Z_space(linear_index_in_Z_space);
		set_state(true);
	}
	[[nodiscard]] uint32_t get_linear_index_in_Z_space() const { return _value & INDEX_MASK; };
	[[nodiscard]] uint32_t get_linear_index_in_container_space() const {
		return get_linear_index_in_Z_space();
	};

	void set_linear_index_in_Z_space(const uint32_t linear_index_in_Z_space) {
		if (linear_index_in_Z_space > MAX_INDEX) {
			throw coretools::TDevError("counter exceeds 31-bit maximum (", MAX_INDEX, ")");
		}
		_value = (linear_index_in_Z_space & INDEX_MASK) | (_value & STATE_MASK);
	}

	[[nodiscard]] inline bool is_one() const { return (_value & STATE_MASK) != 0; }

	void set_state(bool state)
	// we want to set the state this is given
	// so if the state is false, the value should be negative
	// the state is true, we set the state to positive
	{
		_value = state ? (_value | STATE_MASK) : (_value & INDEX_MASK);
	}

	void switch_state() { _value ^= STATE_MASK; }

	bool operator<(const TStorageZ &right) const {
		return get_linear_index_in_Z_space() < right.get_linear_index_in_Z_space();
	}
	bool operator<(const uint32_t right) const { return get_linear_index_in_Z_space() < right; }
	bool operator==(const uint32_t right) const { return get_linear_index_in_Z_space() == right; }
	bool operator==(const TStorageZ &right) const {
		return get_linear_index_in_Z_space() < right.get_linear_index_in_Z_space();
	}
	bool operator!=(const uint32_t right) const { return get_linear_index_in_Z_space() != right; }
};

static_assert(sizeof(TStorageZ) == 4);

#endif // TSTORAGEZ_H
