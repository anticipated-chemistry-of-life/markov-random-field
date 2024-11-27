//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEZ_H
#define TSTORAGEZ_H
#include <cmath>
#include <cstdint>

class TStorageZ {
private:
	int32_t _value = 0;

public:
	TStorageZ() = default;
	explicit TStorageZ(const int32_t linear_index_in_Z_space) {
		set_linear_index_in_Z_space(linear_index_in_Z_space);
		set_state(true);
	}
	[[nodiscard]] uint32_t get_linear_index_in_Z_space() const { return std::abs(_value); };

	void set_linear_index_in_Z_space(const int32_t linear_index_in_Z_space) { _value = linear_index_in_Z_space; }

	[[nodiscard]] bool is_one() const { return !std::signbit(_value); } // this should be ok
	void set_state(const bool state) {
		// we want to set the state this is given
		// so if the state is false, the value should be negative
		// the state is true, we set the state to positive
		_value = (state * 2 - 1) * std::abs(_value);
	}

	void switch_state() { _value = -_value; }

	bool operator<(const TStorageZ &right) const {
		return get_linear_index_in_Z_space() < right.get_linear_index_in_Z_space();
	}
	bool operator<(const uint32_t right) const { return get_linear_index_in_Z_space() < right; }
	bool operator==(const uint32_t right) const { return get_linear_index_in_Z_space() == right; }
	bool operator!=(const uint32_t right) const { return get_linear_index_in_Z_space() != right; }
};

#endif // TSTORAGEZ_H
