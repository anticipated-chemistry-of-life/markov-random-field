//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEZ_H
#define TSTORAGEZ_H
#include <cmath>
#include <cstdint>

class TStorageZ {
private:
	int32_t _coordinate = 0;

public:
	TStorageZ() = default;
	explicit TStorageZ(const int32_t coordinate) {
		set_coordinate(coordinate);
		set_state(true);
	}
	[[nodiscard]] uint32_t get_coordinate() const { return std::abs(_coordinate); };

	void set_coordinate(const int32_t coordinate) { _coordinate = coordinate; }

	[[nodiscard]] bool is_one() const { return !std::signbit(_coordinate); } // this should be ok
	void set_state(const bool state) {
		// we want to set the state this is given
		// so if the state is false, the coordinate should be negative
		// the state is true, we set the state to positive
		_coordinate = (state * 2 - 1) * std::abs(_coordinate);
	}

	bool operator<(const TStorageZ &right) const { return get_coordinate() < right.get_coordinate(); }
	bool operator<(const uint32_t right) const { return get_coordinate() < right; }
	bool operator==(const uint32_t right) const { return get_coordinate() == right; }
	bool operator!=(const uint32_t right) const { return get_coordinate() != right; }
};

#endif // TSTORAGEZ_H
