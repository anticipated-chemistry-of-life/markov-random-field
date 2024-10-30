//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEY_H
#define TSTORAGEY_H
#include <cstdint>
/** TStorage Y is the class to store a single value Y of the Random Markov Field
 * We store in a single 64 bits integer (8 bytes) the following information:
 * - the first 16 bits are the counter of the number of times the element was a one in the MCMC
 * - the 17th (position 16 ) is the current state of the element (0 or 1)
 * - the rest is the position of the element Y in the multidimensional linearized Y vector of leaves
 */
class TStorageY {
private:
	/** The value where we are going to store the state (0/1), the coodinate and the counter. */
	uint64_t _value = 0;

	/** to take the last 47 bits of our 64 bits integer we apply
	 * a mask of 17 zeros and 47 ones. This is done
	 * by shifting 1 to the left by 47 and subtracting 1
	 */
	static constexpr uint64_t _counter_mask    = ~((1ULL << 48) - 1);
	static constexpr uint64_t _state_mask      = (1ULL << 47);
	static constexpr uint64_t _coordinate_mask = (1ULL << 47) - 1;

public:
	TStorageY() = default;
	explicit TStorageY(const uint64_t coordinate) {
		// if we construct a TStorageY, that means that by default we will set the state to 1
		// with coordinate "coordinate"
		set_coordinate(coordinate);
		set_state(true);
	}

	/** @return The value of the 64 bits integer */
	[[nodiscard]] uint64_t value() const { return _value; };

	/** @return The coordinate of the element in the Y vector */
	uint64_t get_coordinate() const { return _value & _coordinate_mask; };

	/** Set the coordinate of the element in the Y vector
	 * @param coordinate The coordinate of the element in the Y vector
	 */
	void set_coordinate(const uint64_t coordinate) {
		_value = (_value & ~_coordinate_mask) | (coordinate & _coordinate_mask);
	}

	/** @return Wheter the element is one or zero. */
	bool is_one() const { return (_value & _state_mask) >> 47; }

	/** The state is stored in the 17th bit of our 64 bits integer.
	 * If set_state(true) then we set the 17th bit to 1.
	 * If set_state(false) then we set the 17th bit to 0.
	 */
	void set_state(const bool state) { _value = (_value & ~_state_mask) | (static_cast<uint64_t>(state)) << 47; }

	/** @return The counter of the number of times the element was a one in the MCMC */
	uint16_t get_counter() const { return _value >> 48; }

	/** Set the counter of the number of times the element was a one in the MCMC
	 * @param counter The counter of the number of times the element was a one in the MCMC
	 */
	void set_counter(const uint16_t counter) {
		_value = (~_counter_mask & _value) | (static_cast<uint64_t>(counter) << 48);
	}

	/** We update the counter if the element is a one */
	void update_counter() {
		if (is_one()) { set_counter(get_counter() + 1); }
	}

	/** We reset the counter to 0 (this will be used after the burn-in period of the MCMC) */
	void reset_counter() { set_counter(0); }

	bool operator<(const TStorageY &right) const { return get_coordinate() < right.get_coordinate(); }
	bool operator<(const uint64_t right) const { return get_coordinate() < right; }
	bool operator==(const uint64_t right) const { return get_coordinate() == right; }
	bool operator!=(const uint64_t right) const { return get_coordinate() != right; }
};

#endif // TSTORAGEY_H
