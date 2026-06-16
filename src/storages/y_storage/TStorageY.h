//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEY_H
#define TSTORAGEY_H

#include "coretools/Main/TError.h"
#include <cstdint>
#include <cstdlib>
/** TStorage Y is the class to store a single value Y of the Random Markov Field
 * We store in a single 64 bits integer (8 bytes) the following information:
 * - the first 16 bits are the counter of the number of times the element was a one in the MCMC
 * - the 17th (position 16 ) is the current state of the element (0 or 1)
 * - the rest is the linear index in the Y space.
 */

class TStorageY {
private:
	uint16_t _value = 0; // bit 15 = state, bits 0..14 = counter

	static constexpr uint16_t STATE_MASK   = 0x8000; // 1000 0000 0000 0000
	static constexpr uint16_t COUNTER_MASK = 0x7FFF; // 0111 1111 1111 1111

public:
	static constexpr uint16_t MAX_COUNTER = COUNTER_MASK; // 32767

	TStorageY() = default;
	explicit TStorageY(bool state) { set_state(state); } // counter starts at 0

	[[nodiscard]] uint16_t value() const { return _value; }

	[[nodiscard]] bool is_one() const { return (_value & STATE_MASK) != 0; }
	void set_state(bool state) { _value = state ? (_value | STATE_MASK) : (_value & COUNTER_MASK); }
	void switch_state() { _value ^= STATE_MASK; }

	[[nodiscard]] uint16_t get_counter() const { return _value & COUNTER_MASK; }
	void set_counter(uint16_t counter) {
		if (counter > MAX_COUNTER) {
			throw coretools::TDevError("counter exceeds 15-bit maximum (", MAX_COUNTER, ")");
		}
		_value = (_value & STATE_MASK) | counter;
	}
	void update_counter() {
		if (is_one()) { set_counter(get_counter() + 1); }
	}
	void reset_counter() { _value &= STATE_MASK; } // clears counter, keeps state

	bool operator==(const TStorageY &other) const { return _value == other._value; }
	bool operator!=(const TStorageY &other) const { return _value != other._value; }
	/// "Empty" == the sentinel the sparse matrix uses for an absent cell:
	/// state == false AND counter == 0 (i.e. _value == 0). Equivalent to *this == TStorageY{}.
	[[nodiscard]] bool is_empty() const { return _value == 0; }
};
static_assert(sizeof(TStorageY) == 2);

#endif // TSTORAGEY_H
