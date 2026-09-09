//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEY_H
#define TSTORAGEY_H

#include "coretools/Main/TError.h"
#include <cstdint>
#include <cstdlib>

/** One cell of the field: the state, and how often the chain counted that cell a one, packed into
 * a single 16-bit word.
 *
 * - bit 15 is the state (0 or 1),
 * - bits 0..14 are the posterior counter, so it holds 32767 counted iterations.
 *
 * The cell knows nothing of where it is. Both fields key their cells by linear index, so the
 * position is the container's business and not the cell's.
 *
 * Both backends hold this cell, which is what makes their posterior fields comparable: a chain of
 * n iterations is thinned to one iteration in ceil(n / 32767) whichever one runs it. The counter
 * was 16 bits wide on the dense side while the two fields held different cells, and a posterior
 * field written then is at twice the resolution of one written now.
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
	/// "Empty" == what a cell no storage holds reads as: state == false AND counter == 0 (i.e.
	/// _value == 0). Equivalent to *this == TStorageY{}.
	[[nodiscard]] bool is_empty() const { return _value == 0; }
};
static_assert(sizeof(TStorageY) == 2);

#endif // TSTORAGEY_H
