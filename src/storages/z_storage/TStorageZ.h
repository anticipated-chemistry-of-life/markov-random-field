//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEZ_H
#define TSTORAGEZ_H

#include <cstdint>
/** A single Z cell: one state bit in one byte.
 *
 * A node state carries no MCMC counter -- only the field tracks a posterior fraction of ones -- so
 * a cell of one collapses to its state. Both node states hold their states as bare bytes and hand
 * this out when a caller asks for their stored cells, which is where a written node-state file
 * gets the state of a row from.
 */
class TStorageZ {
private:
	/// 0 = false (this is also what a cell no storage holds reads as), 1 = true.
	uint8_t _state = 0;

public:
	TStorageZ()  = default;
	~TStorageZ() = default;
	explicit TStorageZ(bool state) { set_state(state); }

	[[nodiscard]] bool is_one() const { return _state != 0; }
	void set_state(bool state) { _state = state ? 1 : 0; }
	void switch_state() { _state ^= 1; }

	/// "Empty" == what a cell no storage holds reads as (state false). Equivalent to
	/// *this == TStorageZ{}.
	[[nodiscard]] bool is_empty() const { return _state == 0; }

	bool operator==(const TStorageZ &other) const { return _state == other._state; }
	bool operator!=(const TStorageZ &other) const { return _state != other._state; }
};

static_assert(sizeof(TStorageZ) == 1);

#endif // TSTORAGEZ_H
