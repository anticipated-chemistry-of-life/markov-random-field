//
// Created by VISANI Marco on 17.10.2024.
//

#ifndef TSTORAGEZ_H
#define TSTORAGEZ_H

#include <cstdint>
/** A single Z cell. Since Z migrated to a TSparseMatrix (like Y), the (row, col)
 * position encodes the linear index, so the element no longer needs to store its own
 * index. Z also has no MCMC counter (only Y tracks a posterior fraction of ones), so
 * the cell collapses to a single state bit stored in one byte.
 */
class TStorageZ {
private:
	/// 0 = false (this is also the "absent" sentinel TSparseMatrix uses), 1 = true.
	uint8_t _state = 0;

public:
	TStorageZ()  = default;
	~TStorageZ() = default;
	explicit TStorageZ(bool state) { set_state(state); }

	[[nodiscard]] bool is_one() const { return _state != 0; }
	void set_state(bool state) { _state = state ? 1 : 0; }
	void switch_state() { _state ^= 1; }

	/// "Empty" == the sentinel the sparse matrix uses for an absent cell (state false).
	/// Equivalent to *this == TStorageZ{}.
	[[nodiscard]] bool is_empty() const { return _state == 0; }

	bool operator==(const TStorageZ &other) const { return _state == other._state; }
	bool operator!=(const TStorageZ &other) const { return _state != other._state; }
};

static_assert(sizeof(TStorageZ) == 1);

#endif // TSTORAGEZ_H
