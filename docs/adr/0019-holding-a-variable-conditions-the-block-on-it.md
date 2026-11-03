# Holding a variable conditions the block on it

_Supersedes in part [ADR-0010](0010-the-block-draws-the-leaf-layer.md): its paragraph on what `--fix_Y` and `--fix_Z` hold. The block update and everything else ADR-0010 decides stand._

**`--fix_Y` holds the field alone, and `--fix_Z` holds every node of both node states, the leaf blocks included.** The block update still draws whatever the flags leave of a leaf pair, given what they hold. **A node state the run supplies (`--set_<tree>_Z`) is kept at the chain start, its leaf block included**, rather than having its tree field overwritten by the field. These are ADR-0008's meanings for the two flags, brought back under the block.

## Why

ADR-0010 said that holding the field held both tree fields too, because one draw owns all three. That does not follow. Holding one variable of the block is conditioning on it: the other two still have a distribution, `p(z_s, z_m | y, parents)`. ADR-0010's reading froze both tree fields at the field for the whole chain, and the chain start put them there even when the run had supplied them.

Issue #106 is what that cost. The independent reference's rung 1 pins the field and both true node states, leaves included, so that the moves on alpha, nu and the branch lengths are scored against the true configuration. Under ADR-0010 the chain held true internal nodes over leaves equal to the field. The field is a corrupted AND of the two tree fields, about 23% ones where the true tree fields are about 46%. Alpha coverage fell to about 25% (bias −0.27), log nu was biased by +0.4, and the leaf branches went to the top bin. An exact grid posterior over each clique's alpha and log nu matched the C++ to a correlation of 0.99 once it was given the same leaves, and covered at 91–94% when given the true ones. The sampler was right; what it conditioned on was not what the rung supplied. Rung 2 (field pinned, node states free) inherited the same frozen leaves.

To reproduce the before: `rung1_pin_field_and_states.sh` from `simulate_independent.py --mass_spec dense` (seed 42) at `9e3ab28`, scored by `validate_independent.py`. Issue #106 records the grid posterior and its numbers.

## How

A held variable enters the eight-state draw as the indicator of the state it holds (`block_update::factors_given_held`). A held field replaces the cell data factor by `{y == 0, y == 1}`. Held tree fields replace both tree factors by their own states. Every block state that would move a held variable then has no mass, so the kernel ADR-0005 derived and `TFieldMath_Tests.cpp` checks draws the conditional with no new arithmetic. Each data source is still told the model's own factors and the state it ended in, so its tally is of the field the update leaves, held or not.

When both are held there is nothing left to draw, and no block update runs. The tally the chain start built stands for the whole chain, as ADR-0010 had it for `--fix_Y`.

## Consequences

**A run that held only the field and ran no block update now runs one.** It costs what an unheld block costs, data passes included. The passes are what keep each data source's tally on the field.

**No flag combination reproduces the old `--fix_Y`.** It held the field and both tree fields at the field while the nodes above them moved. `--fix_Y --fix_Z` with no `--set_<tree>_Z` comes closest: the chain start puts both tree fields at the field and builds the internal nodes from them, and then holds every node too. A run that wants leaves at the field and free internal nodes has no flag for it, which is deliberate: that configuration is the defect #106 found. The constructor's log lines say what each flag holds.

**A node-state file that names only internal nodes now leaves its leaf rows at zero.** The chain start used to fill them from the field. Under `--fix_Z` those zeros are held; otherwise the first block update draws the leaf layer. Such files predate ADR-0005 and are warm starts at most.

## Considered options

**Score rungs 1 and 2 against the leaf layer they actually hold.** Rejected. It would make the reference agree with the C++ by describing a configuration that no unpinned chain ever reaches, and it would leave the ladder unable to isolate the tree process from the link.

**Hold the leaf layer only when a node state is supplied, and keep ADR-0010's `--fix_Y`.** Rejected. Rung 2 would still freeze both tree fields at the field, which is the same defect one rung lower.
