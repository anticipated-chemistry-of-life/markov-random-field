# Commit and Comment Style

Commit messages and code comments use **Simplified Technical English (STE)** in spirit.

A commit message is read while bisecting. A code comment is read while the reader is already lost. Both are read under load, so both are written to be scanned.

## The rules

- **One idea per sentence.** Put the next idea in the next sentence.
- **Active voice.** Name who does the thing. "The block update writes the cell", not "the cell is written".
- **Short sentences.** About 20 words. A sentence that needs a comma to hold two clauses is two sentences.
- **Plain words.** Prefer the short word. Keep the domain terms exactly as `CONTEXT.md` spells them.
- **One word per thing.** The glossary decides which word. Use it every time.
- **Present tense for behaviour.** "The sparse window buffers the insert." Past tense is for what the commit changed.
- **State the defect, the cause, the repair.** In that order, and stop.

## Scope

The style applies to:

- every commit message from here on,
- new comments,
- comments in a file you are already rewriting.

Existing comments elsewhere stay as they are. A repo-wide rewrite is a large diff with no functional change, and it costs review attention that the code has a better claim on.

## Comments that argue belong in a record

A comment that states what the code does belongs in the code. A comment that argues a point — a derivation, a numerical-analysis result, a rejected alternative — belongs in an ADR under `docs/adr/`.

Leave two sentences and a pointer where the argument was. This keeps the header scannable and stops the same argument living in two places.

## Before and after

Taken from `be31746`, which is a good commit written at four times the length it needs.

**Before:**

> `get_fraction_of_ones` could return a number above 1. A cell that was a one at every counted iteration reported 1.000031 on a million-iteration chain, and that value went straight into `<prefix>_Y_posterior.txt`, where a reader is right to stop at a probability greater than one.
>
> The denominator was `n_iterations / thinning factor`, floored, while `add_to_counter` fires `ceil(n / t)` times over iterations 0 to n-1. Whenever `t` does not divide `n` those differ by exactly one. It needs `t > 1` to bite, so the default chain of twenty thousand never showed it.
>
> The obvious repair is to round the denominator up instead, and it is not enough. The counted iterations do not start at zero. …

**After:**

> `get_fraction_of_ones` returned a value above 1. A cell that was a one at every counted iteration reported 1.000031 on a million-iteration chain. The value went into the posterior file.
>
> The denominator was `n_iterations / thinning factor`, floored. `add_to_counter` fires `ceil(n / t)` times. The two differ by one when `t` does not divide `n`. The default chain hides it, because it needs `t > 1`.
>
> A ceiling does not repair it. Burn-in leaves the iteration index above zero, so the number of thinned iterations depends on where the chain started.
>
> `_total_counts` now rises inside `add_to_counter` and resets with the counters. The denominator is the numerator's ceiling by construction. Both backends change, because they must agree.

Same facts. Same order. Half the words, and every sentence carries one of them.

## What is not adopted

The **letter** of ASD-STE100: its approved word list, and the hard word counts.

The list has no entry for *posterior*, *sparse*, *Metropolis*, *clique* or *phylogeny*. STE handles this with a project technical-name list, which each project writes for itself. Writing one for this project is real work, and it is separable from the style above. Adopt the letter when that list exists.
