# Mass-spec scores enter only relative to their feature's best row

The mass-spec file gives every candidate of a feature a calibration score (CONTEXT.md, "Calibration score"). A calibration score is a density under a kernel density estimate, fitted on correct identifications, over the CSI:FingerID score and the feature's ion mass. Its absolute scale depends on the bandwidth and on the ion mass of the feature. Two features' scores are therefore not comparable, and a score alone says nothing.

The model reads a calibration score only as a **relative calibration score**: the score divided by the best score among all of its feature's rows (`src/mass_spec/mass_spec_file.cpp`). The best row of a feature is 1. The unknown molecule gets one fixed relative calibration score, the unknown fraction `c` (`--ms_unknown_fraction`, default 0.1). So `c = 0.1` reads "the unknown molecule is one tenth as likely as the feature's best row", in every feature.

Nothing else is lost. Every assignment move compares two explanations of one feature (`src/mass_spec/msms_run.cpp`). A ratio between two scores of one feature does not change when both are divided by the same reference. The only number the division changes is the one the scores are compared against: the unknown molecule.

## The reference includes skipped candidates

A candidate that is not a molecule leaf is skipped, but its score still sets its feature's reference.

Take a feature whose best row is a molecule outside the molecule tree. That row is evidence that the feature is something the tree does not hold, which is what the unknown molecule stands for. Without that row in the reference, the best leaf candidate would take bin 0. It would then look as good as the true best explanation, and the unknown molecule would lose to it more often than the data says. With the row in the reference, every leaf candidate keeps its true ratio to the best explanation.

The same holds for SIRIUS's own cut. SIRIUS keeps at most 200 candidates per feature, so the file never holds the rest. The loader cannot see them, and the load report says so.

## The 1/16-nat grid

A relative calibration score is stored as a bin `b` in `0 … 255`, and bin `b` stands for a score whose log is `-b/16` (`RELATIVE_CALIBRATION_SCORE`, `src/constants.h`). A candidate's bin is `round(-16 ln(score / reference))`. The unknown molecule's is `round(-16 ln c)`, so `c = 0.1` is bin 37 and `c = 1` is bin 0.

The bin has 8 bits because `TFeatureLikelihood` packs it with a 24-bit molecule index into 4 bytes. 256 bins at 1/16 nat reach down to `e^(-255/16)`, about `1.2e-7`. Rounding moves a ratio by at most 1/32 nat, about 3%. A score below the floor is clamped to bin 255 and counted in the load report. Such a candidate is almost seven orders of magnitude below its feature's best row, so the clamp moves its chance of being assigned by a negligible amount.

Bin 0 is the reference and is valid. The old table (`PHRED_LIKE_PROBA`) mapped bin 0 to probability 0, so `TMassSpecRun` rejected bin 0 for a candidate and for the unknown molecule. Both guards are gone, with the tests that pinned them.

## Considered options

**Raw densities.** The score would enter as it is. It was rejected because a single `c` would then mean a different thing in every feature: a feature with high densities would almost never take the unknown molecule, and one with low densities almost always would. Binning raw densities would also need one global scale, over values that span many orders of magnitude.

**A per-feature unknown column.** The file would carry, per feature, the score of the unknown molecule. It was rejected because the prepared file carries no such number, and producing one is upstream work. One global `c` is one number to choose and to vary, instead of a column to invent.

**A likelihood ratio against an incorrect-identification KDE.** A second density, fitted on incorrect identifications, would divide the first. The ratio would then have a meaning on its own, and the unknown molecule could sit at a principled value. It was rejected for now because the file carries one density, fitted on correct identifications, and a second fit is upstream work. If that fit exists one day, the loader changes and the model does not: the moves already read only ratios within a feature.

## Consequences

`mass_spec_file::load` is a pure function of the file, the leaf names and `c`, and its tests build fixtures and read the runs back.

The feature term of the mass-spec likelihood (`mass_spec_math::feature_term`) is a relative calibration score, not a probability. It cancels from every field draw and from every filter-probability move. The contamination move reads it, through the whole log-likelihood, but it does not depend on the contamination either, so it cancels there too.
