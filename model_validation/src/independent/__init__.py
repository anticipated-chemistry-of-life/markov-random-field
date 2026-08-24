"""An independent reference implementation of the latent field.

This package deliberately does *not* share code with the C++ binary it is used to
check. It re-derives the model from the maths, and only mirrors the C++ where a
file format or an index ordering has to agree exactly (see `indexing` and `io`).

See `docs/adr/0001-neutralise-molecules-dimension-for-validation.md` for why a
single-tree reference is enough to check a two-tree model.
"""
