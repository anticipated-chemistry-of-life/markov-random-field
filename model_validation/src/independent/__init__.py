"""An independent reference implementation of the latent field.

This package deliberately does *not* share code with the C++ binary it is used to
check. It re-derives the model from the maths, and only mirrors the C++ where a
file format or an index ordering has to agree exactly (see `indexing` and `io`).

Both trees are active. Each draws its own node state top-down over every one of
its nodes, and the field is a noisy AND of the two leaf blocks -- see
`docs/adr/0005-each-tree-owns-its-leaf-level-field.md`, which is what makes an
exact reference possible without neutralising anything.
"""
