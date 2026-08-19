"""Separate the two sources of downward bias in log_nu.

Fit a single scalar shift delta with log_nu_c = true_log_nu_c + delta, POOLED over all 128
cliques (so per-clique small-sample bias averages out). Truth is delta = 0.
  * if the product objective still gives delta << 0  -> the objective itself is wrong
  * if the pseudo-likelihood gives delta ~ 0         -> the missing normaliser is the cause
"""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

d = np.load("/tmp/replicated_states.npz")
S, M, sA, sN, mA, mN = d["S"], d["M"], d["sA"], d["sN"], d["mA"], d["mN"]
sPar, mPar, sLeaf, mLeaf, sRoot = d["sPar"], d["mPar"], d["sLeaf"], d["mLeaf"], d["sRoot"]
nS, L = S.shape
sLeafPos = {int(v): i for i, v in enumerate(sLeaf)}
kids = [[] for _ in range(nS)]
for n in range(nS):
    if n not in sRoot: kids[sPar[n]].append(n)

def logmats(nus):
    return np.log(np.array([expm(np.array([[-a*n, a*n], [(1-a)*n, (a-1)*n]]))
                            for a, n in zip(sA, nus)]))
lPm = np.log(np.array([expm(np.array([[-a*n, a*n], [(1-a)*n, (a-1)*n]]))
                       for a, n in zip(mA, mN)]))
cl = np.arange(L)
Si = S.astype(int)
# molecule-parent state for each (species leaf k, clique l)
mo_par = np.array([M[mPar[v]] for v in mLeaf]).astype(int).T   # (k, l)

def neg_prod(delta):
    lm = logmats(sN*np.exp(delta)); tot = 0.0
    for n in range(nS):
        if n in sRoot: continue
        tot += lm[cl, Si[sPar[n]], Si[n]].sum()
    return -tot

def neg_pl(delta):
    lm = logmats(sN*np.exp(delta)); tot = 0.0
    for n in range(nS):
        lp = np.stack([np.log(1-sA), np.log(sA)]) if n in sRoot else lm[cl, Si[sPar[n]], :].T
        lp = lp.copy()
        for c in kids[n]: lp += lm[cl, :, Si[c]].T
        if n in sLeafPos:
            k = sLeafPos[n]
            lp += lPm[k, mo_par[k], :].T
        tot += (lp[Si[n], cl] - np.logaddexp(lp[0], lp[1])).sum()
    return -tot

for name, f in [("product objective (current C++ code)", neg_prod),
                ("pseudo-likelihood (normalised conds)", neg_pl)]:
    r = minimize_scalar(f, bounds=(-4, 3), method="bounded")
    print(f"{name:38s}  delta_hat = {r.x:+.4f}   (truth 0)")
