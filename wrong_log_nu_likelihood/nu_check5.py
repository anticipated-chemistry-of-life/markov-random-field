"""Same small two-tree MRF as before, but now also evaluate the maximum-pseudo-likelihood
estimator (normalised full conditionals) to check it is consistent."""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

rng = np.random.default_rng(7)

def P(alpha, nu, t=1.0):
    L = np.array([[-alpha*nu, alpha*nu], [(1-alpha)*nu, (alpha-1)*nu]])
    return expm(L*t)

def balanced(depth):
    n = 2**(depth+1) - 1
    parent = [-1] + [(i-1)//2 for i in range(1, n)]
    leaves = [i for i in range(n) if 2*i+1 >= n]
    internal = [i for i in range(1, n) if 2*i+1 < n]
    return n, parent, leaves, internal

DEPTH = 5
nS, parS, leafS, intS = balanced(DEPTH)
nM, parM, leafM, intM = balanced(DEPTH)
alpha = 0.5
TRUE = -1.5
Ms_true, Mm_true = P(alpha, np.exp(TRUE)), P(alpha, np.exp(TRUE))
sIdx = {n: k for k, n in enumerate(leafS)}
mIdx = {n: k for k, n in enumerate(leafM)}

Zs = rng.random((nS, len(leafM))) < 0.5
Zm = rng.random((nM, len(leafS))) < 0.5
Y  = rng.random((len(leafS), len(leafM))) < 0.5
get_s = lambda n, l: Y[sIdx[n], l] if n in sIdx else Zs[n, l]
get_m = lambda n, k: Y[k, mIdx[n]] if n in mIdx else Zm[n, k]

for _ in range(600):
    for l in range(len(leafM)):
        for node in [0] + intS:
            lp = np.log([1-alpha, alpha]) if node == 0 else np.log(Ms_true[int(get_s(parS[node], l)), :])
            for ch in (2*node+1, 2*node+2): lp = lp + np.log(Ms_true[:, int(get_s(ch, l))])
            Zs[node, l] = rng.random() < 1/(1+np.exp(lp[0]-lp[1]))
    for k in range(len(leafS)):
        for node in [0] + intM:
            lp = np.log([1-alpha, alpha]) if node == 0 else np.log(Mm_true[int(get_m(parM[node], k)), :])
            for ch in (2*node+1, 2*node+2): lp = lp + np.log(Mm_true[:, int(get_m(ch, k))])
            Zm[node, k] = rng.random() < 1/(1+np.exp(lp[0]-lp[1]))
    for k, sl in enumerate(leafS):
        for l, ml in enumerate(leafM):
            lp = np.log(Ms_true[int(get_s(parS[sl], l)), :]) + np.log(Mm_true[int(get_m(parM[ml], k)), :])
            Y[k, l] = rng.random() < 1/(1+np.exp(lp[0]-lp[1]))

kidsS = {}
for n in range(1, nS): kidsS.setdefault(parS[n], []).append(n)

def neg_prod(x):                       # what the C++ code maximises
    M = np.log(P(alpha, np.exp(x))); tot = 0.0
    for l in range(len(leafM)):
        tot += np.log(alpha if get_s(0, l) else 1-alpha)
        for n in range(1, nS): tot += M[int(get_s(parS[n], l)), int(get_s(n, l))]
    return -tot

def neg_pl(x):                         # normalised full conditionals (MPLE)
    M = np.log(P(alpha, np.exp(x))); tot = 0.0
    for l in range(len(leafM)):
        for n in range(nS):
            lp = np.log([1-alpha, alpha]) if n == 0 else M[int(get_s(parS[n], l)), :].copy()
            for c in kidsS.get(n, []): lp = lp + M[:, int(get_s(c, l))]
            if n in sIdx:              # leaf: add the other tree's factor
                lp = lp + np.log(Mm_true[int(get_m(parM[leafM[l]], sIdx[n])), :])
            tot += lp[int(get_s(n, l))] - np.logaddexp(lp[0], lp[1])
    return -tot

print(f"true log_nu                                  {TRUE:.4f}")
print(f"MLE, unnormalised product (current C++ code) {minimize_scalar(neg_prod, bounds=(-8,4), method='bounded').x:.4f}")
print(f"MPLE, normalised full conditionals           {minimize_scalar(neg_pl,   bounds=(-8,4), method='bounded').x:.4f}")
