"""Small two-tree MRF, simulated exactly the way the C++ code does, then estimated
the way the C++ code does. Sweep the *other* tree's nu to show the bias is caused by
the product-of-two-trees coupling (i.e. the missing normalising constant)."""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

rng = np.random.default_rng(1)

def P(alpha, nu, t=1.0):
    L = np.array([[-alpha*nu, alpha*nu], [(1-alpha)*nu, (alpha-1)*nu]])
    return expm(L*t)

def balanced(depth):
    """returns (n_nodes, parent[], leaves[], internal_nonroot[], root=0)"""
    n = 2**(depth+1) - 1
    parent = [-1] + [(i-1)//2 for i in range(1, n)]
    leaves = [i for i in range(n) if 2*i+1 >= n]
    internal = [i for i in range(1, n) if 2*i+1 < n]
    return n, parent, leaves, internal

DEPTH = 5
nS, parS, leafS, intS = balanced(DEPTH)
nM, parM, leafM, intM = balanced(DEPTH)
alpha = 0.5
true_log_nu_s = -1.5

def simulate(log_nu_m, sweeps=400):
    Ms = P(alpha, np.exp(true_log_nu_s))
    Mm = P(alpha, np.exp(log_nu_m))
    # Zs[i, l] : species node i, molecule leaf l   (internal species nodes + root)
    Zs = rng.random((nS, len(leafM))) < 0.5
    Zm = rng.random((nM, len(leafS))) < 0.5
    Y  = rng.random((len(leafS), len(leafM))) < 0.5
    sIdx = {n: k for k, n in enumerate(leafS)}
    mIdx = {n: k for k, n in enumerate(leafM)}

    def get_s(node, l):   # state of species-tree node `node` in clique `l` (molecule leaf idx)
        return Y[sIdx[node], l] if node in sIdx else Zs[node, l]
    def get_m(node, k):
        return Y[k, mIdx[node]] if node in mIdx else Zm[node, k]

    for _ in range(sweeps):
        # Gibbs on species internal nodes (parent term + children terms, species tree only)
        for l in range(len(leafM)):
            for node in [0] + intS:
                lp = np.zeros(2)
                if node == 0:
                    lp += np.log([1-alpha, alpha])
                else:
                    ps = int(get_s(parS[node], l)); lp += np.log(Ms[ps, :])
                for ch in (2*node+1, 2*node+2):
                    cs = int(get_s(ch, l)); lp += np.log(Ms[:, cs])
                Zs[node, l] = rng.random() < 1/(1+np.exp(lp[0]-lp[1]))
        for k in range(len(leafS)):
            for node in [0] + intM:
                lp = np.zeros(2)
                if node == 0:
                    lp += np.log([1-alpha, alpha])
                else:
                    pm = int(get_m(parM[node], k)); lp += np.log(Mm[pm, :])
                for ch in (2*node+1, 2*node+2):
                    cs = int(get_m(ch, k)); lp += np.log(Mm[:, cs])
                Zm[node, k] = rng.random() < 1/(1+np.exp(lp[0]-lp[1]))
        # Gibbs on Y: product of experts over the two trees
        for k, sl in enumerate(leafS):
            for l, ml in enumerate(leafM):
                ps = int(get_s(parS[sl], l)); pm = int(get_m(parM[ml], k))
                lp = np.log(Ms[ps, :]) + np.log(Mm[pm, :])
                Y[k, l] = rng.random() < 1/(1+np.exp(lp[0]-lp[1]))
    return Zs, Y, sIdx

def mle(Zs, Y, sIdx, skip_leaves=False):
    def nll(x):
        M = P(alpha, np.exp(x)); s = 0.0
        for l in range(Y.shape[1]):
            g = lambda n: Y[sIdx[n], l] if n in sIdx else Zs[n, l]
            s += np.log(alpha if g(0) else 1-alpha)
            for node in range(1, nS):
                if skip_leaves and node in sIdx: continue
                s += np.log(M[int(g(parS[node])), int(g(node))])
        return -s
    return minimize_scalar(nll, bounds=(-8, 4), method="bounded").x

print(f"species tree: {nS} nodes, molecule tree: {nM} nodes, true species log_nu = {true_log_nu_s}")
print(f"{'other tree log_nu':>18} {'MLE all nodes':>14} {'MLE skip leaves':>16}")
for log_nu_m in [-1.5, 0.0, 1.5, 3.0, 5.0]:
    Zs, Y, sIdx = simulate(log_nu_m)
    print(f"{log_nu_m:>18.1f} {mle(Zs, Y, sIdx):>14.3f} {mle(Zs, Y, sIdx, True):>16.3f}")
print("\n(large other-tree log_nu => that tree is ~independent => coupling vanishes)")
