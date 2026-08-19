"""Replicate the C++ simulation in numpy (exact Gibbs on p ∝ f_species · f_molecules, t=1,
true per-clique alpha/log_nu) and compare summary statistics + estimators against the
actual C++ output. Tells us whether the C++ simulated data really is the MRF at t=1."""
import numpy as np, pandas as pd
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

D = "/Users/visanim/work/random-markov-field/metabolite_inference/model_validation/s_balanced_255_m_balanced_255/"
rng = np.random.default_rng(3)

def load_tree(name):
    e = pd.read_csv(D + f"{name}.txt", sep="\t")
    order, seen, ch, pa = [], set(), set(), set()
    for c, p in zip(e["child"], e["parent"]):
        ch.add(c); pa.add(p)
        if p not in seen: order.append(p); seen.add(p)
        if c not in seen: order.append(c); seen.add(c)
    roots = pa - ch
    leaves = [n for n in order if n in ch and n not in pa]
    parent = dict(zip(e["child"], e["parent"]))
    idx = {n: i for i, n in enumerate(order)}
    par = np.array([idx[parent[n]] if n in parent else -1 for n in order])
    kids = [[] for _ in order]
    for c, p in parent.items(): kids[idx[p]].append(idx[c])
    return order, idx, par, kids, [idx[n] for n in leaves], [idx[n] for n in roots]

sOrd, sIx, sPar, sKids, sLeaf, sRoot = load_tree("species")
mOrd, mIx, mPar, mKids, mLeaf, mRoot = load_tree("molecules")
pv = pd.read_csv(D + "acol_input_simulated.txt", sep="\t"); pv = dict(zip(pv["name"], pv["value"].astype(float)))

nS, nM, L = len(sOrd), len(mOrd), len(sLeaf)
sLeafPos = {v: i for i, v in enumerate(sLeaf)}       # species node idx -> row in Y
mLeafPos = {v: i for i, v in enumerate(mLeaf)}

def mats(alphas, nus):
    return np.array([expm(np.array([[-a*n, a*n], [(1-a)*n, (a-1)*n]])) for a, n in zip(alphas, nus)])

sA = np.array([pv[f"species_alpha_{mOrd[c]}"]   for c in mLeaf])
sN = np.exp([pv[f"species_log_nu_{mOrd[c]}"]    for c in mLeaf])
mA = np.array([pv[f"molecules_alpha_{sOrd[c]}"] for c in sLeaf])
mN = np.exp([pv[f"molecules_log_nu_{sOrd[c]}"]  for c in sLeaf])
Ps, Pm = mats(sA, sN), mats(mA, mN)                  # (L,2,2), indexed by clique

S = rng.random((nS, L)) < 0.5     # species node x molecule-leaf clique
M = rng.random((nM, L)) < 0.5     # molecule node x species-leaf clique
def syncY_from_S():
    return np.array([S[v] for v in sLeaf])           # (L species-leaves, L cliques) == Y
Y = syncY_from_S()
for v in mLeaf: M[v] = Y[:, mLeafPos[v]]

lPs, lPm = np.log(Ps), np.log(Pm)
cl = np.arange(L)
for sweep in range(500):
    for n in range(nS):                              # species tree nodes (leaves handled via Y)
        if n in sLeafPos: continue
        lp = np.zeros((2, L))
        if n in sRoot:
            lp += np.log(np.stack([1-sA, sA]))
        else:
            lp += lPs[cl, S[sPar[n]].astype(int), :].T
        for c in sKids[n]:
            lp += lPs[cl, :, S[c].astype(int)].T
        S[n] = rng.random(L) < 1/(1+np.exp(lp[0]-lp[1]))
    for n in range(nM):
        if n in mLeafPos: continue
        lp = np.zeros((2, L))
        if n in mRoot:
            lp += np.log(np.stack([1-mA, mA]))
        else:
            lp += lPm[cl, M[mPar[n]].astype(int), :].T
        for c in mKids[n]:
            lp += lPm[cl, :, M[c].astype(int)].T
        M[n] = rng.random(L) < 1/(1+np.exp(lp[0]-lp[1]))
    # Y: all cells conditionally independent given S(internal) and M(internal)
    sp_par = np.array([S[sPar[v]] for v in sLeaf]).astype(int)          # (k, l)
    mo_par = np.array([M[mPar[v]] for v in mLeaf]).astype(int).T        # (k, l)
    kidx, lidx = np.arange(L)[:, None], np.arange(L)[None, :]
    lp0 = lPs[lidx, sp_par, 0] + lPm[kidx, mo_par, 0]
    lp1 = lPs[lidx, sp_par, 1] + lPm[kidx, mo_par, 1]
    Y = rng.random((L, L)) < 1/(1+np.exp(lp0-lp1))
    for i, v in enumerate(sLeaf): S[v] = Y[i]
    for j, v in enumerate(mLeaf): M[v] = Y[:, j]

def agreement(S):
    ai = at = bi = bt = 0
    for n in range(nS):
        if n in sRoot: continue
        same = (S[n] == S[sPar[n]])
        if n in sLeafPos: ai += same.sum(); at += L
        else:            bi += same.sum(); bt += L
    return bi/bt, ai/at

def fit(skip_leaves=False):
    out = []
    for c in range(40):
        a = sA[c]
        def nll(x):
            lm = np.log(expm(np.array([[-a*np.exp(x), a*np.exp(x)], [(1-a)*np.exp(x), (a-1)*np.exp(x)]])))
            t = np.log(a if S[sRoot[0], c] else 1-a)
            for n in range(nS):
                if n in sRoot: continue
                if skip_leaves and n in sLeafPos: continue
                t += lm[int(S[sPar[n], c]), int(S[n, c])]
            return -t
        out.append(minimize_scalar(nll, bounds=(-8, 4), method="bounded").x)
    return np.mean(out)

gi, gl = agreement(S)
print("            replicated-in-numpy MRF (t=1)   |   actual C++ simulated data")
print(f"internal agreement   {gi:.4f}                    |   0.9578")
print(f"leaf agreement       {gl:.4f}                    |   0.9788")
print(f"MLE (code objective) {fit():.4f}                   |   -2.5424")
print(f"MLE (skip leaves)    {fit(True):.4f}                   |   -2.1917")
print(f"true mean log_nu     {np.log(sN[:40]).mean():.4f}")

# --- MPLE (normalised full conditionals) on the replicated data ---
sKidsArr = sKids
def mple():
    out = []
    for c in range(40):
        a = sA[c]
        mo_par_c = np.array([int(M[mPar[mLeaf[c]], k]) for k in range(L)])   # molecule parent state
        def nll(x):
            lm = np.log(expm(np.array([[-a*np.exp(x), a*np.exp(x)], [(1-a)*np.exp(x), (a-1)*np.exp(x)]])))
            tot = 0.0
            for n in range(nS):
                lp = np.log([1-a, a]) if n in sRoot else lm[int(S[sPar[n], c]), :].copy()
                for ch in sKidsArr[n]: lp = lp + lm[:, int(S[ch, c])]
                if n in sLeafPos:
                    k = sLeafPos[n]
                    lp = lp + lPm[k, mo_par_c[k], :]
                tot += lp[int(S[n, c])] - np.logaddexp(lp[0], lp[1])
            return -tot
        out.append(minimize_scalar(nll, bounds=(-8, 4), method="bounded").x)
    return np.mean(out)
print(f"MPLE (normalised)    {mple():.4f}   <- on replicated data, truth {np.log(sN[:40]).mean():.4f}")

np.savez("/tmp/replicated_states.npz", S=S, M=M, sA=sA, sN=sN, mA=mA, mN=mN,
         sPar=sPar, mPar=mPar, sLeaf=np.array(sLeaf), mLeaf=np.array(mLeaf), sRoot=np.array(sRoot))
print("saved states")
