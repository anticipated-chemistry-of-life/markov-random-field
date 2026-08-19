"""Does replacing the unnormalised product by the *normalised full conditionals*
(maximum pseudo-likelihood, the standard consistent estimator for an MRF)
recover the true log_nu on the real simulated data?"""
import numpy as np, pandas as pd
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

D = "/Users/visanim/work/random-markov-field/metabolite_inference/model_validation/s_balanced_255_m_balanced_255/"

def P(alpha, nu, t=1.0):
    L = np.array([[-alpha*nu, alpha*nu], [(1-alpha)*nu, (alpha-1)*nu]])
    return expm(L*t)

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
    kids = {}
    for c, p in parent.items(): kids.setdefault(p, []).append(c)
    return order, roots, leaves, parent, kids

order, roots, leaves, parS, kidsS = load_tree("species")
morder, mroots, mleaves, parM, kidsM = load_tree("molecules")
pv = pd.read_csv(D + "acol_input_simulated.txt", sep="\t")
pv = dict(zip(pv["name"], pv["value"].astype(float)))

Y  = pd.read_csv(D + "acol_simulated_Y.txt", sep="\t")
Zs = pd.read_csv(D + "acol_simulated_Z_species.txt", sep="\t")
Zm = pd.read_csv(D + "acol_simulated_Z_molecules.txt", sep="\t")
st = {}
for s, m, v in zip(Y["species"], Y["molecules"], Y["Y_state"]): st[(s, m)] = int(v)
for s, m, v in zip(Zs["species"], Zs["molecules"], Zs["Z_state"]): st[(s, m)] = int(v)
mst = dict(st)
for s, m, v in zip(Zm["species"], Zm["molecules"], Zm["Z_state"]): mst[(s, m)] = int(v)

rows = []
for clique in mleaves[:40]:
    alpha = pv[f"species_alpha_{clique}"]
    true_ln = pv[f"species_log_nu_{clique}"]
    # fixed molecule-tree factor entering each species leaf of this clique
    mol = {}
    for sl in leaves:
        Mm = P(pv[f"molecules_alpha_{sl}"], np.exp(pv[f"molecules_log_nu_{sl}"]))
        pm = mst[(sl, parM[clique])]
        mol[sl] = np.log(Mm[pm, :])

    def neg_pl(x):
        M = np.log(P(alpha, np.exp(x)))
        tot = 0.0
        for n in order:
            lp = np.zeros(2)
            if n in roots: lp += np.log([1-alpha, alpha])
            else:          lp += M[st[(parS[n], clique)], :]
            for c in kidsS.get(n, []): lp += M[:, st[(c, clique)]]
            if n in leaves: lp += mol[n]              # the other tree's factor
            tot += lp[st[(n, clique)]] - np.logaddexp(lp[0], lp[1])
        return -tot

    rows.append((true_ln, minimize_scalar(neg_pl, bounds=(-8, 4), method="bounded").x))

df = pd.DataFrame(rows, columns=["true", "mple"])
print(f"mean true log_nu                                   {df.true.mean():.4f}")
print(f"mean MPLE (normalised full conditionals)           {df.mple.mean():.4f}   bias {(df.mple-df.true).mean():+.4f}")
print(f"mean MLE  (current code: unnormalised product)     -2.5424   bias -1.0249")
print(f"mean MLE  (--skip_leaves_in_alpha_nu_update)       -2.1917   bias -0.6742")
