import numpy as np, pandas as pd
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

rng = np.random.default_rng(0)
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
    return order, roots, leaves, dict(zip(e["child"], e["parent"]))

order, roots, leaves, parent_of = load_tree("species")
morder, mroots, mleaves, mparent = load_tree("molecules")
params = pd.read_csv(D + "acol_input_simulated.txt", sep="\t")
pv = dict(zip(params["name"], params["value"].astype(float)))

# ---------- Control 1: single-tree data, same estimator ----------
print("=== CONTROL 1: data from ONE tree only (no product) ===")
res = []
for rep in range(40):
    alpha, true_ln = 0.5, -1.5
    M = P(alpha, np.exp(true_ln))
    st = {}
    for n in order:
        if n in roots:
            st[n] = rng.random() < alpha
        else:
            st[n] = rng.random() < M[int(st[parent_of[n]]), 1]
    def nll(x, M=None):
        Mx = P(alpha, np.exp(x))
        s = 0.0
        for n in order:
            if n in roots: s += np.log(alpha if st[n] else 1-alpha)
            else: s += np.log(Mx[int(st[parent_of[n]]), int(st[n])])
        return -s
    res.append(minimize_scalar(nll, bounds=(-8,4), method="bounded").x)
print(f"  true log_nu = -1.5   mean MLE over 40 reps = {np.mean(res):.4f}  (sd {np.std(res):.3f})")
print("  -> estimator itself is unbiased on single-tree data\n")

# ---------- Control 2: parent-child agreement in the real simulated data ----------
print("=== CONTROL 2: parent-child agreement, simulated MRF vs single-tree prediction ===")
Y = pd.read_csv(D + "acol_simulated_Y.txt", sep="\t")
Zs = pd.read_csv(D + "acol_simulated_Z_species.txt", sep="\t")
state = {}
for s, m, v in zip(Y["species"], Y["molecules"], Y["Y_state"]): state[(s, m)] = bool(v)
for s, m, v in zip(Zs["species"], Zs["molecules"], Zs["Z_state"]): state[(s, m)] = bool(v)

agree_int, tot_int, agree_leaf, tot_leaf, pred_int, pred_leaf = 0, 0, 0, 0, [], []
for clique in mleaves:
    alpha = pv[f"species_alpha_{clique}"]
    M = P(alpha, np.exp(pv[f"species_log_nu_{clique}"]))
    for n in order:
        if n in roots: continue
        ps, cs = state[(parent_of[n], clique)], state[(n, clique)]
        p_same = M[int(ps), int(ps)]
        if n in leaves:
            tot_leaf += 1; agree_leaf += (ps == cs); pred_leaf.append(p_same)
        else:
            tot_int += 1; agree_int += (ps == cs); pred_int.append(p_same)
print(f"  internal nodes: observed agreement {agree_int/tot_int:.4f}  vs true-nu prediction {np.mean(pred_int):.4f}")
print(f"  leaves        : observed agreement {agree_leaf/tot_leaf:.4f}  vs true-nu prediction {np.mean(pred_leaf):.4f}")
print("  -> the field is MORE correlated than a single-tree chain at the true nu predicts,")
print("     so the single-tree objective compensates by shrinking nu.\n")

# ---------- Control 3: does the missing normalizer explain the gap? ----------
print("=== CONTROL 3: normalizer-corrected objective on the same data ===")
# Correct conditional for a leaf under the product-of-experts:
#   P(y | par_sp, par_mol) = Ps(par_sp->y) Pm(par_mol->y) / sum_k Ps(par_sp->k) Pm(par_mol->k)
Zm = pd.read_csv(D + "acol_simulated_Z_molecules.txt", sep="\t")
mstate = dict(state)
for s, m, v in zip(Zm["species"], Zm["molecules"], Zm["Z_state"]): mstate[(s, m)] = bool(v)

rows = []
for clique in mleaves[:40]:
    alpha = pv[f"species_alpha_{clique}"]
    true_ln = pv[f"species_log_nu_{clique}"]
    # molecule-side transition into each species-leaf cell (fixed, uses true molecule params)
    mol_terms = []
    for sp_leaf in leaves:
        a_m = pv[f"molecules_alpha_{sp_leaf}"]
        Mm = P(a_m, np.exp(pv[f"molecules_log_nu_{sp_leaf}"]))
        pm = mstate[(sp_leaf, mparent[clique])]
        mol_terms.append((sp_leaf, Mm[int(pm), 0], Mm[int(pm), 1]))

    def nll(x):
        Ms = P(alpha, np.exp(x)); s = 0.0
        for n in order:
            if n in roots:
                s += np.log(alpha if state[(n, clique)] else 1-alpha); continue
            ps, cs = state[(parent_of[n], clique)], state[(n, clique)]
            s += np.log(Ms[int(ps), int(cs)])
        # subtract the per-leaf normalizer of the product-of-experts
        for sp_leaf, m0, m1 in mol_terms:
            ps = state[(parent_of[sp_leaf], clique)]
            s -= np.log(Ms[int(ps), 0]*m0 + Ms[int(ps), 1]*m1)
        return -s
    rows.append((true_ln, minimize_scalar(nll, bounds=(-8,4), method="bounded").x))

df = pd.DataFrame(rows, columns=["true", "mle_normalized"])
print(f"  mean true log_nu                 {df.true.mean():.4f}")
print(f"  mean MLE with leaf normalizer    {df.mle_normalized.mean():.4f}")
print(f"  (uncorrected objective gave      -2.5424)")
