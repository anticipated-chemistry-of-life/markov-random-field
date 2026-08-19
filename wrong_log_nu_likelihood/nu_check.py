import numpy as np, pandas as pd
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
import sys

D = "/Users/visanim/work/random-markov-field/metabolite_inference/model_validation/s_balanced_255_m_balanced_255/"

tree_name = "species"
other = "molecules"

edges = pd.read_csv(D + f"{tree_name}.txt", sep="\t")
# replicate C++ node ordering
node_order, node_set, is_child, is_parent = [], set(), set(), set()
for c, p in zip(edges["child"], edges["parent"]):
    is_child.add(c); is_parent.add(p)
    if p not in node_set: node_order.append(p); node_set.add(p)
    if c not in node_set: node_order.append(c); node_set.add(c)
roots = is_parent - is_child
leaves = [n for n in node_order if n in is_child and n not in is_parent]
parent_of = dict(zip(edges["child"], edges["parent"]))
print("n_nodes", len(node_order), "n_roots", len(roots), "n_leaves", len(leaves))

params = pd.read_csv(D + "acol_input_simulated.txt", sep="\t")
pv = dict(zip(params["name"], params["value"].astype(float)))

Y = pd.read_csv(D + "acol_simulated_Y.txt", sep="\t")
Zs = pd.read_csv(D + f"acol_simulated_Z_{tree_name}.txt", sep="\t")

# state lookup: (node_name_in_this_tree, other_leaf_name) -> state
state = {}
for s, m, v in zip(Y[tree_name], Y[other], Y["Y_state"]):
    state[(s, m)] = bool(v)
for s, m, v in zip(Zs[tree_name], Zs[other], Zs["Z_state"]):
    state[(s, m)] = bool(v)

n_bins = 100
delta = 2.0 / (n_bins + 1.0)
t_eff = delta * (n_bins / 2 + 0.5)
print("t_eff", t_eff)

def P(alpha, nu, t):
    L = np.array([[-alpha * nu, alpha * nu], [(1 - alpha) * nu, (alpha - 1) * nu]])
    return expm(L * t)

# cliques of the species tree are indexed by molecule leaves
mol_edges = pd.read_csv(D + f"{other}.txt", sep="\t")
mo, ms_, mis_c, mis_p = [], set(), set(), set()
for c, p in zip(mol_edges["child"], mol_edges["parent"]):
    mis_c.add(c); mis_p.add(p)
    if p not in ms_: mo.append(p); ms_.add(p)
    if c not in ms_: mo.append(c); ms_.add(c)
mol_leaves = [n for n in mo if n in mis_c and n not in mis_p]
print("n cliques (molecule leaves)", len(mol_leaves))

def LL(log_nu, clique, alpha, skip_leaves):
    nu = np.exp(log_nu)
    M = P(alpha, nu, t_eff)
    tot = 0.0
    for n in node_order:
        if skip_leaves and n in leaves:
            continue
        st = state[(n, clique)]
        if n in roots:
            tot += np.log(alpha if st else 1 - alpha)
        else:
            ps = state[(parent_of[n], clique)]
            tot += np.log(M[int(ps), int(st)])
    return tot

rows = []
for clique in mol_leaves[:40]:
    alpha = pv[f"{tree_name}_alpha_{clique}"]
    true_log_nu = pv[f"{tree_name}_log_nu_{clique}"]
    res_all = minimize_scalar(lambda x: -LL(x, clique, alpha, False), bounds=(-8, 4), method="bounded")
    res_skip = minimize_scalar(lambda x: -LL(x, clique, alpha, True), bounds=(-8, 4), method="bounded")
    rows.append((clique, true_log_nu, res_all.x, res_skip.x))

df = pd.DataFrame(rows, columns=["clique", "true_log_nu", "mle_all_nodes", "mle_skip_leaves"])
print(df.head(15).to_string())
print()
print("mean true      ", df.true_log_nu.mean())
print("mean MLE all   ", df.mle_all_nodes.mean())
print("mean MLE skip  ", df.mle_skip_leaves.mean())
