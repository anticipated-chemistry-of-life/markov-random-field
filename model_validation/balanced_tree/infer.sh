../../build/acol infer \
    --out ./branch_given_but_not_fixed/acol \
    --tree_species species.txt \
    --tree_molecules molecules.txt \
    --lotus acol_simulated_lotus.tsv \
    --numThreads 1 \
    --write_Y_trace \
    --molecules_branch_lengths acol_input_simulated.txt \
    --species_branch_lengths acol_input_simulated.txt \
    --molecules_mu_1 acol_input_simulated.txt \
    --molecules_mu_0 acol_input_simulated.txt \
    --species_mu_1 acol_input_simulated.txt \
    --species_mu_0 acol_input_simulated.txt \
    --writeTrace \
    --writeBurnin \
    --Z.update false \
    --set_species_Z acol_simulated_Z_species.txt \
    --set_molecules_Z acol_simulated_Z_molecules.txt \
    --species_mu_1.update false \
    --molecules_mu_0.update false \
    --molecules_mu_1.update false \
    --molecules_branch_lengths.update false \
    --species_branch_lengths.update false \
    --species_mu_0.update false \
    # --Y.update false \
    # --set_Y acol_simulated_Y.txt \