../../build/acol infer \
    --out ./test_out/acol \
    --tree_species species.txt \
    --tree_molecules molecules.txt \
    --lotus acol_simulated_lotus.tsv \
    --iterations 20000 \
    --numThreads 1 \
    --writeTrace \
    --writeBurnin \
    --write_Y_trace \
    --set_Y acol_simulated_Y.txt \
    --Y.update false \
    --write_Z_trace \
    --molecules_branch_lengths acol_input_simulated.txt \
    --molecules_branch_lengths.update false \
    --molecules_log_nu acol_input_simulated.txt \
    --molecules_log_nu.update false \
    --molecules_mean_log_nu acol_input_simulated.txt \
    --molecules_mean_log_nu.update false \
    --molecules_var_log_nu acol_input_simulated.txt \
    --molecules_var_log_nu.update false \
    --molecules_alpha acol_input_simulated.txt \
    --molecules_alpha.update false \
    --species_var_log_nu acol_input_simulated.txt \
    --species_var_log_nu.update false \
    --species_mean_log_nu acol_input_simulated.txt \
    --species_mean_log_nu.update false \
    --species_log_nu acol_simulated.txt \
    --species_log_nu.update false \
    --species_alpha acol_input_simulated.txt \
    --species_alpha.update false \
    --species_branch_lengths acol_species_simulated.txt \
    --species_branch_lengths.update false \
    # --set_species_Z acol_simulated_Z_species.txt \
    # --set_molecules_Z acol_simulated_Z_molecules.txt \
    # --Z.update false \