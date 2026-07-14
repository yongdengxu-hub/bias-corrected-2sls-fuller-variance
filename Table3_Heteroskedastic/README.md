# Table 3 (`Table_Est_HC`): heteroskedasticity-robust (HC0) Monte Carlo

Run `main_HC_sizeci_delta78.m` (T=200, ρ=0.279, δ≈0.78, R=20,000, fixed seed `rng(20260628)`),
then `extract_table_est_hc78.m` to print the exact LaTeX table body (Bias, Var,
$\hat V_{HC}/Var$, $\hat V_{HC}^{BC}/Var$ for β and γ, both panels, L=0,2,4).

The correction implemented in `doMonteCarlosSimulation_YX_Hetro.m`
(`fuller_variance_BC_HC`) is the two-branch rule described in Supplementary Appendix S.4.2.3:

- **L=0**: HC0 sandwich + Theorem-4 bias, gated multiplicatively (`cpr:posBC`, window
  (0,10) in $\widehat{CP}^*$). This is the only case in the whole paper that uses the gate —
  see S.4.2.3 for why (the Theorem-4 bias estimate has a genuinely wild tail at L=0, and
  every non-gated alternative tested underperforms it there).
- **L≥1**: conventional $k$-class variance + its own bias, applied additively with a simple
  per-coordinate revert-to-uncorrected rule if negative — the same rule used everywhere else
  in the paper (Table 1/2, Table 3 Panel A, Table 4).

2SLS (`TSLS_variance_BC_HC`) uses the additive+revert rule at every L; it has no finite
population variance at L=0 (Supplementary Appendix S.9), so its L=0 entries are large and not
exactly reproducible across seeds (see the top-level README's reproducibility notes).
