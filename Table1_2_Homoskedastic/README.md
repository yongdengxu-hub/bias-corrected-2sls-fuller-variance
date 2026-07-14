# Table 1 (`Table_Est_b`) and Table 2 (`tab:size_main`): homoskedastic Monte Carlo

Run `main_table1_table2.m` (T=200, ρ=0.279, β=0.2, γ=0.6, CP={8,10,15,25,70}, L={0,2,4},
R=20,000, fixed seed `rng(20260628)`), then `extract_table1_table2.m` to print both LaTeX
table bodies from the same three saved `.mat` files — a single Monte Carlo run supplies both
the bias/variance statistics (Table 1) and the size/coverage statistics (Table 2).

The correction implemented in `doMonteCarlosSimulation_YX.m` (`TSLS_variance_BC`,
`fuller_variance_BC`) is the paper's common rule, used identically here and in Table 3
Panel A/L≥1 and Table 4: additive correction, reverting to the uncorrected variance
coordinate-by-coordinate wherever that would be negative.

**Reproducibility note**: at L=0, 2SLS has no finite population variance (Supplementary
Appendix S.4.1), so its sample Bias/Var/ratio entries in this regime are inherently seed-
dependent — verified directly this session (a fixed-seed rerun gave materially different
numbers at CP=15 than an earlier `rng('shuffle')` run: 0.097 vs. 0.307 for `Var`, both
"valid" since the population quantity being estimated does not exist). This run uses a fixed
seed for exact reproducibility of *this* run; a different seed will give different L=0
2SLS values while preserving the qualitative finding.
