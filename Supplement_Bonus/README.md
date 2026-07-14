# Supplement scripts (bonus — lower verification bar)

These are provided for completeness, "in case anyone is interested" in the Supplement's
tables. Unlike the four main-text tables (see the top-level README), the scripts below have
**not** all been re-verified this session to reproduce the currently-published Supplement
numbers exactly. Confidence level is noted per folder.

## `S8_S10_HighEndog_DeltaRobustness/` — HIGH confidence, used and verified this session

- `main_HC_sizeci_delta78_rho08.m` + `extract_het78_rho08.m` → Supplement Tables `tab:het_var08`,
  `tab:het_sizeci08` (S.8.2, high endogeneity ρ=0.8).
- `extract_het78_sizeci.m` → Supplement Table `tab:het_sizeci` (S.8.1); shares the driver
  `main_HC_sizeci_delta78.m` already in `../Table3_Heteroskedastic/` — run that first.
- `main_HC_sizeci_delta95.m` + `extract_table_est_hc95.m`, `main_HC_sizeci_delta40.m` +
  `extract_table_est_hc40.m` → Supplement Tables `tab:het_delta95`, `tab:het_delta40` (S.10,
  δ-robustness). All four also need `doMonteCarlosSimulation_YX_Hetro.m`, `genY_Hetro.m`,
  `FirstStageTest.m`, `getApproximation.m`, `getstat.m`, `conv_power.m`, `X2.mat` from
  `../Table3_Heteroskedastic/`.

## `S10_NonNormalErrors/` — MEDIUM confidence, located but not rerun this session

`run_CI_nn.m` (with its own docstring identifying it precisely: fat-tailed $t_5$ and
right-skewed $\chi^2_4$ errors, ρ∈{0.279,0.8}, fixed seed) and `make_tables_nn.m` (table
body builder) appear to be the scripts behind Supplement Tables `tab:CI_t5`, `tab:CI_chi2`
(S.10, non-normal-error robustness). Located via keyword search this session; not rerun to
confirm exact reproduction of the published numbers.

## `S6_S7_VarRatio_Size_HighEndog_L20/` — LOW confidence, requires manual code editing

`run_homo_all.m` and `run_size_rho.m` appear to target Supplement Tables `tab:VR_rho08`,
`tab:VR_L20`, `tab:CI_L20` (S.6) and `tab:size_rho08`, `tab:CI_rho0279`, `tab:CI_rho08` (S.7).
**Important**: `run_homo_all.m`'s own header comment states it "REQUIRES
`doMonteCarlosSimulation_YX_Hetro.m` switch in HOMO mode" — i.e., the shared heteroskedastic
driver's dispatch (currently wired to call the heteroskedasticity-robust `TSLS_variance_BC_HC`
/ `fuller_variance_BC_HC`) would need to be manually edited back to call the plain
homoskedastic `TSLS_variance_BC` / `fuller_variance_BC` local functions (which still exist in
the file, just unused by the current dispatch) before this script will produce the intended
homoskedastic output. This was not attempted or verified this session.

## Not included

- **S.5 (`tab:diagnostic`)**: a closed-form percentage-reduction formula, not a simulation —
  no script needed; see the Supplement text for the formula.
- **S.3.1 (`tab:contraction`)**: an analytical derivation table, not a simulation.
