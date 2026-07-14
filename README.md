# Replication package: "Almost Unbiased Variance Estimation in Instrumental Variable Regression"

MATLAB R2024b/R2026a. Each subfolder is self-contained: run the single named driver script
from within that folder (adds no external path dependencies beyond the files included).

## Main-paper tables

| Table | Folder | Driver script | Notes |
|---|---|---|---|
| Table 1 (`Table_Est_b`) | `Table1_2_Homoskedastic/` | `main_table1_table2.m` | homoskedastic 2SLS/Fuller bias & variance |
| Table 2 (`tab:size_main`) | `Table1_2_Homoskedastic/` | (same run) | test size & CI coverage, from the same simulation |
| Table 3 (`Table_Est_HC`) | `Table3_Heteroskedastic/` | `main_HC_sizeci_delta78.m` | heteroskedasticity-robust (HC0) bias & variance |
| Table 4, Panel A | `Table4_Empirical/AJR_2001/` | `AJR_reverttable.m` | Acemoglu–Johnson–Robinson (2001) |
| Table 4, Panel B | `Table4_Empirical/AngristKrueger_1991/` | `replicate_AK1991.m` | Angrist–Krueger (1991) |

Each driver saves its own output (`.mat` for the two simulations, `.txt` for the two
empirical applications); each folder's `extract_*.m` (where present) turns the saved output
into the exact LaTeX table rows used in the paper.

### Reproducibility notes

- **Table 1/2, 2SLS at L=0**: the 2SLS estimator has no finite population variance at
  L=0 (Section S.4.1 of the Supplement). Its sample variance in this regime is therefore
  not exactly reproducible across simulation seeds — rerunning with a different seed will
  give a different, equally valid, extreme value. The qualitative finding (order-of-magnitude
  instability driven by non-existent moments) is robust; the specific digits are not.
- **Table 4, Panel B, specification III** (Angrist–Krueger, 29 instruments, age/age² controls
  nearly collinear with year-of-birth dummies): the first-stage design matrix is close to
  numerically singular (condition number ~1e-26 at double precision). Its reported figures
  can shift in the third decimal place across MATLAB versions or hardware; the qualitative
  finding (weak instruments, bias correction reverses significance) is robust.
- **Table 1/2 vs. Table 3's `L>=1` Fuller correction**: as of the current draft, the Fuller
  heteroskedastic correction uses a gated construction only at `L=0`; at `L>=1` it uses the
  same additive-correction-with-revert rule as every other case in the paper (Table 1/2,
  Table 3 Panel A, Table 4). See Supplementary Appendix S.4.2.3 for the full derivation and
  the evidence behind this choice.
- All simulation drivers use a fixed seed (`rng(20260628)`) for exact reproducibility, except
  where noted above.

## Supplement tables

Scripts and data for the Supplement's tables (S.5–S.10) are in `Supplement_Bonus/`, provided
for completeness. These have not been re-verified to the same standard as the four main-text
tables above; see that folder's own notes for what is and is not confirmed to reproduce the
currently published Supplement numbers exactly.
