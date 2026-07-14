# Table 4, Panel A: Acemoglu, Johnson, and Robinson (2001)

Run `AJR_reverttable.m` directly (loads `Table7_Data.xlsx`, no simulation, deterministic).
Prints the LaTeX rows for both 2SLS and Fuller, β and γ, specifications (1)–(9), matching
Table 4 Panel A exactly, plus a `revert` diagnostic row showing which cells reverted.

Uses the paper's homoskedastic bias correction (Sections S.1–S.2) applied directly to the
real data, with the same per-coordinate revert rule as every other non-Fuller-heteroskedastic
case in the paper (`Var_BC = Var - bias`, reverting to the uncorrected `Var` where negative).

**Note on specifications (1) and (5), Fuller β**: the estimated bias there is large but still
positive (not negative), so revert never fires and the correction is applied in full,
producing se_BC = 0.095 (spec 1) and 0.036 (spec 5) — the correction removes 81% and 97% of
the estimated variance respectively. This is a deliberate, documented feature of using the
same rule everywhere (see the main paper's Section 6 discussion and Supplementary Appendix
S.4.2.3), not an error.
