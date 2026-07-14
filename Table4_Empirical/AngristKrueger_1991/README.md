# Table 4, Panel B: Angrist and Krueger (1991)

Run `replicate_AK1991.m` directly (loads `AK1991_3039.xlsx`, deterministic, no simulation).
Prints 2SLS and Fuller estimates for β (Educ) and γ (YOB30), specifications I–III.

**Reproducibility caveat — specification III**: the first-stage design matrix in this
specification (30 quarter-of-birth × year-of-birth interaction instruments, plus age and
age² controls that are nearly collinear with the year-of-birth dummies) is numerically close
to singular (condition number on the order of 1e-26 at double precision, confirmed by direct
check this session). Because of this, spec III's reported figures — including the first-stage
$\hat F$ — can shift in the second or third decimal place across MATLAB versions, LAPACK/BLAS
backends, or hardware, even though the calculation itself is fully deterministic given fixed
data. This is a known feature of this specification (it is included in the literature
specifically as an illustration of a near-degenerate weak-instrument design), not a bug in
this code. Specifications I and II do not have this issue.

Fuller's bias is negative at L=29 (specs II, III) for both coefficients, so the additive
correction only ever *increases* the reported variance there and the revert rule never fires;
2SLS is likewise well-behaved (no reverts triggered in any specification).
