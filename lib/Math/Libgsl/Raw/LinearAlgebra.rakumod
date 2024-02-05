use v6;

unit module Math::Libgsl::Raw::LinearAlgebra:ver<0.0.4>:auth<zef:FRITH>;

use NativeCall;
use Math::Libgsl::Raw::Complex :ALL;
use Math::Libgsl::Raw::Matrix :ALL;
use Math::Libgsl::Raw::Matrix::Complex64 :ALL;
use Math::Libgsl::Raw::Permutation :ALL;

constant GSLHELPER = %?RESOURCES<libraries/gslhelper>;

sub LIB {
  run('/sbin/ldconfig', '-p', :chomp, :out)
    .out
    .slurp(:close)
    .split("\n")
    .grep(/^ \s+ libgsl\.so\. \d+ /)
    .sort
    .head
    .comb(/\S+/)
    .head;
}

# LU Decomposition
sub gsl_linalg_LU_decomp(gsl_matrix $A, gsl_permutation $p, int32 $signum is rw --> int32) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_complex_LU_decomp(gsl_matrix_complex $A, gsl_permutation $p, int32 $signum is rw --> int32) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_LU_solve(gsl_matrix $LU, gsl_permutation $p, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_complex_LU_solve(gsl_matrix_complex $A, gsl_permutation $p, gsl_vector_complex $b, gsl_vector_complex $x --> int32) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_LU_svx(gsl_matrix $LU, gsl_permutation $p, gsl_vector $x --> int32) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_complex_LU_svx(gsl_matrix_complex $LU, gsl_permutation $p, gsl_vector_complex $x --> int32) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_LU_refine(gsl_matrix $A, gsl_matrix $LU, gsl_permutation $p, gsl_vector $b, gsl_vector $x, gsl_vector $work  --> int32) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_complex_LU_refine(gsl_matrix_complex $A, gsl_matrix_complex $LU, gsl_permutation $p, gsl_vector_complex $b, gsl_vector_complex $x, gsl_vector_complex $work --> int32) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_LU_invert(gsl_matrix $LU, gsl_permutation $p, gsl_matrix $inverse --> int32) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_complex_LU_invert(gsl_matrix_complex $LU, gsl_permutation $p, gsl_matrix_complex $inverse --> int32) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_LU_det(gsl_matrix $LU, int32 $signum --> num64) is native(&LIB) is export(:LU) { * }
sub mgsl_linalg_complex_LU_det(gsl_matrix_complex $LU, int32 $signum, gsl_complex $res) is native(GSLHELPER) is export(:LU) { * }
sub gsl_linalg_LU_lndet(gsl_matrix $LU --> num64) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_complex_LU_lndet(gsl_matrix_complex $LU --> num64) is native(&LIB) is export(:LU) { * }
sub gsl_linalg_LU_sgndet(gsl_matrix $LU, int32 $signum --> int32) is native(&LIB) is export(:LU) { * }
sub mgsl_linalg_complex_LU_sgndet(gsl_matrix_complex $LU, int32 $signum, gsl_complex $res) is native(GSLHELPER) is export(:LU) { * }
# QR Decomposition
sub gsl_linalg_QR_decomp(gsl_matrix $A, gsl_vector $tau --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_solve(gsl_matrix $QR, gsl_vector $tau, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_svx(gsl_matrix $QR, gsl_vector $tau, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_lssolve(gsl_matrix $QR, gsl_vector $tau, gsl_vector $b, gsl_vector $x, gsl_vector $residual --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_QTvec(gsl_matrix $QR, gsl_vector $tau, gsl_vector $v --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_Qvec(gsl_matrix $QR, gsl_vector $tau, gsl_vector $v --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_QTmat(gsl_matrix $QR, gsl_vector $tau, gsl_matrix $B --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_Rsolve(gsl_matrix $QR, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_Rsvx(gsl_matrix $QR, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_unpack(gsl_matrix $QR, gsl_vector $tau, gsl_matrix $Q, gsl_matrix $R --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_QRsolve(gsl_matrix $Q, gsl_matrix $R, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QR_update(gsl_matrix $Q, gsl_matrix $R, gsl_vector $w, gsl_vector $v --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_R_solve(gsl_matrix $R, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_R_svx(gsl_matrix $R, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
# QR Decomposition with Column Pivoting
sub gsl_linalg_QRPT_decomp(gsl_matrix $A, gsl_vector $tau, gsl_permutation $p, int32 $signum is rw, gsl_vector $norm --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_decomp2(gsl_matrix $A, gsl_matrix $q, gsl_matrix $r, gsl_vector $tau, gsl_permutation $p, int32 $signum is rw, gsl_vector $norm --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_solve(gsl_matrix $QR, gsl_vector $tau, gsl_permutation $p, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_svx(gsl_matrix $QR, gsl_vector $tau, gsl_permutation $p, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_lssolve(gsl_matrix $QR, gsl_vector $tau, gsl_permutation $p, gsl_vector $b, gsl_vector $x, gsl_vector $residual --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_lssolve2(gsl_matrix $QR, gsl_vector $tau, gsl_permutation $p, gsl_vector $b, size_t $rank, gsl_vector $x, gsl_vector $residual --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_QRsolve(gsl_matrix $Q, gsl_matrix $R, gsl_permutation $p, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_update(gsl_matrix $Q, gsl_matrix $R, gsl_permutation $p, gsl_vector $w, gsl_vector $v --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_Rsolve(gsl_matrix $QR, gsl_permutation $p, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_Rsvx(gsl_matrix $QR, gsl_permutation $p, gsl_vector $x --> int32) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_rank(gsl_matrix $QR, num64 $tol --> size_t) is native(&LIB) is export(:QR) { * }
sub gsl_linalg_QRPT_rcond(gsl_matrix $QR, num64 $rcond is rw, gsl_vector $work --> size_t) is native(&LIB) is export(:QR) { * }
# LQ Decomposition
sub gsl_linalg_LQ_decomp(gsl_matrix $A, gsl_vector $tau --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_LQ_solve_T(gsl_matrix $LQ, gsl_vector $tau, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_LQ_svx_T(gsl_matrix $LQ, gsl_vector $tau, gsl_vector $x --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_LQ_lssolve_T(gsl_matrix $LQ, gsl_vector $tau, gsl_vector $b, gsl_vector $x, gsl_vector $residual --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_LQ_Lsolve_T(gsl_matrix $LQ, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_LQ_Lsvx_T(gsl_matrix $LQ, gsl_vector $x --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_L_solve_T(gsl_matrix $L, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_LQ_vecQ(gsl_matrix $LQ, gsl_vector $tau, gsl_vector $v --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_LQ_vecQT(gsl_matrix $LQ, gsl_vector $tau, gsl_vector $v --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_LQ_unpack(gsl_matrix $LQ, gsl_vector $tau, gsl_matrix $Q, gsl_matrix $L --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_LQ_update(gsl_matrix $Q, gsl_matrix $R, gsl_vector $v, gsl_vector $w --> int32) is native(&LIB) is export(:LQ) { * }
sub gsl_linalg_LQ_LQsolve(gsl_matrix $Q, gsl_matrix $L, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:LQ) { * }
# Complete Orthogonal Decomposition
sub gsl_linalg_COD_decomp(gsl_matrix $A, gsl_vector $tau_Q, gsl_vector $tau_Z, gsl_permutation $p, size_t $rank is rw, gsl_vector $work --> int32) is native(&LIB) is export(:COD) { * }
sub gsl_linalg_COD_decomp_e(gsl_matrix $A, gsl_vector $tau_Q, gsl_vector $tau_Z, gsl_permutation $p, num64 $tol, size_t $rank is rw, gsl_vector $work --> int32) is native(&LIB) is export(:COD) { * }
sub gsl_linalg_COD_lssolve(gsl_matrix $QRZT, gsl_vector $tau_Q, gsl_vector $tau_Z, gsl_permutation $p, size_t $rank, gsl_vector $b, gsl_vector $x, gsl_vector $residual --> int32) is native(&LIB) is export(:COD) { * }
sub gsl_linalg_COD_lssolve2(num64 $lambda, gsl_matrix $QRZT, gsl_vector $tau_Q, gsl_vector $tau_Z, gsl_permutation $p, size_t $rank, gsl_vector $b, gsl_vector $x, gsl_vector $residual, gsl_matrix $S, gsl_vector $work --> int32) is native(&LIB) is export(:COD) { * }
sub gsl_linalg_COD_unpack(gsl_matrix $QRZT, gsl_vector $tau_Q, gsl_vector $tau_Z, size_t $rank, gsl_matrix $Q, gsl_matrix $R, gsl_matrix $Z --> int32) is native(&LIB) is export(:COD) { * }
sub gsl_linalg_COD_matZ(gsl_matrix $QRZT, gsl_vector $tau_Z, size_t $rank, gsl_matrix $A, gsl_vector $work --> int32) is native(&LIB) is export(:COD) { * }
# Singular Value Decomposition
sub gsl_linalg_SV_decomp(gsl_matrix $A, gsl_matrix $V, gsl_vector $S, gsl_vector $work --> int32) is native(&LIB) is export(:SV) { * }
sub gsl_linalg_SV_decomp_mod(gsl_matrix $A, gsl_matrix $X, gsl_matrix $V, gsl_vector $S, gsl_vector $work --> int32) is native(&LIB) is export(:SV) { * }
sub gsl_linalg_SV_decomp_jacobi(gsl_matrix $A, gsl_matrix $V, gsl_vector $S --> int32) is native(&LIB) is export(:SV) { * }
sub gsl_linalg_SV_solve(gsl_matrix $U, gsl_matrix $V, gsl_vector $S, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:SV) { * }
sub gsl_linalg_SV_leverage(gsl_matrix $U, gsl_vector $h --> int32) is native(&LIB) is export(:SV) { * }
# Cholesky Decomposition
sub gsl_linalg_cholesky_decomp1(gsl_matrix $A --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_complex_cholesky_decomp(gsl_matrix_complex $A --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_solve(gsl_matrix $cholesky, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_complex_cholesky_solve(gsl_matrix_complex $cholesky, gsl_vector_complex $b, gsl_vector_complex $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_solve_mat(gsl_matrix $cholesky, gsl_matrix $B, gsl_matrix $X --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_svx(gsl_matrix $cholesky, gsl_vector $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_svx_mat(gsl_matrix $cholesky, gsl_matrix $X --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_complex_cholesky_svx(gsl_matrix_complex $cholesky, gsl_vector_complex $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_invert(gsl_matrix $cholesky --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_complex_cholesky_invert(gsl_matrix_complex $cholesky --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_decomp2(gsl_matrix $A, gsl_vector $S --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_solve2(gsl_matrix $cholesky, gsl_vector $S, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_svx2(gsl_matrix $cholesky, gsl_vector $S, gsl_vector $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_decomp_unit(gsl_matrix $A, gsl_vector $D --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_scale(gsl_matrix $A, gsl_vector $S --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_scale_apply(gsl_matrix $A, gsl_vector $S --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_cholesky_rcond(gsl_matrix $cholesky, num64 $rcond is rw, gsl_vector $work --> int32) is native(&LIB) is export(:cholesky) { * }
# Pivoted Cholesky Decomposition
sub gsl_linalg_pcholesky_decomp(gsl_matrix $A, gsl_permutation $p --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_pcholesky_solve(gsl_matrix $LDLT, gsl_permutation $p, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_pcholesky_svx(gsl_matrix $LDLT, gsl_permutation $p, gsl_vector $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_pcholesky_decomp2(gsl_matrix $A, gsl_permutation $p, gsl_vector $S --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_pcholesky_solve2(gsl_matrix $LDLT, gsl_permutation $p, gsl_vector $S, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_pcholesky_svx2(gsl_matrix $LDLT, gsl_permutation $p, gsl_vector $S, gsl_vector $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_pcholesky_invert(gsl_matrix $LDLT, gsl_permutation $p, gsl_matrix $Ainv --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_pcholesky_rcond(gsl_matrix $LDLT, gsl_permutation $p, num64 $rcond is rw, gsl_vector $work --> int32) is native(&LIB) is export(:cholesky) { * }
# Modified Cholesky Decomposition
sub gsl_linalg_mcholesky_decomp(gsl_matrix $A, gsl_permutation $p, gsl_vector $E --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_mcholesky_solve(gsl_matrix $LDLT, gsl_permutation $p, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_mcholesky_svx(gsl_matrix $LDLT, gsl_permutation $p, gsl_vector $x --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_mcholesky_rcond(gsl_matrix $LDLT, gsl_permutation $p, num64 $rcond is rw, gsl_vector $work --> int32) is native(&LIB) is export(:cholesky) { * }
sub gsl_linalg_mcholesky_invert(gsl_matrix $LDLT, gsl_permutation $p, gsl_matrix $Ainv --> int32) is native(&LIB) is export(:cholesky) { * }
# LDLT Decomposition (v2.6)
sub gsl_linalg_ldlt_decomp(gsl_matrix $A --> int32) is native(&LIB) is export(:LDLT) { * }
sub gsl_linalg_ldlt_solve(gsl_matrix $LDLT, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:LDLT) { * }
sub gsl_linalg_ldlt_svx(gsl_matrix $LDLT, gsl_vector $x --> int32) is native(&LIB) is export(:LDLT) { * }
sub gsl_linalg_ldlt_rcond(gsl_matrix $LDLT, num64 $rcond is rw, gsl_vector $work --> int32) is native(&LIB) is export(:LDLT) { * }
# Tridiagonal Decomposition of Real Symmetric Matrices
sub gsl_linalg_symmtd_decomp(gsl_matrix $A, gsl_vector $tau --> int32) is native(&LIB) is export(:tridiag) { * }
sub gsl_linalg_symmtd_unpack(gsl_matrix $A, gsl_vector $tau, gsl_matrix $Q, gsl_vector $diag, gsl_vector $subdiag --> int32) is native(&LIB) is export(:tridiag) { * }
sub gsl_linalg_symmtd_unpack_T(gsl_matrix $A, gsl_vector $diag, gsl_vector $subdiag --> int32) is native(&LIB) is export(:tridiag) { * }
# Tridiagonal Decomposition of Hermitian Matrices
sub gsl_linalg_hermtd_decomp(gsl_matrix_complex $A, gsl_vector_complex $tau --> int32) is native(&LIB) is export(:tridiag) { * }
sub gsl_linalg_hermtd_unpack(gsl_matrix_complex $A, gsl_vector_complex $tau, gsl_matrix_complex $U, gsl_vector $diag, gsl_vector $subdiag --> int32) is native(&LIB) is export(:tridiag) { * }
sub gsl_linalg_hermtd_unpack_T(gsl_matrix_complex $A, gsl_vector $diag, gsl_vector $subdiag --> int32) is native(&LIB) is export(:tridiag) { * }
# Hessenberg Decomposition of Real Matrices
sub gsl_linalg_hessenberg_decomp(gsl_matrix $A, gsl_vector $tau --> int32) is native(&LIB) is export(:hessenberg) { * }
sub gsl_linalg_hessenberg_unpack(gsl_matrix $H, gsl_vector $tau, gsl_matrix $U --> int32) is native(&LIB) is export(:hessenberg) { * }
sub gsl_linalg_hessenberg_unpack_accum(gsl_matrix $H, gsl_vector $tau, gsl_matrix $V --> int32) is native(&LIB) is export(:hessenberg) { * }
sub gsl_linalg_hessenberg_set_zero(gsl_matrix $H --> int32) is native(&LIB) is export(:hessenberg) { * }
# Hessenberg-Triangular Decomposition of Real Matrices
sub gsl_linalg_hesstri_decomp(gsl_matrix $A, gsl_matrix $B, gsl_matrix $U, gsl_matrix $V, gsl_vector $work --> int32) is native(&LIB) is export(:hessenberg) { * }
# Bidiagonalization
sub gsl_linalg_bidiag_decomp(gsl_matrix $A, gsl_vector $tau_U, gsl_vector $tau_V --> int32) is native(&LIB) is export(:bidiag) { * }
sub gsl_linalg_bidiag_unpack(gsl_matrix $A, gsl_vector $tau_U, gsl_matrix $U, gsl_vector $tau_V, gsl_matrix $V, gsl_vector $diag, gsl_vector $superdiag --> int32) is native(&LIB) is export(:bidiag) { * }
sub gsl_linalg_bidiag_unpack2(gsl_matrix $A, gsl_vector $tau_U, gsl_vector $tau_V, gsl_matrix $V --> int32) is native(&LIB) is export(:bidiag) { * }
sub gsl_linalg_bidiag_unpack_B(gsl_matrix $A, gsl_vector $diag, gsl_vector $superdiag --> int32) is native(&LIB) is export(:bidiag) { * }
# Givens Rotations
sub gsl_linalg_givens(num64 $a, num64 $b, num64 $c is rw, num64 $s is rw) is native(&LIB) is export(:givens) { * }
sub gsl_linalg_givens_gv(gsl_vector $v, size_t $i, size_t $j, num64 $c, num64 $s) is native(&LIB) is export(:givens) { * }
# Householder Transformations
sub gsl_linalg_householder_transform(gsl_vector $w --> num64) is native(&LIB) is export(:householder) { * }
sub mgsl_linalg_complex_householder_transform(gsl_vector_complex $w, gsl_complex $res) is native(GSLHELPER) is export(:householder) { * }
sub gsl_linalg_householder_hm(num64 $tau, gsl_vector $v, gsl_matrix $A --> int32) is native(&LIB) is export(:householder) { * }
sub mgsl_linalg_complex_householder_hm(gsl_complex $tau, gsl_vector_complex $v, gsl_matrix_complex $A --> int32) is native(GSLHELPER) is export(:householder) { * }
sub gsl_linalg_householder_mh(num64 $tau, gsl_vector $v, gsl_matrix $A --> int32) is native(&LIB) is export(:householder) { * }
sub mgsl_linalg_complex_householder_mh(gsl_complex $tau, gsl_vector_complex $v, gsl_matrix_complex $A --> int32) is native(GSLHELPER) is export(:householder) { * }
sub gsl_linalg_householder_hv(num64 $tau, gsl_vector $v, gsl_vector $w --> int32) is native(&LIB) is export(:householder) { * }
sub mgsl_linalg_complex_householder_hv(gsl_complex $tau, gsl_vector_complex $v, gsl_vector_complex $w --> int32) is native(GSLHELPER) is export(:householder) { * }
# Householder solver for linear systems
sub gsl_linalg_HH_solve(gsl_matrix $A, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:householder) { * }
sub gsl_linalg_HH_svx(gsl_matrix $A, gsl_vector $x --> int32) is native(&LIB) is export(:householder) { * }
# Tridiagonal Systems
sub gsl_linalg_solve_tridiag(gsl_vector $diag, gsl_vector $e, gsl_vector $f, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:tridiag) { * }
sub gsl_linalg_solve_symm_tridiag(gsl_vector $diag, gsl_vector $e, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:tridiag) { * }
sub gsl_linalg_solve_cyc_tridiag(gsl_vector $diag, gsl_vector $e, gsl_vector $f, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:tridiag) { * }
sub gsl_linalg_solve_symm_cyc_tridiag(gsl_vector $diag, gsl_vector $e, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:tridiag) { * }
# Triangular Systems
# v2.6
sub gsl_linalg_tri_invert(int32 $Uplo, int32 $Diag, gsl_matrix $T --> int32) is native(&LIB) is export(:triang) { * }
sub gsl_linalg_complex_tri_invert(int32 $Uplo, int32 $Diag, gsl_matrix_complex $T --> int32) is native(&LIB) is export(:triang) { * }
sub gsl_linalg_tri_LTL(gsl_matrix $L --> int32) is native(&LIB) is export(:triang) { * }
sub gsl_linalg_complex_tri_LHL(gsl_matrix_complex $L --> int32) is native(&LIB) is export(:triang) { * }
sub gsl_linalg_tri_UL(gsl_matrix $LU --> int32) is native(&LIB) is export(:triang) { * }
sub gsl_linalg_complex_tri_UL(gsl_matrix_complex $LU --> int32) is native(&LIB) is export(:triang) { * }
sub gsl_linalg_tri_rcond(int32 $Uplo, gsl_matrix $A, num64 $rcond is rw, gsl_vector $work --> int32) is native(&LIB) is export(:triang) { * }
# v2.5 & v2.6
sub gsl_linalg_tri_upper_rcond(gsl_matrix $A, num64 $rcond is rw, gsl_vector $work --> int32) is native(&LIB) is export(:banlance) { * }
sub gsl_linalg_tri_lower_rcond(gsl_matrix $A, num64 $rcond is rw, gsl_vector $work --> int32) is native(&LIB) is export(:banlance) { * }
sub gsl_linalg_tri_upper_invert(gsl_matrix $T --> int32) is native(&LIB) is export(:triang) { * }
sub gsl_linalg_tri_lower_invert(gsl_matrix $T --> int32) is native(&LIB) is export(:triang) { * }
sub gsl_linalg_tri_upper_unit_invert(gsl_matrix $T --> int32) is native(&LIB) is export(:triang) { * }
sub gsl_linalg_tri_lower_unit_invert(gsl_matrix $T --> int32) is native(&LIB) is export(:triang) { * }
# Banded Cholesky Decomposition v2.6
sub gsl_linalg_cholesky_band_decomp(gsl_matrix $A --> int32) is native(&LIB) is export(:band) { * }
sub gsl_linalg_cholesky_band_solve(gsl_matrix $LLT, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:band) { * }
sub gsl_linalg_cholesky_band_svx(gsl_matrix $LLT, gsl_vector $x --> int32) is native(&LIB) is export(:band) { * }
sub gsl_linalg_cholesky_band_invert(gsl_matrix $LLT, gsl_matrix $Ainv --> int32) is native(&LIB) is export(:band) { * }
sub gsl_linalg_cholesky_band_unpack(gsl_matrix $LLT, gsl_matrix $L --> int32) is native(&LIB) is export(:band) { * }
sub gsl_linalg_cholesky_band_rcond(gsl_matrix $LLT, num64 $rcond is rw, gsl_vector $work --> int32) is native(&LIB) is export(:band) { * }
# Banded LDLT Decomposition v2.6
sub gsl_linalg_ldlt_band_decomp(gsl_matrix $A --> int32) is native(&LIB) is export(:band) { * }
sub gsl_linalg_ldlt_band_solve(gsl_matrix $LDLT, gsl_vector $b, gsl_vector $x --> int32) is native(&LIB) is export(:band) { * }
sub gsl_linalg_ldlt_band_svx(gsl_matrix $LDLT, gsl_vector $x --> int32) is native(&LIB) is export(:band) { * }
sub gsl_linalg_ldlt_band_unpack(gsl_matrix $LDLT, gsl_matrix $L, gsl_vector $D --> int32) is native(&LIB) is export(:band) { * }
sub gsl_linalg_ldlt_band_rcond(gsl_matrix $LDLT, num64 $rcond is rw, gsl_vector $work --> int32) is native(&LIB) is export(:band) { * }
# Balancing
sub gsl_linalg_balance_matrix(gsl_matrix $A, gsl_vector $D --> int32) is native(&LIB) is export(:banlance) { * }
