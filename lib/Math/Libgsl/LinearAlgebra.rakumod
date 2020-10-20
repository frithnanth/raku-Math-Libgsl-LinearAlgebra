use v6.c;

unit class Math::Libgsl::LinearAlgebra:ver<0.0.2>:auth<cpan:FRITH>;

use Math::Libgsl::Raw::LinearAlgebra :ALL;
use NativeCall;
use Math::Libgsl::Vector;
use Math::Libgsl::Matrix;
use Math::Libgsl::Permutation;
use Math::Libgsl::Exception;
use Math::Libgsl::Constants;

# LU Decomposition
sub LU-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> List) is export
{
  my Math::Libgsl::Permutation $p .= new: $A.matrix.size1;
  my int32 $s = 0;
  my $ret = gsl_linalg_LU_decomp($A.matrix, $p.p, $s);
  fail X::Libgsl.new: errno => $ret, error => "Error in LU-decomp" if $ret ≠ GSL_SUCCESS;
  return $s, $p;
}

sub LU-solve(Math::Libgsl::Matrix $LU where *.matrix.size1 == $LU.matrix.size2,
             Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1,
             Math::Libgsl::Vector $b where *.vector.size == $LU.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $LU.matrix.size2;
  my $ret = gsl_linalg_LU_solve($LU.matrix, $p.p, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in LU-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub LU-svx(Math::Libgsl::Matrix $LU where *.matrix.size1 == $LU.matrix.size2,
           Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1,
           Math::Libgsl::Vector $x where *.size == $LU.matrix.size1 --> Int) is export
{
  gsl_linalg_LU_svx($LU.matrix, $p.p, $x.vector);
}

sub LU-refine(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
              Math::Libgsl::Matrix $LU where { $LU.matrix.size1 == $LU.matrix.size2 && $A.matrix.size1 == $LU.matrix.size2 },
              Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1,
              Math::Libgsl::Vector $b where *.vector.size == $LU.matrix.size1,
              Math::Libgsl::Vector $x where *.vector.size == $LU.matrix.size1 --> Int) is export
{
  my Math::Libgsl::Vector $work .= new($LU.matrix.size1);
  gsl_linalg_LU_refine($A.matrix, $LU.matrix, $p.p, $b.vector, $x.vector, $work.vector);
}

sub LU-invert(Math::Libgsl::Matrix $LU,
              Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1 --> Math::Libgsl::Matrix) is export
{
  my Math::Libgsl::Matrix $inverse .= new($LU.matrix.size1, $LU.matrix.size1);
  my $ret = gsl_linalg_LU_invert($LU.matrix, $p.p, $inverse.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in LU-invert" if $ret ≠ GSL_SUCCESS;
  return $inverse;
}

sub LU-det(Math::Libgsl::Matrix $LU, Int $signum where * ~~ -1|1 --> Num) is export
{
  gsl_linalg_LU_det($LU.matrix, $signum);
}

sub LU-lndet(Math::Libgsl::Matrix $LU --> Num) is export
{
  gsl_linalg_LU_lndet($LU.matrix);
}

sub LU-sgndet(Math::Libgsl::Matrix $LU, Int $signum where * ~~ -1|1 --> Int) is export
{
  gsl_linalg_LU_sgndet($LU.matrix, $signum);
}

# QR Decomposition
sub QR-decomp(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $tau .= new: min($A.matrix.size1, $A.matrix.size2);
  my $ret = gsl_linalg_QR_decomp($A.matrix, $tau.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QR-decomp" if $ret ≠ GSL_SUCCESS;
  return $tau;
}

sub QR-solve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2,
             Math::Libgsl::Vector $tau where *.vector.size == $QR.matrix.size2,
             Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $QR.matrix.size2;
  my $ret = gsl_linalg_QR_solve($QR.matrix, $tau.vector, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QR-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub QR-svx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2,
           Math::Libgsl::Vector $tau where *.vector.size == $QR.matrix.size2,
           Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size1 --> Int) is export
{
  gsl_linalg_QR_svx($QR.matrix, $tau.vector, $x.vector);
}

sub QR-lssolve(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2,
               Math::Libgsl::Vector $tau where *.vector.size == $QR.matrix.size2,
               Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> List) is export
{
  my Math::Libgsl::Vector $x        .= new: $QR.matrix.size2;
  my Math::Libgsl::Vector $residual .= new: $QR.matrix.size1;
  my $ret = gsl_linalg_QR_lssolve($QR.matrix, $tau.vector, $b.vector, $x.vector, $residual.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QR-lssolve" if $ret ≠ GSL_SUCCESS;
  return $x, $residual;
}

sub QR-QTvec(Math::Libgsl::Matrix $QR,
             Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2),
             Math::Libgsl::Vector $v where *.vector.size == $QR.matrix.size1 --> Int) is export
{
  gsl_linalg_QR_QTvec($QR.matrix, $tau.vector, $v.vector);
}

sub QR-Qvec(Math::Libgsl::Matrix $QR,
            Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2),
            Math::Libgsl::Vector $v where *.vector.size == $QR.matrix.size1 --> Int) is export
{
  gsl_linalg_QR_Qvec($QR.matrix, $tau.vector, $v.vector);
}

sub QR-QTmat(Math::Libgsl::Matrix $QR,
             Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2),
             Math::Libgsl::Matrix $B where *.matrix.size1 == $QR.matrix.size1 --> Int) is export
{
  gsl_linalg_QR_QTmat($QR.matrix, $tau.vector, $B.matrix);
}

sub QR-Rsolve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2,
              Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $QR.matrix.size2;
  my $ret = gsl_linalg_QR_Rsolve($QR.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QR-Rsolve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub QR-Rsvx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2,
            Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size2 --> Int) is export
{
  gsl_linalg_QR_Rsvx($QR.matrix, $x.vector);
}

sub QR-unpack(Math::Libgsl::Matrix $QR,
              Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2)
              --> List) is export
{
  my Math::Libgsl::Matrix $Q .= new: $QR.matrix.size1, $QR.matrix.size1;
  my Math::Libgsl::Matrix $R .= new: $QR.matrix.size1, $QR.matrix.size2;
  my $ret = gsl_linalg_QR_unpack($QR.matrix, $tau.vector, $Q.matrix, $R.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in QR-unpack" if $ret ≠ GSL_SUCCESS;
  return $Q, $R;
}

sub QR-QRsolve(Math::Libgsl::Matrix $Q,
               Math::Libgsl::Matrix $R where { $R.matrix.size1 == $R.matrix.size2 && $Q.matrix.size1 == $R.matrix.size1 },
               Math::Libgsl::Vector $b where *.vector.size == $R.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $Q.matrix.size1;
  my $ret = gsl_linalg_QR_QRsolve($Q.matrix, $R.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QR-QRsolve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub QR-update(Math::Libgsl::Matrix $Q,
              Math::Libgsl::Matrix $R where { $Q.matrix.size1 == $R.matrix.size1 && $Q.matrix.size2 == $R.matrix.size1 },
              Math::Libgsl::Vector $w where *.vector.size == $R.matrix.size1,
              Math::Libgsl::Vector $v where *.vector.size == $R.matrix.size2 --> Int) is export
{
  gsl_linalg_QR_update($Q.matrix, $R.matrix, $w.vector, $v.vector);
}

sub R-solve(Math::Libgsl::Matrix $R where *.matrix.size1 == $R.matrix.size2,
            Math::Libgsl::Vector $b where *.vector.size == $R.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $R.matrix.size2;
  my $ret = gsl_linalg_R_solve($R.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in R-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub R-svx(Math::Libgsl::Matrix $R where *.matrix.size1 == $R.matrix.size2,
          Math::Libgsl::Vector $x where *.vector.size == $R.matrix.size2 --> Int) is export
{
  gsl_linalg_R_svx($R.matrix, $x.vector);
}

# QR Decomposition with Column Pivoting
sub QRPT-decomp(Math::Libgsl::Matrix $A --> List) is export
{
  my $M = $A.matrix.size1;
  my $N = $A.matrix.size2;
  my Math::Libgsl::Vector $tau    .= new: min($M, $N);
  my Math::Libgsl::Vector $norm   .= new: $N;
  my Math::Libgsl::Permutation $p .= new: $N;
  my int32 $signum;
  my $ret = gsl_linalg_QRPT_decomp($A.matrix, $tau.vector, $p.p, $signum, $norm.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QRPT-decomp" if $ret ≠ GSL_SUCCESS;
  return $tau, $p, $signum;
}

sub QRPT-decomp2(Math::Libgsl::Matrix $A --> List) is export
{
  my $M = $A.matrix.size1;
  my $N = $A.matrix.size2;
  my Math::Libgsl::Matrix $Q      .= new: $M, $M;
  my Math::Libgsl::Matrix $R      .= new: $M, $N;
  my Math::Libgsl::Vector $tau    .= new: min($M, $N);
  my Math::Libgsl::Vector $norm   .= new: $N;
  my Math::Libgsl::Permutation $p .= new: $N;
  my int32 $signum;
  my $ret = gsl_linalg_QRPT_decomp2($A.matrix, $Q.matrix, $R.matrix, $tau.vector, $p.p, $signum, $norm.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QRPT-decomp2" if $ret ≠ GSL_SUCCESS;
  return $Q, $R, $p, $signum;
}

sub QRPT-solve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2,
               Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2),
               Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size1,
               Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $QR.matrix.size2;
  my $ret = gsl_linalg_QRPT_solve($QR.matrix, $tau.vector, $p.p, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QRPT-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub QRPT-svx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2,
             Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2),
             Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size1,
             Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size2 --> Int) is export
{
  gsl_linalg_QRPT_svx($QR.matrix, $tau.vector, $p.p, $x.vector);
}

sub QRPT-lssolve(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2,
               Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2),
               Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2,
               Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> List) is export
{
  my Math::Libgsl::Vector $x        .= new: $QR.matrix.size2;
  my Math::Libgsl::Vector $residual .= new: $QR.matrix.size1;
  my $ret = gsl_linalg_QRPT_lssolve($QR.matrix, $tau.vector, $p.p, $b.vector, $x.vector, $residual.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QRPT-lssolve" if $ret ≠ GSL_SUCCESS;
  return $x, $residual;
}

sub QRPT-lssolve2(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2,
               Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2),
               Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2,
               Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1,
               Int $rank where 0 < * ≤ $QR.matrix.size2 --> List) is export
{
  my Math::Libgsl::Vector $x        .= new: $QR.matrix.size2;
  my Math::Libgsl::Vector $residual .= new: $QR.matrix.size1;
  my $ret = gsl_linalg_QRPT_lssolve2($QR.matrix, $tau.vector, $p.p, $b.vector, $rank, $x.vector, $residual.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QRPT-lssolve2" if $ret ≠ GSL_SUCCESS;
  return $x, $residual;
}

sub QRPT-QRsolve(Math::Libgsl::Matrix $Q where *.matrix.size1 == $Q.matrix.size2,
                 Math::Libgsl::Matrix $R where {$R.matrix.size1 == $R.matrix.size2 && $R.matrix.size1 == $Q.matrix.size1},
                 Math::Libgsl::Permutation $p where *.p.size == $Q.matrix.size1,
                 Math::Libgsl::Vector $b where *.vector.size == $Q.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $Q.matrix.size1;
  my $ret = gsl_linalg_QRPT_QRsolve($Q.matrix, $R.matrix, $p.p, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QRPT-QRsolve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub QRPT-update(Math::Libgsl::Matrix $Q,
                Math::Libgsl::Matrix $R where { $Q.matrix.size1 == $R.matrix.size1 && $Q.matrix.size2 == $R.matrix.size1 },
                Math::Libgsl::Permutation $p where *.p.size == $R.matrix.size1,
                Math::Libgsl::Vector $w where *.vector.size == $R.matrix.size1,
                Math::Libgsl::Vector $v where *.vector.size == $R.matrix.size2 --> Int) is export
{
  gsl_linalg_QRPT_update($Q.matrix, $R.matrix, $p.p, $w.vector, $v.vector);
}

sub QRPT-Rsolve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2,
                Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2,
                Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $QR.matrix.size2;
  my $ret = gsl_linalg_QRPT_Rsolve($QR.matrix, $p.p, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in QRPT-Rsolve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub QRPT-Rsvx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2,
              Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2,
              Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size2 --> Int) is export
{
  gsl_linalg_QRPT_Rsvx($QR.matrix, $p.p, $x.vector);
}

sub QRPT-rank(Math::Libgsl::Matrix $QR, Num() $tolerance --> Int) is export
{
  my num64 $tol = $tolerance;
  gsl_linalg_QRPT_rank($QR.matrix, $tol);
}

sub QRPT-rcond(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2 --> Num) is export
{
  my num64 $rcond;
  my Math::Libgsl::Vector $work .= new($QR.matrix.size2 * 3);
  my $ret = gsl_linalg_QRPT_rcond($QR.matrix, $rcond, $work);
  fail X::Libgsl.new: errno => $ret, error => "Error in QRPT-rcond" if $ret ≠ GSL_SUCCESS;
  return $rcond;
}

# LQ Decomposition
sub LQ-decomp(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $tau .= new: min($A.matrix.size1, $A.matrix.size2);
  my $ret = gsl_linalg_LQ_decomp($A.matrix, $tau.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in LQ-decomp" if $ret ≠ GSL_SUCCESS;
  return $tau;
}

sub LQ-solve-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2,
               Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2),
               Math::Libgsl::Vector $b where *.vector.size == $LQ.matrix.size2 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $LQ.matrix.size1;
  my $ret = gsl_linalg_LQ_solve_T($LQ.matrix, $tau.vector, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in LQ-solve-T" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub LQ-svx-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2,
               Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2),
               --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $LQ.matrix.size1;
  my $ret = gsl_linalg_LQ_svx_T($LQ.matrix, $tau.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in LQ-svx-T" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub LQ-lssolve-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 ≥ $LQ.matrix.size2,
               Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2),
               Math::Libgsl::Vector $b where *.vector.size == $LQ.matrix.size1 --> List) is export
{
  my Math::Libgsl::Vector $x        .= new: $LQ.matrix.size2;
  my Math::Libgsl::Vector $residual .= new: $LQ.matrix.size1;
  my $ret = gsl_linalg_LQ_lssolve_T($LQ.matrix, $tau.vector, $b.vector, $x.vector, $residual.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in LQ-lssolve-T" if $ret ≠ GSL_SUCCESS;
  return $x, $residual;
}

sub LQ-Lsolve-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2,
                Math::Libgsl::Vector $b where *.vector.size == $LQ.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $LQ.matrix.size1;
  my $ret = gsl_linalg_LQ_Lsolve_T($LQ.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in LQ-Lsolve-T" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub LQ-Lsvx-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $LQ.matrix.size1;
  my $ret = gsl_linalg_LQ_Lsvx_T($LQ.matrix, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in LQ-Lsvx-T" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub L-solve-T(Math::Libgsl::Matrix $L where *.matrix.size1 == $L.matrix.size2,
              Math::Libgsl::Vector $b where *.vector.size == $L.matrix.size2 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $L.matrix.size1;
  my $ret = gsl_linalg_L_solve_T($L.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in L-solve-T" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub LQ-vecQ(Math::Libgsl::Matrix $LQ,
            Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2),
            Math::Libgsl::Vector $v where *.vector.size == $LQ.matrix.size1 --> Int) is export
{
  gsl_linalg_LQ_vecQ($LQ.matrix, $tau.vector, $v.vector);
}

sub LQ-vecQT(Math::Libgsl::Matrix $LQ,
             Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2),
             Math::Libgsl::Vector $v where *.vector.size == $LQ.matrix.size1 --> Int) is export
{
  gsl_linalg_LQ_vecQT($LQ.matrix, $tau.vector, $v.vector);
}

sub LQ-unpack(Math::Libgsl::Matrix $LQ,
            Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2)
            --> List) is export
{
  my Math::Libgsl::Matrix $Q .= new: $LQ.matrix.size2, $LQ.matrix.size2;
  my Math::Libgsl::Matrix $L .= new: $LQ.matrix.size1, $LQ.matrix.size2;
  my $ret = gsl_linalg_LQ_unpack($LQ.matrix, $tau.vector, $Q.matrix, $L.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in LQ-unpack" if $ret ≠ GSL_SUCCESS;
  return $Q, $L;
}

sub LQ-update(Math::Libgsl::Matrix $Q,
              Math::Libgsl::Matrix $L where { $L.matrix.size2 == $Q.matrix.size1 && $L.matrix.size2 == $Q.matrix.size2 },
              Math::Libgsl::Vector $v where *.vector.size == $L.matrix.size1,
              Math::Libgsl::Vector $w where *.vector.size == $L.matrix.size2
              --> Int) is export
{
  gsl_linalg_LQ_update($Q.matrix, $L.matrix, $v.vector, $w.vector);
}

sub LQ-LQsolve(Math::Libgsl::Matrix $Q,
               Math::Libgsl::Matrix $L where { $L.matrix.size1 == $L.matrix.size2 && $Q.matrix.size1 == $L.matrix.size2 },
               Math::Libgsl::Vector $b where *.vector.size == $L.matrix.size2 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $L.matrix.size2;
  my $ret = gsl_linalg_LQ_LQsolve($Q.matrix, $L.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in LQ-LQsolve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

# Complete Orthogonal Decomposition
sub COD-decomp(Math::Libgsl::Matrix $A --> List) is export
{
  my $M = $A.matrix.size1;
  my $N = $A.matrix.size2;
  my Math::Libgsl::Vector $tau-Q .= new: min($M, $N);
  my Math::Libgsl::Vector $tau-Z .= new: min($M, $N);
  my Math::Libgsl::Vector $work   .= new: $N;
  my Math::Libgsl::Permutation $p .= new: $N;
  my size_t $rank;
  my $ret = gsl_linalg_COD_decomp($A.matrix, $tau-Q.vector, $tau-Z.vector, $p.p, $rank, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in COD-decomp" if $ret ≠ GSL_SUCCESS;
  return $tau-Q, $tau-Z, $p, $rank;
}

sub COD-decomp-e(Math::Libgsl::Matrix $A, Num() $tolerance --> List) is export
{
  my $M = $A.matrix.size1;
  my $N = $A.matrix.size2;
  my Math::Libgsl::Vector $tau-Q .= new: min($M, $N);
  my Math::Libgsl::Vector $tau-Z .= new: min($M, $N);
  my Math::Libgsl::Vector $work   .= new: $N;
  my Math::Libgsl::Permutation $p .= new: $N;
  my Int $rank;
  my num64 $tol = $tolerance;
  my $ret = gsl_linalg_COD_decomp_e($A.matrix, $tau-Q.vector, $tau-Z.vector, $p.p, $tol, $rank, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in COD-decomp-e" if $ret ≠ GSL_SUCCESS;
  return $tau-Q, $tau-Z, $p, $rank;
}

sub COD-lssolve(Math::Libgsl::Matrix $QRZT where *.matrix.size1 ≥ $QRZT.matrix.size2,
                Math::Libgsl::Vector $tau-Q where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2),
                Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2),
                Math::Libgsl::Permutation $p where *.p.size == $QRZT.matrix.size2,
                Int $rank where * ≤ min($QRZT.matrix.size1, $QRZT.matrix.size2),
                Math::Libgsl::Vector $b where *.vector.size == $QRZT.matrix.size1 --> List) is export
{
  my Math::Libgsl::Vector $x        .= new: $QRZT.matrix.size2;
  my Math::Libgsl::Vector $residual .= new: $QRZT.matrix.size1;
  my $ret = gsl_linalg_COD_lssolve($QRZT.matrix, $tau-Q.vector, $tau-Z.vector, $p.p, $rank,
                                   $b.vector, $x.vector, $residual.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in COD-lssolve" if $ret ≠ GSL_SUCCESS;
  return $x, $residual;
}

sub COD-lssolve2(Math::Libgsl::Matrix $QRZT where *.matrix.size1 ≥ $QRZT.matrix.size2,
                 Math::Libgsl::Vector $tau-Q where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2),
                 Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2),
                 Math::Libgsl::Permutation $p where *.p.size == $QRZT.matrix.size2,
                 Int $rank where * ≤ min($QRZT.matrix.size1, $QRZT.matrix.size2),
                 Math::Libgsl::Vector $b where *.vector.size == $QRZT.matrix.size1,
                 Num() $lambda --> List) is export
{
  my Math::Libgsl::Vector $x        .= new: $QRZT.matrix.size2;
  my Math::Libgsl::Vector $residual .= new: $QRZT.matrix.size1;
  my Math::Libgsl::Matrix $S        .= new: $rank, $rank;
  my Math::Libgsl::Vector $work     .= new: $rank;
  my $ret = gsl_linalg_COD_lssolve2($lambda, $QRZT.matrix, $tau-Q.vector, $tau-Z.vector, $p.p, $rank,
                                    $b.vector, $x.vector, $residual.vector, $S.matrix, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in COD-lssolve2" if $ret ≠ GSL_SUCCESS;
  return $x, $residual;
}

sub COD-unpack(Math::Libgsl::Matrix $QRZT,
               Math::Libgsl::Vector $tau-Q where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2),
               Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2),
               Int $rank where * ≤ min($QRZT.matrix.size1, $QRZT.matrix.size2) --> List) is export
{
  my $M = $QRZT.matrix.size1;
  my $N = $QRZT.matrix.size2;
  my Math::Libgsl::Matrix $Q .= new: $M, $M;
  my Math::Libgsl::Matrix $R .= new: $M, $N;
  my Math::Libgsl::Matrix $Z .= new: $N, $N;
  my $ret = gsl_linalg_COD_unpack($QRZT.matrix, $tau-Q.vector, $tau-Z.vector, $rank, $Q.matrix, $R.matrix, $Z.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in COD-unpack" if $ret ≠ GSL_SUCCESS;
  return $Q, $R, $Z;
}

sub COD-matZ(Math::Libgsl::Matrix $QRZT,
             Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2),
             Math::Libgsl::Matrix $A where *.matrix.size2 == $QRZT.matrix.size2,
             Int $rank --> Int) is export
{
  my Math::Libgsl::Vector $work .= new: $A.matrix.size1;
  gsl_linalg_COD_matZ($QRZT.matrix, $tau-Z.vector, $rank, $A.matrix, $work.vector);
}

# Singular Value Decomposition
sub SV-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List) is export
{
  my $M = $A.matrix.size1;
  my $N = $A.matrix.size2;
  my Math::Libgsl::Matrix $V .= new: $N, $N;
  my Math::Libgsl::Vector $S .= new: $N;
  my Math::Libgsl::Vector $work   .= new: $N;
  my $ret = gsl_linalg_SV_decomp($A.matrix, $V.matrix, $S.vector, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in SV-decomp" if $ret ≠ GSL_SUCCESS;
  return $V, $S;
}

sub SV-decomp-mod(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List) is export
{
  my $M = $A.matrix.size1;
  my $N = $A.matrix.size2;
  my Math::Libgsl::Matrix $V    .= new: $N, $N;
  my Math::Libgsl::Matrix $X    .= new: $N, $N;
  my Math::Libgsl::Vector $S    .= new: $N;
  my Math::Libgsl::Vector $work .= new: $N;
  my $ret = gsl_linalg_SV_decomp_mod($A.matrix, $X.matrix, $V.matrix, $S.vector, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in SV-decomp-mod" if $ret ≠ GSL_SUCCESS;
  return $V, $S;
}

sub SV-decomp-jacobi(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List) is export
{
  my $M = $A.matrix.size1;
  my $N = $A.matrix.size2;
  my Math::Libgsl::Matrix $V    .= new: $N, $N;
  my Math::Libgsl::Vector $S    .= new: $N;
  my $ret = gsl_linalg_SV_decomp_jacobi($A.matrix, $V.matrix, $S.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in SV-decomp-jacobi" if $ret ≠ GSL_SUCCESS;
  return $V, $S;
}

sub SV-solve(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2,
             Math::Libgsl::Matrix $V where { $V.matrix.size1 == $V.matrix.size2 && $V.matrix.size1 == $A.matrix.size2 },
             Math::Libgsl::Vector $S where *.vector.size == $A.matrix.size2,
             Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1
             --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $A.matrix.size2;
  my $ret = gsl_linalg_SV_solve($A.matrix, $V.matrix, $S.vector, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in SV-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub SV-leverage(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $h .= new: $A.matrix.size1;
  my $ret = gsl_linalg_SV_leverage($A.matrix, $h.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in SV-leverage" if $ret ≠ GSL_SUCCESS;
  return $h;
}

# Cholesky Decomposition
sub cholesky-decomp1(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int) is export
{
  gsl_linalg_cholesky_decomp1($A.matrix);
}

sub cholesky-solve(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
                   Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $A.matrix.size2;
  my $ret = gsl_linalg_cholesky_solve($A.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub cholesky-svx(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
                 Math::Libgsl::Vector $x where *.vector.size == $A.matrix.size1 --> Int) is export
{
  gsl_linalg_cholesky_svx($A.matrix, $x.vector);
}

sub cholesky-solve-mat(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
                       Math::Libgsl::Matrix $B where *.matrix.size1 == $A.matrix.size1
                       --> Math::Libgsl::Matrix) is export
{
  my Math::Libgsl::Matrix $X .= new: $B.matrix.size1, $B.matrix.size2;
  my $ret = gsl_linalg_cholesky_solve_mat($A.matrix, $B.matrix, $X.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-solve-mat" if $ret ≠ GSL_SUCCESS;
  return $X;
}

sub cholesky-svx-mat(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
                     Math::Libgsl::Matrix $X where *.matrix.size1 == $A.matrix.size1 --> Int) is export
{
  gsl_linalg_cholesky_svx_mat($A.matrix, $X.matrix);
}

sub cholesky-invert(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int) is export
{
  gsl_linalg_cholesky_invert($A.matrix);
}

sub cholesky-decomp2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $s .= new: $A.matrix.size1;
  my $ret = gsl_linalg_cholesky_decomp2($A.matrix, $s.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-decomp2" if $ret ≠ GSL_SUCCESS;
  return $s;
}

sub cholesky-solve2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
                    Math::Libgsl::Vector $s where *.vector.size == $A.matrix.size1,
                    Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $A.matrix.size2;
  my $ret = gsl_linalg_cholesky_solve2($A.matrix, $s.vector, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-solve2" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub cholesky-svx2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
                  Math::Libgsl::Vector $s where *.size == $A.matrix.size2,
                  Math::Libgsl::Vector $x where *.vector.size == $A.matrix.size2 --> Int) is export
{
  gsl_linalg_cholesky_svx2($A.matrix, $s.vector, $x.vector);
}

sub cholesky-decomp-unit(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $d .= new: $A.matrix.size1;
  my $ret = gsl_linalg_cholesky_decomp_unit($A.matrix, $d.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-decomp-unit" if $ret ≠ GSL_SUCCESS;
  return $d;
}

sub cholesky-scale(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $s .= new: $A.matrix.size2;
  my $ret = gsl_linalg_cholesky_scale($A.matrix, $s.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-scale" if $ret ≠ GSL_SUCCESS;
  return $s;
}

sub cholesky-scale-apply(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
                         Math::Libgsl::Vector $s where *.vector.size == $A.matrix.size2 --> Int) is export
{
  gsl_linalg_cholesky_scale_apply($A.matrix, $s.vector);
}

sub cholesky-rcond(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num) is export
{
  my Math::Libgsl::Vector $work .= new: 3 * $A.matrix.size2;
  my num64 $rcond;
  my $ret = gsl_linalg_cholesky_rcond($A.matrix, $rcond, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-rcond" if $ret ≠ GSL_SUCCESS;
  return $rcond;
}

# Pivoted Cholesky Decomposition
sub pcholesky-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2
                     --> Math::Libgsl::Permutation) is export
{
  my Math::Libgsl::Permutation $p .= new: $A.matrix.size1;
  my $ret = gsl_linalg_pcholesky_decomp($A.matrix, $p.p);
  fail X::Libgsl.new: errno => $ret, error => "Error in pcholesky-decomp" if $ret ≠ GSL_SUCCESS;
  return $p;
}

sub pcholesky-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                    Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1,
                    Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1
                    --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $LDLT.matrix.size2;
  my $ret = gsl_linalg_pcholesky_solve($LDLT.matrix, $p.p, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in pcholesky-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub pcholesky-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                  Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1,
                  Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int) is export
{
  gsl_linalg_pcholesky_svx($LDLT.matrix, $p.p, $x.vector);
}

sub pcholesky-decomp2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> List) is export
{
  my Math::Libgsl::Permutation $p .= new: $A.matrix.size1;
  my Math::Libgsl::Vector $s .= new: $A.matrix.size2;
  my $ret = gsl_linalg_pcholesky_decomp2($A.matrix, $p.p, $s.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in pcholesky-decomp2" if $ret ≠ GSL_SUCCESS;
  return $p, $s;
}

sub pcholesky-solve2(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                     Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1,
                     Math::Libgsl::Vector $s where *.vector.size == $LDLT.matrix.size1,
                     Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1
                     --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $LDLT.matrix.size2;
  my $ret = gsl_linalg_pcholesky_solve2($LDLT.matrix, $p.p, $s.vector, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in pcholesky-solve2" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub pcholesky-svx2(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                   Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1,
                   Math::Libgsl::Vector $s where *.vector.size == $LDLT.matrix.size1,
                   Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int) is export
{
  gsl_linalg_pcholesky_svx2($LDLT.matrix, $p.p, $s.vector, $x.vector);
}

sub pcholesky-invert(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                     Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1
                     --> Math::Libgsl::Matrix) is export
{
  my Math::Libgsl::Matrix $Ainv .= new: $LDLT.matrix.size1, $LDLT.matrix.size2;
  my $ret = gsl_linalg_pcholesky_invert($LDLT.matrix, $p.p, $Ainv.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in pcholesky-invert" if $ret ≠ GSL_SUCCESS;
  return $Ainv;
}

sub pcholesky-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                    Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1
                    --> Num) is export
{
  my Math::Libgsl::Vector $work .= new: 3 * $LDLT.matrix.size2;
  my num64 $rcond;
  my $ret = gsl_linalg_pcholesky_rcond($LDLT.matrix, $p.p, $rcond, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in pcholesky-rcond" if $ret ≠ GSL_SUCCESS;
  return $rcond;
}

# Modified Cholesky Decomposition
sub mcholesky-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
                     Bool :$perturbation = True
                     --> List) is export
{
  my Math::Libgsl::Permutation $p .= new: $A.matrix.size1;
  my Math::Libgsl::Vector $E;
  $E .= new($A.matrix.size1) if $perturbation;
  my $ret = gsl_linalg_mcholesky_decomp($A.matrix, $p.p, $perturbation ?? $E.vector !! Math::Libgsl::Vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in mcholesky-decomp" if $ret ≠ GSL_SUCCESS;
  return $p, $E;
}

sub mcholesky-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                    Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1,
                    Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1
                    --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $LDLT.matrix.size2;
  my $ret = gsl_linalg_mcholesky_solve($LDLT.matrix, $p.p, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in mcholesky-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub mcholesky-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                  Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1,
                  Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int) is export
{
  gsl_linalg_mcholesky_svx($LDLT.matrix, $p.p, $x.vector);
}

sub mcholesky-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                    Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1
                    --> Num) is export
{
  my Math::Libgsl::Vector $work .= new: 3 * $LDLT.matrix.size2;
  my num64 $rcond;
  my $ret = gsl_linalg_mcholesky_rcond($LDLT.matrix, $p.p, $rcond, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in mcholesky-rcond" if $ret ≠ GSL_SUCCESS;
  return $rcond;
}

sub mcholesky-invert(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                     Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1
                     --> Math::Libgsl::Matrix) is export
{
  my Math::Libgsl::Matrix $Ainv .= new: $LDLT.matrix.size1, $LDLT.matrix.size2;
  my $ret = gsl_linalg_mcholesky_invert($LDLT.matrix, $p.p, $Ainv.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in mcholesky-invert" if $ret ≠ GSL_SUCCESS;
  return $Ainv;
}

# LDLT Decomposition
sub ldlt-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-decomp: version < v2.6" if $gsl-version < 2.6;
  gsl_linalg_ldlt_decomp($A.matrix);
}

sub ldlt-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
               Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1
               --> Math::Libgsl::Vector) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-solve: version < v2.6" if $gsl-version < 2.6;
  my Math::Libgsl::Vector $x .= new: $LDLT.matrix.size2;
  my $ret = gsl_linalg_ldlt_solve($LDLT.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in ldlt-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub ldlt-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
             Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-svx: version < v2.6" if $gsl-version < 2.6;
  gsl_linalg_ldlt_svx($LDLT.matrix, $x.vector);
}

sub ldlt-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2 --> Num) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-rcond: version < v2.6" if $gsl-version < 2.6;
  my Math::Libgsl::Vector $work .= new: 3 * $LDLT.matrix.size2;
  my num64 $rcond;
  my $ret = gsl_linalg_ldlt_rcond($LDLT.matrix, $rcond, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in ldlt-rcond" if $ret ≠ GSL_SUCCESS;
  return $rcond;
}

# Tridiagonal Decomposition of Real Symmetric Matrices
sub symmtd-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $tau .= new: $A.matrix.size1 - 1;
  my $ret = gsl_linalg_symmtd_decomp($A.matrix, $tau.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in symmtd-decomp" if $ret ≠ GSL_SUCCESS;
  return $tau;
}

sub symmtd-unpack(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
                  Math::Libgsl::Vector $tau where *.vector.size == $A.matrix.size1 - 1 --> List) is export
{
  my $M = $A.matrix.size1;
  my Math::Libgsl::Matrix $Q     .= new: $M, $M;
  my Math::Libgsl::Vector $diag  .= new: $A.matrix.size1;
  my Math::Libgsl::Vector $sdiag .= new: $A.matrix.size1 - 1;
  my $ret = gsl_linalg_symmtd_unpack($A.matrix, $tau.vector, $Q.matrix, $diag.vector, $sdiag.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in symmtd-unpack" if $ret ≠ GSL_SUCCESS;
  return $Q, $diag, $sdiag;
}

sub symmtd-unpack-T(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> List) is export
{
  my $M = $A.matrix.size1;
  my Math::Libgsl::Vector $diag  .= new: $A.matrix.size1;
  my Math::Libgsl::Vector $sdiag .= new: $A.matrix.size1 - 1;
  my $ret = gsl_linalg_symmtd_unpack_T($A.matrix, $diag.vector, $sdiag.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in symmtd-unpack-T" if $ret ≠ GSL_SUCCESS;
  return $diag, $sdiag;
}

# Hessenberg Decomposition of Real Matrices
sub hessenberg-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $tau .= new: $A.matrix.size1;
  my $ret = gsl_linalg_hessenberg_decomp($A.matrix, $tau.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in hessenberg-decomp" if $ret ≠ GSL_SUCCESS;
  return $tau;
}

sub hessenberg-unpack(Math::Libgsl::Matrix $H where *.matrix.size1 == $H.matrix.size2,
                      Math::Libgsl::Vector $tau where *.vector.size == $H.matrix.size1
                      --> Math::Libgsl::Matrix) is export
{
  my $M = $H.matrix.size1;
  my Math::Libgsl::Matrix $U .= new: $M, $M;
  my $ret = gsl_linalg_hessenberg_unpack($H.matrix, $tau.vector, $U.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in hessenberg-unpack" if $ret ≠ GSL_SUCCESS;
  return $U;
}

sub hessenberg-unpack-accum(Math::Libgsl::Matrix $H where *.matrix.size1 == $H.matrix.size2,
                            Math::Libgsl::Vector $tau where *.vector.size == $H.matrix.size1
                            --> Math::Libgsl::Matrix) is export
{
  my $M = $H.matrix.size1;
  my Math::Libgsl::Matrix $V .= new: $M, $M;
  my $ret = gsl_linalg_hessenberg_unpack_accum($H.matrix, $tau.vector, $V.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in hessenberg-unpack-accum" if $ret ≠ GSL_SUCCESS;
  return $V;
}

sub hessenberg-set-zero(Math::Libgsl::Matrix $H where *.matrix.size1 == $H.matrix.size2 --> Int) is export
{
  gsl_linalg_hessenberg_set_zero($H.matrix);
}

# Hessenberg-Triangular Decomposition of Real Matrices
sub hesstri-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2,
                   Math::Libgsl::Matrix $B where { $B.matrix.size1 == $A.matrix.size1 && $B.matrix.size2 == $A.matrix.size2},
                   Bool :$similarity = True
                   --> List) is export
{
  my $M = $A.matrix.size1;
  my Math::Libgsl::Vector $work .= new: $M;
  my Math::Libgsl::Matrix $U;
  my Math::Libgsl::Matrix $V;
  if $similarity {
    $U .= new: $M, $M;
    $V .= new: $M, $M;
  }
  my $ret = gsl_linalg_hesstri_decomp($A.matrix,
                                      $B.matrix,
                                      $similarity ?? $U.matrix !! Math::Libgsl::Matrix,
                                      $similarity ?? $V.matrix !! Math::Libgsl::Matrix,
                                      $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in hesstri-decomp" if $ret ≠ GSL_SUCCESS;
  return $U, $V;
}

# Bidiagonalization
sub bidiag-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List) is export
{
  my Math::Libgsl::Vector $tau_U .= new: $A.matrix.size2;
  my Math::Libgsl::Vector $tau_V .= new: $A.matrix.size2 - 1;
  my $ret = gsl_linalg_bidiag_decomp($A.matrix, $tau_U.vector, $tau_V.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in bidiag-decomp" if $ret ≠ GSL_SUCCESS;
  return $tau_U, $tau_V;
}

sub bidiag-unpack(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2,
                  Math::Libgsl::Vector $tau_U where *.vector.size == $A.matrix.size2,
                  Math::Libgsl::Vector $tau_V where *.vector.size == $A.matrix.size2 - 1
                  --> List) is export
{
  my $M = $A.matrix.size1;
  my $N = $A.matrix.size2;
  my Math::Libgsl::Matrix $U     .= new: $M, $N;
  my Math::Libgsl::Matrix $V     .= new: $N, $N;
  my Math::Libgsl::Vector $diag  .= new: $N;
  my Math::Libgsl::Vector $sdiag .= new: $N - 1;
  my $ret = gsl_linalg_bidiag_unpack($A.matrix, $tau_U.vector, $U.matrix, $tau_V.vector, $V.matrix, $diag.vector, $sdiag.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in bidiag-unpack" if $ret ≠ GSL_SUCCESS;
  return $U, $V, $diag, $sdiag;
}

sub bidiag-unpack2(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2,
                   Math::Libgsl::Vector $tau_U where *.vector.size == $A.matrix.size2,
                   Math::Libgsl::Vector $tau_V where *.vector.size == $A.matrix.size2 - 1
                   --> Math::Libgsl::Matrix) is export
{
  my $M = $A.matrix.size1;
  my $N = $A.matrix.size2;
  my Math::Libgsl::Matrix $V .= new: $N, $N;
  my $ret = gsl_linalg_bidiag_unpack2($A.matrix, $tau_U.vector, $tau_V.vector, $V.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in bidiag-unpack2" if $ret ≠ GSL_SUCCESS;
  return $V;
}

sub bidiag-unpack-B(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List) is export
{
  my $N = $A.matrix.size2;
  my Math::Libgsl::Vector $diag  .= new: $N;
  my Math::Libgsl::Vector $sdiag .= new: $N - 1;
  my $ret = gsl_linalg_bidiag_unpack_B($A.matrix, $diag.vector, $sdiag.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in bidiag-unpack-B" if $ret ≠ GSL_SUCCESS;
  return $diag, $sdiag;
}

# Givens Rotations
sub givens(Num() $a, Num() $b --> List) is export
{
  my num64 $A = $a;
  my num64 $B = $b;
  my num64 ($c, $s);
  gsl_linalg_givens($A, $B, $c, $s);
  return $c, $s;
}

sub givens-gv(Math::Libgsl::Vector $v, Int, $i, Int $j, Num() $c, Num() $s) is export
{
  my num64 $C = $c;
  my num64 $S = $s;
  gsl_linalg_givens_gv($v.vector, $i, $j, $c, $s);
}

# Householder Transformations
sub householder-transform(Math::Libgsl::Vector $w --> Num) is export
{
  gsl_linalg_householder_transform($w.vector);
}

sub householder-hm(Num() $tau, Math::Libgsl::Vector $v, Math::Libgsl::Matrix $A --> Int) is export
{
  my num64 $τ = $tau;
  gsl_linalg_householder_hm($τ, $v.vector, $A.matrix);
}

sub householder-mh(Num() $tau, Math::Libgsl::Vector $v, Math::Libgsl::Matrix $A --> Int) is export
{
  my num64 $τ = $tau;
  gsl_linalg_householder_mh($τ, $v.vector, $A.matrix);
}

sub householder-hv(Num() $tau, Math::Libgsl::Vector $v, Math::Libgsl::Vector $w --> Int) is export
{
  my num64 $τ = $tau;
  gsl_linalg_householder_hv($τ, $v.vector, $w.vector);
}

# Householder solver for linear systems
sub HH-solve(Math::Libgsl::Matrix $A where *.matrix.size1 ≤ $A.matrix.size2,
             Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1
             --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $A.matrix.size2;
  my $ret = gsl_linalg_HH_solve($A.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in HH-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub HH-svx(Math::Libgsl::Matrix $A where *.matrix.size1 ≤ $A.matrix.size2,
           Math::Libgsl::Vector $x where *.vector.size == $A.matrix.size2 --> Int) is export
{
  gsl_linalg_HH_svx($A.matrix, $x.vector);
}

# Tridiagonal Systems
sub tridiag-solve(Math::Libgsl::Vector $diag,
                  Math::Libgsl::Vector $abovediag where *.vector.size == $diag.vector.size - 1,
                  Math::Libgsl::Vector $belowdiag where *.vector.size == $diag.vector.size - 1,
                  Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size
                  --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $diag.vector.size;
  my $ret = gsl_linalg_solve_tridiag($diag.vector, $abovediag.vector, $belowdiag.vector, $rhs.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in tridiag-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub tridiag-symm-solve(Math::Libgsl::Vector $diag,
                       Math::Libgsl::Vector $offdiag where *.vector.size == $diag.vector.size - 1,
                       Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size
                       --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $diag.vector.size;
  my $ret = gsl_linalg_solve_symm_tridiag($diag.vector, $offdiag.vector, $rhs.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in tridiag-symm-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub tridiag-cyc-solve(Math::Libgsl::Vector $diag where *.vector.size ≥ 3,
                      Math::Libgsl::Vector $abovediag where *.vector.size == $diag.vector.size,
                      Math::Libgsl::Vector $belowdiag where *.vector.size == $diag.vector.size,
                      Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size
                      --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $diag.vector.size;
  my $ret = gsl_linalg_solve_cyc_tridiag($diag.vector, $abovediag.vector, $belowdiag.vector, $rhs.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in tridiag-cyc-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub tridiag-symm-cyc-solve(Math::Libgsl::Vector $diag where *.vector.size ≥ 3,
                           Math::Libgsl::Vector $offdiag where *.vector.size == $diag.vector.size,
                           Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size
                           --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $x .= new: $diag.vector.size;
  my $ret = gsl_linalg_solve_symm_cyc_tridiag($diag.vector, $offdiag.vector, $rhs.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in tridiag-symm-cyc-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

# Triangular Systems
# v2.5 & v2.6
sub tri-upper-rcond(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num) is export
{
  my Math::Libgsl::Vector $work .= new: 3 * $A.matrix.size2;
  my num64 $rcond;
  my $ret = gsl_linalg_tri_upper_rcond($A.matrix, $rcond, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in tri-upper-rcond" if $ret ≠ GSL_SUCCESS;
  return $rcond;
}

sub tri-lower-rcond(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num) is export
{
  my Math::Libgsl::Vector $work .= new: 3 * $A.matrix.size2;
  my num64 $rcond;
  my $ret = gsl_linalg_tri_lower_rcond($A.matrix, $rcond, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in tri-lower-rcond" if $ret ≠ GSL_SUCCESS;
  return $rcond;
}

sub tri-upper-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int) is export
{
  gsl_linalg_tri_upper_invert($T.matrix);
}

sub tri-lower-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int) is export
{
  gsl_linalg_tri_lower_invert($T.matrix);
}

sub tri-upper-unit-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int) is export
{
  gsl_linalg_tri_upper_unit_invert($T.matrix);
}

sub tri-lower-unit-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int) is export
{
  gsl_linalg_tri_lower_unit_invert($T.matrix);
}

# v2.6
sub tri-invert(Int $Uplo, Int $Diag, Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in tri-invert: version < v2.6" if $gsl-version < 2.6;
  gsl_linalg_tri_invert($Uplo, $Diag, $T.matrix);
}

sub tri-LTL(Math::Libgsl::Matrix $L where *.matrix.size1 == $L.matrix.size2 --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in tri-LTL: version < v2.6" if $gsl-version < 2.6;
  gsl_linalg_tri_LTL($L.matrix);
}

sub tri-UL(Math::Libgsl::Matrix $LU where *.matrix.size1 == $LU.matrix.size2 --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in tri-UL: version < v2.6" if $gsl-version < 2.6;
  gsl_linalg_tri_UL($LU.matrix);
}

sub tri-rcond(Int $Uplo, Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in tri-rcond: version < v2.6" if $gsl-version < 2.6;
  my Math::Libgsl::Vector $work .= new: 3 * $A.matrix.size2;
  my num64 $rcond;
  my $ret = gsl_linalg_tri_rcond($Uplo, $A.matrix, $rcond, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in tri-rcond" if $ret ≠ GSL_SUCCESS;
  return $rcond;
}

# Banded Cholesky Decomposition v2.6
sub cholesky-band-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in cholesky-band-decomp: version < v2.6" if $gsl-version < 2.6;
  gsl_linalg_cholesky_band_decomp($A.matrix);
}

sub cholesky-band-solve(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2,
                        Math::Libgsl::Vector $b where *.vector.size == $LLT.matrix.size1
                        --> Math::Libgsl::Vector) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in cholesky-band-solve: version < v2.6" if $gsl-version < 2.6;
  my Math::Libgsl::Vector $x .= new: $b.vector.size;
  my $ret = gsl_linalg_cholesky_band_solve($LLT.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-band-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub cholesky-band-svx(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2,
                      Math::Libgsl::Vector $x where *.vector.size == $LLT.matrix.size1
                      --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in cholesky-band-svx: version < v2.6" if $gsl-version < 2.6;
  gsl_linalg_cholesky_band_svx($LLT.matrix, $x.vector);
}

sub cholesky-band-invert(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2
                         --> Math::Libgsl::Matrix) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in cholesky-band-invert: version < v2.6" if $gsl-version < 2.6;
  my Math::Libgsl::Matrix $Ainv .= new: $LLT.matrix.size1, $LLT.matrix.size2;
  my $ret = gsl_linalg_cholesky_band_invert($LLT.matrix, $Ainv.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-band-invert" if $ret ≠ GSL_SUCCESS;
  return $Ainv;
}

sub cholesky-band-unpack(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2
                         --> Math::Libgsl::Matrix) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in cholesky-band-unpack: version < v2.6" if $gsl-version < 2.6;
  my Math::Libgsl::Matrix $L .= new: $LLT.matrix.size1, $LLT.matrix.size2;
  my $ret = gsl_linalg_cholesky_band_unpack($LLT.matrix, $L.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-band-unpack" if $ret ≠ GSL_SUCCESS;
  return $L;
}

sub cholesky-band-rcond(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2 --> Num) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in cholesky-band-rcond: version < v2.6" if $gsl-version < 2.6;
  my Math::Libgsl::Vector $work .= new: 3 * $LLT.matrix.size2;
  my num64 $rcond;
  my $ret = gsl_linalg_cholesky_band_rcond($LLT.matrix, $rcond, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-band-rcond" if $ret ≠ GSL_SUCCESS;
  return $rcond;
}

# Banded LDLT Decomposition v2.6
sub ldlt-band-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-band-decomp: version < v2.6" if $gsl-version < 2.6;
  gsl_linalg_ldlt_band_decomp($A.matrix);
}

sub ldlt-band-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                    Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1
                    --> Math::Libgsl::Vector) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-band-solve: version < v2.6" if $gsl-version < 2.6;
  my Math::Libgsl::Vector $x .= new: $b.vector.size;
  my $ret = gsl_linalg_ldlt_band_solve($LDLT.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in ldlt-band-solve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub ldlt-band-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2,
                  Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1
                  --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-band-svx: version < v2.6" if $gsl-version < 2.6;
  gsl_linalg_ldlt_band_svx($LDLT.matrix, $x.vector);
}

sub ldlt-band-unpack(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2 --> List) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-band-unpack: version < v2.6" if $gsl-version < 2.6;
  my Math::Libgsl::Matrix $L .= new: $LDLT.matrix.size1, $LDLT.matrix.size2;
  my Math::Libgsl::Vector $D .= new: $LDLT.matrix.size1;
  my $ret = gsl_linalg_ldlt_band_unpack($LDLT.matrix, $L.matrix, $D.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in ldlt-band-unpack" if $ret ≠ GSL_SUCCESS;
  return $L, $D;
}

sub ldlt-band-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2 --> Num) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-band-rcond: version < v2.6" if $gsl-version < 2.6;
  my Math::Libgsl::Vector $work .= new: 3 * $LDLT.matrix.size2;
  my num64 $rcond;
  my $ret = gsl_linalg_ldlt_band_rcond($LDLT.matrix, $rcond, $work.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in ldlt-band-rcond" if $ret ≠ GSL_SUCCESS;
  return $rcond;
}

# Balancing
sub balance-matrix(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector) is export
{
  my Math::Libgsl::Vector $D .= new: $A.matrix.size1;
  my $ret = gsl_linalg_balance_matrix($A.matrix, $D.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in balance-matrix" if $ret ≠ GSL_SUCCESS;
  return $D;
}

=begin pod

=head1 NAME

Math::Libgsl::LinearAlgebra - An interface to libgsl, the Gnu Scientific Library - Linear Algebra.

=head1 SYNOPSIS

=begin code :lang<perl6>

use Math::Libgsl::LinearAlgebra;

=end code

=head1 DESCRIPTION

Math::Libgsl::LinearAlgebra is an interface to the linear algebra functions of libgsl, the GNU Scientific Library.
This package provides both the low-level interface to the C library (Raw) and a more comfortable interface layer for the Raku programmer.

This module provides functions for Num and Complex data types.

=head2 Num

=head3 LU-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> List)

This function factorizes the matrix A into the LU decomposition PA = LU.
The factorization is done in place, so on output the diagonal and upper triangular (or trapezoidal) part of the input matrix A contain the matrix U. The lower triangular (or trapezoidal) part of the input matrix (excluding the diagonal) contains L. The diagonal elements of L are unity, and are not stored.
The return value is a List: the sign of the permutation and a permutation object, which encodes the permutation matrix P.
In case of error a failure object is returned.

=head3 LU-solve(Math::Libgsl::Matrix $LU where *.matrix.size1 == $LU.matrix.size2, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $LU.matrix.size1 --> Math::Libgsl::Vector)

This function solves the square system Ax = b using the LU decomposition of A into (LU, p) given by the output of LU-decomp.
In case of error a failure object is returned.

=head3 LU-svx(Math::Libgsl::Matrix $LU where *.matrix.size1 == $LU.matrix.size2, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector $x where *.size == $LU.matrix.size1 --> Int)

This function solves the square system Ax = b in-place using the precomputed LU decomposition of A into (LU, p). On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LU-refine(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Matrix $LU where { $LU.matrix.size1 == $LU.matrix.size2 && $A.matrix.size1 == $LU.matrix.size2 }, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $LU.matrix.size1, Math::Libgsl::Vector $x where *.vector.size == $LU.matrix.size1 --> Int)

This function applies an iterative improvement to x, the solution of Ax = b, from the precomputed LU decomposition of A into (LU, p).
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LU-invert(Math::Libgsl::Matrix $LU, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1 --> Math::Libgsl::Matrix)

This function computes the inverse of a matrix A from its LU decomposition (LU, p), returning the matrix inverse.
In case of error a failure object is returned.

=head3 LU-det(Math::Libgsl::Matrix $LU, Int $signum where * ~~ -1|1 --> Num)

This function computes the determinant of a matrix A from its LU decomposition, $LU, and the sign of the permutation, $signum.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LU-lndet(Math::Libgsl::Matrix $LU --> Num)

This function computes the determinant the logarithm of the absolute value of the determinant of a matrix A, ln |det(A)|, from its LU decomposition, $LU.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LU-sgndet(Math::Libgsl::Matrix $LU, Int $signum where * ~~ -1|1 --> Int)

This function computes the sign or phase factor of the determinant of a matrix A, det(A)/|det(A)| from its LU decomposition, $LU.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 QR-decomp(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector)

This function factorizes the M-by-N matrix $A into the QR decomposition A = QR. On output the diagonal and upper triangular part of the input matrix contain the matrix R. The returned vector and the columns of the lower triangular part of the matrix A contain the Householder coefficients and Householder vectors which encode the orthogonal matrix Q.
In case of error a failure object is returned.

=head3 QR-solve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector)

This function solves the square system Ax = b using the QR decomposition of A held in ($QR, $tau) which must have been computed previously with QR-decomp().
In case of error a failure object is returned.

=head3 QR-svx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $QR.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size1 --> Int)

This function solves the square system Ax = b in-place using the QR decomposition of A held in ($QR, $tau) which must have been computed previously by QR-decomp(). On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 QR-lssolve(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> List)

This function finds the least squares solution to the overdetermined system Ax = b where the matrix A has more rows than columns. The least squares solution minimizes the Euclidean norm of the residual, ||Ax − b||. The routine requires as input the QR decomposition of A into ($QR, $tau) given by QR-decomp().
The function returns a List of two Math::Libgsl::Vector objects: the solution x and the residual.
In case of error a failure object is returned.

=head3 QR-QTvec(Math::Libgsl::Matrix $QR, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Vector $v where *.vector.size == $QR.matrix.size1 --> Int)

These function applies the matrix T(Q) encoded in the decomposition (QR, tau) to the vector $v, storing the result T(Q) v in $v. The matrix multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix T(Q).
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 QR-Qvec(Math::Libgsl::Matrix $QR, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Vector $v where *.vector.size == $QR.matrix.size1 --> Int)

This function applies the matrix Q encoded in the decomposition ($QR, $tau) to the vector $v, storing the result Qv in $v. The matrix multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix Q.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 QR-QTmat(Math::Libgsl::Matrix $QR, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Matrix $B where *.matrix.size1 == $QR.matrix.size1 --> Int)

This function applies the matrix T(Q) encoded in the decomposition ($QR, $tau) to the M-by-K matrix $B, storing the result T(Q) B in $B. The matrix multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix T(Q).
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 QR-Rsolve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector)

This function solves the triangular system Rx = b and returns the Math::Libgsl::Vector object $x.
In case of error a failure object is returned.

=head3 QR-Rsvx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size2 --> Int)

This function solves the triangular system Rx = b for x in-place. On input $x should contain the right-hand side b and is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 QR-unpack(Math::Libgsl::Matrix $QR, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2) --> List)

This function unpacks the encoded QR decomposition ($QR, $tau) into the matrices Q and R.
The function returns a List of two Math::Libgsl::Matrix objects: $Q which is M-by-M and $R which is M-by-N.
In case of error a failure object is returned.

=head3 QR-QRsolve(Math::Libgsl::Matrix $Q, Math::Libgsl::Matrix $R where { $R.matrix.size1 == $R.matrix.size2 && $Q.matrix.size1 == $R.matrix.size1 }, Math::Libgsl::Vector $b where *.vector.size == $R.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system Rx = T(Q) b for x. It can be used when the QR decomposition of a matrix is available in unpacked form as ($Q, $R).
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 QR-update(Math::Libgsl::Matrix $Q, Math::Libgsl::Matrix $R where { $Q.matrix.size1 == $R.matrix.size1 && $Q.matrix.size2 == $R.matrix.size1 }, Math::Libgsl::Vector $w where *.vector.size == $R.matrix.size1, Math::Libgsl::Vector $v where *.vector.size == $R.matrix.size2 --> Int)

This function performs a rank-1 update wT(v) of the QR decomposition ($Q, $R). The update is given by Q'R' = Q(R+wT(v)) where the output matrices $Q and $R are also orthogonal and right triangular. Note that $w is destroyed by the update.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 R-solve(Math::Libgsl::Matrix $R where *.matrix.size1 == $R.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $R.matrix.size1 --> Math::Libgsl::Vector)

This function solves the triangular system Rx = b for the N-by-N matrix $R.
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 R-svx(Math::Libgsl::Matrix $R where *.matrix.size1 == $R.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $R.matrix.size2 --> Int)

This function solves the triangular system Rx = b in-place. On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 QRPT-decomp(Math::Libgsl::Matrix $A --> List)

This function factorizes the M-by-N matrix $A into the QRT(P) decomposition A = QRT(P).
On output the diagonal and upper triangular part of the input matrix contain the matrix R.
The function's output is a List of three objects: the Math::Libgsl::Vector $tau, the Math::Libgsl::Permutation $p and the sign of the permutation Int $signum.
In case of error a failure object is returned.

=head3 QRPT-decomp2(Math::Libgsl::Matrix $A --> List)

This function factorizes the matrix $A into the decomposition A = QRT(P) without modifying $A itself.
The function returns a List: the Math::Libgsl::Matrix $Q, the Math::Libgsl::Matrix $R, the Math::Libgsl::Permutation $p, and the sign of the permutation Int $signum.
In case of error a failure object is returned.

=head3 QRPT-solve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector)

This function solves the square system Ax = b using the QRT(P) decomposition of A held in ($QR, $tau, $p) which must have been computed previously by QRPT-decomp.
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 QRPT-svx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size1, Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size2 --> Int)

This function solves the square system Ax = b in-place using the QRT(P) decomposition of A held in ($QR, $tau, $p). On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 QRPT-lssolve(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> List)

This function finds the least squares solution to the overdetermined system Ax = b where the matrix A has more rows than columns and is assumed to have full rank. The least squares solution minimizes the Euclidean norm of the residual, ||b − Ax||. The routine requires as input the QR decomposition of A into ($QR, $tau, $p) given by QRPT-decomp.
The function returns a List of two Math::Libgsl::Vector objects: the solution x and the residual.
In case of error a failure object is returned.

=head3 QRPT-lssolve2(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1, Int $rank where 0 < * ≤ $QR.matrix.size2 --> List)

This function finds the least squares solution to the overdetermined system Ax = b where the matrix A has more rows than columns and has rank given by the input rank. If the user does not know the rank of A, it may be estimated by calling QRPT-rank.
The routine requires as input the QR decomposition of A into ($QR, $tau, $p) given by QRPT-decomp.
The function returns a List of two Math::Libgsl::Vector objects: the solution x and the residual.
In case of error a failure object is returned.

=head3 QRPT-QRsolve(Math::Libgsl::Matrix $Q where *.matrix.size1 == $Q.matrix.size2, Math::Libgsl::Matrix $R where {$R.matrix.size1 == $R.matrix.size2 && $R.matrix.size1 == $Q.matrix.size1}, Math::Libgsl::Permutation $p where *.p.size == $Q.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $Q.matrix.size1 --> Math::Libgsl::Vector)

This function solves the square system RT(P) x = T(Q) b for x. It can be used when the QR decomposition of a matrix is available in unpacked form as ($Q, $R).
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 QRPT-update(Math::Libgsl::Matrix $Q, Math::Libgsl::Matrix $R where { $Q.matrix.size1 == $R.matrix.size1 && $Q.matrix.size2 == $R.matrix.size1 }, Math::Libgsl::Permutation $p where *.p.size == $R.matrix.size1, Math::Libgsl::Vector $w where *.vector.size == $R.matrix.size1, Math::Libgsl::Vector $v where *.vector.size == $R.matrix.size2 --> Int)

This function performs a rank-1 update wT(v) of the QRT(P) decomposition ($Q, $R, $p). The update is given by Q' R' = Q(R + wT(v) P) where the output matrices Q' and R' are also orthogonal and right triangular. Note that $w is destroyed by the update.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 QRPT-Rsolve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector)

This function solves the triangular system RT(P) x = b for the N-by-N matrix R contained in $QR.
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 QRPT-Rsvx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size2 --> Int)

This function solves the triangular system RT(P) x = b in-place for the N-by-N matrix R contained in $QR. On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 QRPT-rank(Math::Libgsl::Matrix $QR, Num() $tolerance --> Int)

This function returns the rank of the triangular matrix R contained in $QR.

=head3 QRPT-rcond(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2 --> Num)

This function returns the reciprocal condition number (using the 1-norm) of the R factor, stored in the upper triangle of $QR.
In case of error a failure object is returned.

=head3 LQ-decomp(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector)

This function factorizes the M-by-N matrix $A into the LQ decomposition A = LQ. On output the diagonal and lower trapezoidal part of the input matrix contain the matrix L.
This function returns the Math::Libgsl::Vector $tau.
The vector $tau and the elements above the diagonal of the matrix $A contain the Householder coefficients and Householder vectors which encode the orthogonal matrix Q.
In case of error a failure object is returned.

=head3 LQ-solve-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2), Math::Libgsl::Vector $b where *.vector.size == $LQ.matrix.size2 --> Math::Libgsl::Vector)

This function finds the solution to the system Ax = b. The routine requires as input the LQ decomposition of A into ($LQ, $tau) given by LQ-decomp.
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 LQ-svx-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2), --> Math::Libgsl::Vector)

This function finds the solution to the system Ax = b. The routine requires as input the LQ decomposition of A into ($LQ, $tau) given by LQ-decomp. On input $x should contain the right-hand side b, which is replaced by the solution on output.
In case of error a failure object is returned.

=head3 LQ-lssolve-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 ≥ $LQ.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2), Math::Libgsl::Vector $b where *.vector.size == $LQ.matrix.size1 --> List)

This function finds the minimum norm least squares solution to the underdetermined system Ax = b, where the M-by-N matrix A has M ≤ N. The routine requires as input the LQ decomposition of A into ($LQ, $tau) given by LQ-decomp.
The function returns a List of two Math::Libgsl::Vector objects: the solution x and the residual.
In case of error a failure object is returned.

=head3 LQ-Lsolve-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $LQ.matrix.size1 --> Math::Libgsl::Vector)

The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 LQ-Lsvx-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2 --> Math::Libgsl::Vector)

On input $x should contain the right-hand side b, which is replaced by the solution on output.
In case of error a failure object is returned.

=head3 L-solve-T(Math::Libgsl::Matrix $L where *.matrix.size1 == $L.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $L.matrix.size2 --> Math::Libgsl::Vector)

The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 LQ-vecQ(Math::Libgsl::Matrix $LQ, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2), Math::Libgsl::Vector $v where *.vector.size == $LQ.matrix.size1 --> Int)

This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LQ-vecQT(Math::Libgsl::Matrix $LQ, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2), Math::Libgsl::Vector $v where *.vector.size == $LQ.matrix.size1 --> Int)

This function applies T(Q) to the vector v, storing the result T(Q) v in $v on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LQ-unpack(Math::Libgsl::Matrix $LQ, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2) --> List)

This function unpacks the encoded LQ decomposition ($LQ, $tau).
The function outputs a List: the Math::Libgsl::Matrix $Q and the Math::Libgsl::Matrix $L.
In case of error a failure object is returned.

=head3 LQ-update(Math::Libgsl::Matrix $Q, Math::Libgsl::Matrix $L where { $L.matrix.size2 == $Q.matrix.size1 && $L.matrix.size2 == $Q.matrix.size2 }, Math::Libgsl::Vector $v where *.vector.size == $L.matrix.size1, Math::Libgsl::Vector $w where *.vector.size == $L.matrix.size2 --> Int)

This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LQ-LQsolve(Math::Libgsl::Matrix $Q, Math::Libgsl::Matrix $L where { $L.matrix.size1 == $L.matrix.size2 && $Q.matrix.size1 == $L.matrix.size2 }, Math::Libgsl::Vector $b where *.vector.size == $L.matrix.size2 --> Math::Libgsl::Vector)

In case of error a failure object is returned.

=head3 COD-decomp(Math::Libgsl::Matrix $A --> List)
=head3 COD-decomp-e(Math::Libgsl::Matrix $A, Num() $tolerance --> List)

These functions factor the M-by-N matrix $A into the decomposition A = QRZT(P).
On output the matrix R₁₁ is stored in the upper rank-by-rank block of $A.
The matrices Q and Z are encoded in packed storage in $A on output.
This function outputs a List: two Math::Libgsl::Vector objects which contain the Householder scalars corresponding to the matrices Q and Z respectively, a Math::Libgsl::Permutation object which contain the permutation matrix P, and the rank of $A.
In case of error a failure object is returned.

=head3 COD-lssolve(Math::Libgsl::Matrix $QRZT where *.matrix.size1 ≥ $QRZT.matrix.size2, Math::Libgsl::Vector $tau-Q where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QRZT.matrix.size2, Int $rank where * ≤ min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Vector $b where *.vector.size == $QRZT.matrix.size1 --> List)

This function finds the unique minimum norm least squares solution to the overdetermined system Ax = b where the matrix $A has more rows than columns. The least squares solution minimizes the Euclidean norm of the residual, ||b − Ax|| as well as the norm of the solution ||x||. The routine requires as input the QRZT decomposition of $A into ($QRZT, $tau_Q, $tau_Z, $p, $rank) given by COD-decomp.
The function outputs a List: a Math::Libgsl::Vector object which is the solution x, and a Math::Libgsl::Vector object which stores the residual b − Ax.
In case of error a failure object is returned.

=head3 COD-lssolve2(Math::Libgsl::Matrix $QRZT where *.matrix.size1 ≥ $QRZT.matrix.size2, Math::Libgsl::Vector $tau-Q where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QRZT.matrix.size2, Int $rank where * ≤ min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Vector $b where *.vector.size == $QRZT.matrix.size1, Num() $lambda --> List)

This function finds the solution to the regularized least squares problem in Tikhonov standard form, minₓ||b − Ax||² + λ²||x||². The routine requires as input the QRZT decomposition of A into ($QRZT, $tau_Q, $tau_Z, $p, $rank) given by COD-decomp. The parameter λ is supplied in $lambda.
The function outputs a List: a Math::Libgsl::Vector object which is the solution x, and a Math::Libgsl::Vector object which stores the residual b − Ax.
In case of error a failure object is returned.

=head3 COD-unpack(Math::Libgsl::Matrix $QRZT, Math::Libgsl::Vector $tau-Q where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Int $rank where * ≤ min($QRZT.matrix.size1, $QRZT.matrix.size2) --> List)

This function unpacks the encoded QRZT decomposition ($QRZT, $tau_Q, $tau_Z, $rank).
The function returns a List of three Math::Libgsl::Matrix objects: $Q, $R, $Z.
In case of error a failure object is returned.

=head3 COD-matZ(Math::Libgsl::Matrix $QRZT, Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Matrix $A where *.matrix.size2 == $QRZT.matrix.size2, Int $rank --> Int)

This function multiplies the input matrix $A on the right by Z, A’ = AZ using the encoded QRZT decomposition ($QRZT, $tau_Z, $rank).
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 SV-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List)

This function factorizes the M-by-N matrix $A into the singular value decomposition A = UST(V) for M ≥ N. On output the matrix $A is replaced by U.
The function returns a List: a Math::Libgsl::Matrix object which contains the elements of V in untransposed form, and a Math::Libgsl::Vector object which contains the diagonal elements of the singular value matrix S. The singular values are non-negative and form a non-increasing sequence from S₁ to Sₙ.
In case of error a failure object is returned.

=head3 SV-decomp-mod(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List)

This function computes the SVD using the modified Golub-Reinsch algorithm, which is faster for M ≫ N.
The function returns a List: a Math::Libgsl::Matrix object which contains the elements of V in untransposed form, and a Math::Libgsl::Vector object which contains the diagonal elements of the singular value matrix S. The singular values are non-negative and form a non-increasing sequence from S₁ to Sₙ.
In case of error a failure object is returned.

=head3 SV-decomp-jacobi(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List)

This function computes the SVD of the M-by-N matrix A using one-sided Jacobi orthogonalization for M ≥ N.
The function returns a List: a Math::Libgsl::Matrix object which contains the elements of V in untransposed form, and a Math::Libgsl::Vector object which contains the diagonal elements of the singular value matrix S. The singular values are non-negative and form a non-increasing sequence from S₁ to Sₙ.
In case of error a failure object is returned.

=head3 SV-solve(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2, Math::Libgsl::Matrix $V where { $V.matrix.size1 == $V.matrix.size2 && $V.matrix.size1 == $A.matrix.size2 }, Math::Libgsl::Vector $S where *.vector.size == $A.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system Ax = b using the singular value decomposition (U, S, V) of A which must have been computed previously with COD-decomp.
Only non-zero singular values are used in computing the solution. The parts of the solution corresponding to singular values of zero are ignored. Other singular values can be edited out by setting them to zero before calling this function.
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 SV-leverage(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector)

This function computes the statistical leverage values hᵢ of a matrix A using its singular value decomposition (U, S, V) previously computed with COD-decomp.
The function returns a Math::Libgsl::Vector object which stores the diagonal values of the matrix A(T(A)A)⁻¹T(A) and depend only on the matrix U which is the input to this function.
In case of error a failure object is returned.

=head3 cholesky-decomp1(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int)

This function factorizes the symmetric, positive-definite square matrix $A into the Cholesky decomposition A = LT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used (the upper triangular part is ignored). On output the diagonal and lower triangular part of the input matrix $A contain the matrix L, while the upper triangular part contains the original matrix. If the matrix is not positive-definite then the decomposition will fail, returning the error code GSL_EDOM.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-solve(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector)

These functions solve the system Ax = b using the Cholesky decomposition of $A which must have been previously computed by cholesky-decomp1.
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 cholesky-svx(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $A.matrix.size1 --> Int)

This function solves the system Ax = b in-place using the Cholesky decomposition of $A which must have been previously computed by cholesky-decomp1.
On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-solve-mat(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Matrix $B where *.matrix.size1 == $A.matrix.size1 --> Math::Libgsl::Matrix)

In case of error a failure object is returned.

=head3 cholesky-svx-mat(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Matrix $X where *.matrix.size1 == $A.matrix.size1 --> Int)

This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-invert(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int)

These functions compute the inverse of a matrix from its Cholesky decomposition $A, which must have been previously computed by cholesky-decomp1.
On output, the inverse is stored in-place in $A.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-decomp2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector)

This function calculates a diagonal scaling transformation S for the symmetric, positive-definite square matrix $A, and then computes the Cholesky decomposition SAS = LT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used (the upper triangular part is ignored). On output the diagonal and lower triangular part of the input matrix $A contain the matrix L, while the upper triangular part of the input matrix is overwritten with T(L) (the diagonal terms being identical for both L and L T). If the matrix is not positive-definite then the decomposition will fail, returning the error code GSL_EDOM.
The function returns a Math::Libgsl::Vector object which stores the diagonal scale factors.
In case of error a failure object is returned.

=head3 cholesky-solve2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $s where *.vector.size == $A.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system (SAS)(S⁻¹ x) = Sb using the Cholesky decomposition of SAS held in the matrix $A which must have been previously computed by cholesky-decomp2.
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 cholesky-svx2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $s where *.size == $A.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $A.matrix.size2 --> Int)

This function solves the system (SAS)(S⁻¹ x) = Sb using the Cholesky decomposition of SAS held in the matrix $A which must have been previously computed by cholesky-decomp2.
On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-decomp-unit(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector)

The function returns the Math::Libgsl::Vector $d.
In case of error a failure object is returned.

=head3 cholesky-scale(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector)

This function calculates a diagonal scaling transformation of the symmetric, positive definite matrix $A, such that SAS has a condition number within a factor of N of the matrix of smallest possible condition number over all possible diagonal scalings.
The function outputs a Math::Libgsl::Vector object which contains the scale factors.
In case of error a failure object is returned.

=head3 cholesky-scale-apply(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $s where *.vector.size == $A.matrix.size2 --> Int)

This function applies the scaling transformation S to the matrix $A. On output, $A is replaced by SAS.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-rcond(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num)

This function estimates the reciprocal condition number (using the 1-norm) of the symmetric positive definite matrix A, using its Cholesky decomposition provided in $A.
The function returns a Math::Libgsl::Vector object which stores the reciprocal condition number estimate, defined as 1/(||A||₁ · ||A⁻¹||₁).
In case of error a failure object is returned.

=head3 pcholesky-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Permutation)

This function factors the symmetric, positive-definite square matrix $A into the Pivoted Cholesky decomposition PAT(P) = LDT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used to construct the factorization. On output the diagonal of the input matrix $A stores the diagonal elements of D, and the lower triangular portion of $A contains the matrix L. Since L has ones on its diagonal these do not need to be explicitely stored. The upper triangular portion of $A is unmodified.
The function returns the Math::Libgsl::Permutation object which stores the permutation matrix P.
In case of error a failure object is returned.

=head3 pcholesky-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system Ax = b using the Pivoted Cholesky decomposition of A held in the matrix $LDLT and permutation $p which must have been previously computed by pcholesky-decomp.
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 pcholesky-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int)

This function solves the system Ax = b using the Pivoted Cholesky decomposition of A held in the matrix $LDLT and permutation $p which must have been previously computed by pcholesky-decomp.
On input, $x contains the right hand side vector b which is replaced by the solution vector on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 pcholesky-decomp2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> List)

This function computes the pivoted Cholesky factorization of the matrix SAS, where the input matrix $A is symmetric and positive definite, and the diagonal scaling matrix S is computed to reduce the condition number of $A as much as possible.
On input, the values from the diagonal and lower-triangular part of the matrix $A are used to construct the factorization. On output the diagonal of the input matrix $A stores the diagonal elements of D, and the lower triangular portion of $A contains the matrix L. Since L has ones on its diagonal these do not need to be explicitely stored. The upper triangular portion of $A is unmodified.
The function returns a List: the permutation matrix P is stored in a Math::Libgsl::Permutation object, the diagonal scaling transformation is stored in a Math::Libgsl::Vector object.
In case of error a failure object is returned.

=head3 pcholesky-solve2(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $s where *.vector.size == $LDLT.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system (SAS)(S⁻¹x) = Sb using the Pivoted Cholesky decomposition of SAS held in the matrix $LDLT, permutation $p, and vector $S, which must have been previously computed by pcholesky-decomp2.
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 pcholesky-svx2(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $s where *.vector.size == $LDLT.matrix.size1, Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int)

This function solves the system (SAS)(S⁻¹x) = Sb using the Pivoted Cholesky decomposition of SAS held in the matrix $LDLT, permutation $p, and vector $S, which must have been previously computed by pcholesky-decomp2.
On input, $x contains the right hand side vector b which is replaced by the solution vector on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 pcholesky-invert(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1 --> Math::Libgsl::Matrix)

This function computes the inverse of the matrix A, using the Pivoted Cholesky decomposition stored in $LDLT and $p.
The function returns the Math::Libgsl::Matrix A⁻¹.
In case of error a failure object is returned.

=head3 pcholesky-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1 --> Num)

This function estimates the reciprocal condition number (using the 1-norm) of the symmetric positive definite matrix A, using its pivoted Cholesky decomposition provided in $LDLT.
The function returns the reciprocal condition number estimate, defined as 1/(||A||₁ · ||A⁻¹||₁).
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 mcholesky-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Bool :$perturbation = True --> List)

This function factors the symmetric, indefinite square matrix $A into the Modified Cholesky decomposition P(A + E)T(P) = LDT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used to construct the factorization. On output the diagonal of the input matrix $A stores the diagonal elements of D, and the lower triangular portion of $A contains the matrix L. Since L has ones on its diagonal these do not need to be explicitely stored. The upper triangular portion of $A is unmodified.
The function returns a List: the permutation matrix P, stored in a Math::Libgsl::Permutation object and the diagonal perturbation matrix, stored in a Math::Libgsl::Vector object.
In case of error a failure object is returned.

=head3 mcholesky-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1 --> Math::Libgsl::Vector)

This function solves the perturbed system (A + E)x = b using the Cholesky decomposition of A + E held in the matrix $LDLT and permutation $p which must have been previously computed by mcholesky-decomp.
The function returns the Math::Libgsl::Vector $x.
In case of error a failure object is returned.

=head3 mcholesky-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int)

This function solves the perturbed system (A + E)x = b using the Cholesky decomposition of A + E held in the matrix $LDLT and permutation $p which must have been previously computed by mcholesky-decomp.
On input, $x contains the right hand side vector b which is replaced by the solution vector on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 mcholesky-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1 --> Num)

This function estimates the reciprocal condition number (using the 1-norm) of the perturbed matrix A+E, using its pivoted Cholesky decomposition provided in $LDLT.
The function returns the reciprocal condition number estimate, defined as 1/(||A + E||₁ · ||(A + E)⁻¹||₁).
In case of error a failure object is returned.

=head3 mcholesky-invert(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1 --> Math::Libgsl::Matrix)

The function returns the inverted matrix.
In case of error a failure object is returned.

=head3 ldlt-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int)

This function factorizes the symmetric, non-singular square matrix $A into the decomposition A = LDT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used. The upper triangle of $A is used as temporary workspace. On output the diagonal of $A contains the matrix D and the lower triangle of $A contains the unit lower triangular matrix L. The matrix 1-norm, ||A||₁ is stored in the upper right corner on output
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.
This function is available only in the C library starting from v2.6.

=head3 ldlt-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system Ax = b using the LDT(L) decomposition of A held in the matrix $LDLT which must have been previously computed by ldlt-decomp.
In case of error a failure object is returned.
This function is available only in the C library starting from v2.6.

=head3 ldlt-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int)

This function solves the system Ax = b using the LDT(L) decomposition of A held in the matrix $LDLT which must have been previously computed by ldlt-decomp.
On input $x should contain the right-hand side b, which is replaced by the solution on output.
In case of error a failure object is returned.
This function is available only in the C library starting from v2.6.

=head3 ldlt-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2 --> Num)

This function estimates the reciprocal condition number (using the 1-norm) of the symmetric nonsingular matrix A, using its LDT(L) decomposition provided in $LDLT.
The function returns the reciprocal condition number estimate, defined as 1/(||A + E||₁ · ||(A + E)⁻¹||₁).
In case of error a failure object is returned.
This function is available only in the C library starting from v2.6.

=head3 symmtd-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector)

This function factorizes the symmetric square matrix $A into the symmetric tridiagonal decomposition QTT(Q). On output the diagonal and subdiagonal part of the input matrix $A contain the tridiagonal matrix T. The remaining lower triangular part of the input matrix contains the Householder vectors which, together with the Householder coefficients tau returned as a Math::Libgsl::Vector object, encode the orthogonal matrix Q. The upper triangular part of $A is not referenced.
In case of error a failure object is returned.

=head3 symmtd-unpack(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $A.matrix.size1 - 1 --> List)

This function unpacks the encoded symmetric tridiagonal decomposition ($A, $tau) obtained from symmtd-decomp.
The function returns a List: the orthogonal matrix $Q as a Math::Libgsl::Matrix, the vector of diagonal elements $diag as a Math::Libgsl::Vector, and the vector of subdiagonal elements $subdiag as a Math::Libgsl::Vector.
In case of error a failure object is returned.

=head3 symmtd-unpack-T(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> List)

This function unpacks the diagonal and subdiagonal of the encoded symmetric tridiagonal decomposition ($A, $tau) obtained from symmtd-decomp.
The function returns a List of two Math::Libgsl::Vector: $diag and $subdiag.
In case of error a failure object is returned.

=head3 hessenberg-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector)

This function computes the Hessenberg decomposition of the matrix $A by applying the similarity transformation H = T(U)AU. On output, H is stored in the upper portion of $A. The information required to construct the matrix U is stored in the lower triangular portion of $A. U is a product of N−2 Householder matrices. The Householder vectors are stored in the lower portion of $A (below the subdiagonal) and the Householder coefficients are returned as a Math::Libgsl::Vector.
In case of error a failure object is returned.

=head3 hessenberg-unpack(Math::Libgsl::Matrix $H where *.matrix.size1 == $H.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $H.matrix.size1 --> Math::Libgsl::Matrix)

This function constructs the orthogonal matrix U from the information stored in the Hessenberg matrix $H along with the vector $tau. $H and $tau are outputs from hessenberg-decomp.
The function returns the Math::Libgsl::Matrix which stores $U.
In case of error a failure object is returned.

=head3 hessenberg-unpack-accum(Math::Libgsl::Matrix $H where *.matrix.size1 == $H.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $H.matrix.size1 --> Math::Libgsl::Matrix)

This function is similar to hessenberg-unpack, except it accumulates the matrix U into V, so that V′ = VU.
The function returns the Math::Libgsl::Matrix which stores $V.
In case of error a failure object is returned.

=head3 hessenberg-set-zero(Math::Libgsl::Matrix $H where *.matrix.size1 == $H.matrix.size2 --> Int)

This function sets the lower triangular portion of $H, below the subdiagonal, to zero. It is useful for clearing out the Householder vectors after calling hessenberg-decomp.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 hesstri-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Matrix $B where { $B.matrix.size1 == $A.matrix.size1 && $B.matrix.size2 == $A.matrix.size2}, Bool :$similarity = True --> List)

This function computes the Hessenberg-Triangular decomposition of the matrix pair ($A, $B). On output, H is stored in $A, and R is stored in $B.
If the $similarity Bool parameter is True (default), then the similarity transformations are returned as a List of Math::Libgsl::Matrix objects.
In case of error a failure object is returned.

=head3 bidiag-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List)

This function factorizes the M-by-N matrix $A into bidiagonal form UBT(V). The diagonal and superdiagonal of the matrix B are stored in the diagonal and superdiagonal of $A. The orthogonal matrices U and V are stored as compressed Householder vectors in the remaining elements of $A.
This function returns a List of Math::Libgsl::Vector objects: the two Householder coefficients $tau_U and $tau_V.
In case of error a failure object is returned.

=head3 bidiag-unpack(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2, Math::Libgsl::Vector $tau_U where *.vector.size == $A.matrix.size2, Math::Libgsl::Vector $tau_V where *.vector.size == $A.matrix.size2 - 1 --> List)

This function unpacks the bidiagonal decomposition of $A produced by bidiag-decomp, ($A, $tau_U, $tau_V) into the separate orthogonal matrices U, V and the diagonal vector diag and superdiagonal superdiag. Note that U is stored as a compact M-by-N orthogonal matrix satisfying T(U)U = I for efficiency.
The function returns a List of four objects: the Math::Libgsl::Matrix $U and $V, and the Math::Libgsl::Vector $diag and $sdiag.
In case of error a failure object is returned.

=head3 bidiag-unpack2(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2, Math::Libgsl::Vector $tau_U where *.vector.size == $A.matrix.size2, Math::Libgsl::Vector $tau_V where *.vector.size == $A.matrix.size2 - 1 --> Math::Libgsl::Matrix)

This function unpacks the bidiagonal decomposition of $A produced by bidiag-decomp, ($A, $tau_U, $tau_V) into the separate orthogonal matrices U, V and the diagonal vector diag and superdiagonal superdiag. The matrix U is stored in-place in $A.
The function returns the Math::Libgsl::Matrix object $V.
In case of error a failure object is returned.

=head3 bidiag-unpack-B(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List)

This function unpacks the diagonal and superdiagonal of the bidiagonal decomposition of $A from bidiag-decomp.
The function returns a List of two objects: the Math::Libgsl::Vector $diag and $sdiag.
In case of error a failure object is returned.

=head3 givens(Num() $a, Num() $b --> List)

This function computes c = cos θ and s = sin θ so that the Givens matrix G(θ) acting on the vector (a, b) produces (r, 0), with r = √ a² + b².
The function returns a List of Num: the c and s elements of the Givens matrix.

=head3 givens-gv(Math::Libgsl::Vector $v, Int, $i, Int $j, Num() $c, Num() $s)

This function applies the Givens rotation defined by c = cos θ and s = sin θ to the i and j elements of v. On output, (v(i), v(j)) ← G(θ)(v(i), v(j)).
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 householder-transform(Math::Libgsl::Vector $w --> Num)

This function prepares a Householder transformation H = I − τvT(v) which can be used to zero all the elements of the input vector $w except the first. On output the Householder vector v is stored in $w and the scalar τ is returned. The householder vector v is normalized so that v[0] = 1, however this 1 is not stored in the output vector. Instead, $w[0] is set to the first element of the transformed vector, so that if u = Hw, w[0] = u[0] on output and the remainder of u is zero.
This function returns a Num: $τ.

=head3 householder-hm(Num() $tau, Math::Libgsl::Vector $v, Math::Libgsl::Matrix $A --> Int)

This function applies the Householder matrix H defined by the scalar $tau and the vector $v to the left-hand side of the matrix $A. On output the result HA is stored in $A.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 householder-mh(Num() $tau, Math::Libgsl::Vector $v, Math::Libgsl::Matrix $A --> Int)

This function applies the Householder matrix H defined by the scalar $tau and the vector $v to the right-hand side of the matrix $A. On output the result AH is stored in $A.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 householder-hv(Num() $tau, Math::Libgsl::Vector $v, Math::Libgsl::Vector $w --> Int)

This function applies the Householder transformation H defined by the scalar $tau and the vector $v to the vector $w. On output the result Hw is stored in w.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 HH-solve(Math::Libgsl::Matrix $A where *.matrix.size1 ≤ $A.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system Ax = b directly using Householder transformations.
$b is not modified. The matrix $A is destroyed by the Householder transformations.
The function returns a Math::Libgsl::Vector object: the solution $x.
In case of error a failure object is returned.

=head3 HH-svx(Math::Libgsl::Matrix $A where *.matrix.size1 ≤ $A.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $A.matrix.size2 --> Int)

This function solves the system Ax = b in-place using Householder transformations. On input $x should contain the right-hand side b, which is replaced by the solution on output.
The matrix $A is destroyed by the Householder transformations.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tridiag-solve(Math::Libgsl::Vector $diag, Math::Libgsl::Vector $abovediag where *.vector.size == $diag.vector.size - 1, Math::Libgsl::Vector $belowdiag where *.vector.size == $diag.vector.size - 1, Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size --> Math::Libgsl::Vector)

This function solves the general N-by-N system Ax = b where A is tridiagonal (N ≥ 2). The super-diagonal and sub-diagonal vectors must be one element shorter than the diagonal vector diag.
The function returns the Math::Libgsl::Vector solution $x.
In case of error a failure object is returned.

=head3 tridiag-symm-solve(Math::Libgsl::Vector $diag, Math::Libgsl::Vector $offdiag where *.vector.size == $diag.vector.size - 1, Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size --> Math::Libgsl::Vector)

This function solves the general N-by-N system Ax = b where A is symmetric tridiagonal (N ≥ 2). The off-diagonal vector $offdiag must be one element shorter than the diagonal vector $diag.
The function returns the Math::Libgsl::Vector solution $x.
In case of error a failure object is returned.

=head3 tridiag-cyc-solve(Math::Libgsl::Vector $diag where *.vector.size ≥ 3, Math::Libgsl::Vector $abovediag where *.vector.size == $diag.vector.size, Math::Libgsl::Vector $belowdiag where *.vector.size == $diag.vector.size, Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size --> Math::Libgsl::Vector)

This function solves the general N-by-N system Ax = b where A is cyclic tridiagonal (N ≥ 3). The cyclic super-diagonal and sub-diagonal vectors $abovediag and $belowdiag must have the same number of elements as the diagonal vector $diag.
The function returns the Math::Libgsl::Vector solution $x.
In case of error a failure object is returned.

=head3 tridiag-symm-cyc-solve(Math::Libgsl::Vector $diag where *.vector.size ≥ 3, Math::Libgsl::Vector $offdiag where *.vector.size == $diag.vector.size, Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size --> Math::Libgsl::Vector)

This function solves the general N-by-N system Ax = b where A is symmetric cyclic tridiagonal (N ≥ 3). The cyclic off-diagonal vector $offdiag must have the same number of elements as the diagonal vector $diag.
The function returns the Math::Libgsl::Vector solution $x.
In case of error a failure object is returned.

=head3 tri-upper-rcond(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num)

This function estimates the reciprocal condition number.
In case of error a failure object is returned.

=head3 tri-lower-rcond(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num)

This function estimates the reciprocal condition number.
In case of error a failure object is returned.

=head3 tri-upper-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int)

This function computes the in-place inverse of the triangular matrix $T, stored in the upper triangle.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tri-lower-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int)

This function computes the in-place inverse of the triangular matrix $T, stored in the lower triangle.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tri-upper-unit-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int)

This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tri-lower-unit-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int)

This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tri-invert(Int $Uplo, Int $Diag, Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int)

This function is available only from the C library v2.6.
This function computes the in-place inverse of the triangular matrix $T, stored in the lower triangle when $Uplo = CblasLower and upper triangle when $Uplo = CblasUpper. The parameter $Diag = CblasUnit, CblasNonUnit specifies whether the matrix is unit triangular.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tri-LTL(Math::Libgsl::Matrix $L where *.matrix.size1 == $L.matrix.size2 --> Int)

This function is available only from the C library v2.6.
This function computes the product LT(L) in-place and stores it in the lower triangle of $L on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tri-UL(Math::Libgsl::Matrix $LU where *.matrix.size1 == $LU.matrix.size2 --> Int)

This function is available only from the C library v2.6.
This function compute the product $UL where U is upper triangular and L is unit lower triangular, stored in $LU, as computed by LU-decomp.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tri-rcond(Int $Uplo, Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num)

This function is available only from the C library v2.6.
This function estimates the 1-norm reciprocal condition number of the triangular matrix $A, using the lower triangle when $Uplo is CblasLower and upper triangle when $Uplo is CblasUpper.
The function returns the reciprocal condition number 1/(||A||₁||A⁻¹||₁).
In case of error a failure object is returned.

=head3 cholesky-band-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int)

This function is available only from the C library v2.6.
This function factorizes the symmetric, positive-definite square matrix $A into the Cholesky decomposition A = LT(L). The input matrix $A is given in symmetric banded format, and has dimensions N-by-(p + 1), where p is the lower bandwidth of the matrix. On output, the entries of $A are replaced by the entries of the matrix L in the same format. In addition, the lower right element of $A is used to store the matrix 1-norm, used later by cholesky-band-rcond() to calculate the reciprocal condition number.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-band-solve(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $LLT.matrix.size1 --> Math::Libgsl::Vector)

This function is available only from the C library v2.6.
This function solves the symmetric banded system Ax = b using the Cholesky decomposition of A held in the matrix $LLT which must have been previously computed by cholesky-band-decomp.
The function returns the Math::Libgsl::Vector solution $x.
In case of error a failure object is returned.

=head3 cholesky-band-svx(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $LLT.matrix.size1 --> Int)

This function is available only from the C library v2.6.
This function solves the symmetric banded system Ax = b using the Cholesky decomposition of A held in the matrix $LLT which must have been previously computed by cholesky-band-decomp.
On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-band-invert(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2 --> Math::Libgsl::Matrix)

This function is available only from the C library v2.6.
This function computes the inverse of a symmetric banded matrix from its Cholesky decomposition $LLT, which must have been previously computed by cholesky-band-decomp.
The function returns the inverse matrix as a Math::Libgsl::Matrix solution object.
In case of error a failure object is returned.

=head3 cholesky-band-unpack(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2 --> Math::Libgsl::Matrix)

This function is available only from the C library v2.6.
This function unpacks the lower triangular Cholesky factor from $LLT, which returns as a Math::Libgsl::Matrix object.
In case of error a failure object is returned.

=head3 cholesky-band-rcond(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2 --> Num)

This function is available only from the C library v2.6.
This function estimates the reciprocal condition number (using the 1-norm) of the symmetric banded positive definite matrix A, using its Cholesky decomposition provided in $LLT. The reciprocal condition number estimate, defined as 1/(||A||₁ · ||A⁻¹||₁), is returned.

=head3 ldlt-band-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int)

This function is available only from the C library v2.6.
This function factorizes the symmetric, non-singular square matrix $A into the decomposition A = LDT(L). The input matrix $A is given in symmetric banded format, and has dimensions N-by-(p + 1), where p is the lower bandwidth of the matrix. On output, the entries of $A are replaced by the entries of the matrices D and L in the same format.
If the matrix is singular then the decomposition will fail, returning the error code GSL_EDOM .
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 ldlt-band-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1 --> Math::Libgsl::Vector)

This function is available only from the C library v2.6.
This function solves the symmetric banded system Ax = b using the LDT(L) decomposition of A held in the matrix $LDLT which must have been previously computed by ldlt-band-decomp.
The function returns the Math::Libgsl::Vector solution $x.
In case of error a failure object is returned.

=head3 ldlt-band-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int)

This function is available only from the C library v2.6.
This function solves the symmetric banded system Ax = b using the LDT(L) decomposition of A held in the matrix $LDLT which must have been previously computed by ldlt-band-decomp.
On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 ldlt-band-unpack(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2 --> List)

This function is available only from the C library v2.6.
This function unpacks the unit lower triangular factor L from $LDLT.
This function returns a Math::Libgsl::Matrix object, which holds the lower triangular portion of the matrix L, and a Math::Libgsl::Vector object, which contains the diagonal matrix D.
In case of error a failure object is returned.

=head3 ldlt-band-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2 --> Num)

This function is available only from the C library v2.6.
This function estimates the reciprocal condition number (using the 1-norm) of the symmetric banded nonsingular matrix A, using its LDT(L) decomposition provided in $LDLT. The reciprocal condition number estimate, defined as 1/(||A||₁ · ||A⁻¹||₁), is returned.

=head3 balance-matrix(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector)

This function replaces the matrix $A with its balanced counterpart and returns the diagonal elements of the similarity transformation as a Math::Libgsl::Vector object.
In case of error a failure object is returned.

=head2 Complex

=head3 LU-cdecomp(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> List)

This function factorizes the matrix A into the LU decomposition PA = LU.
The factorization is done in place, so on output the diagonal and upper triangular (or trapezoidal) part of the input matrix A contain the matrix U. The lower triangular (or trapezoidal) part of the input matrix (excluding the diagonal) contains L. The diagonal elements of L are unity, and are not stored.
The return value is a List: the sign of the permutation and a permutation object, which encodes the permutation matrix P.
In case of error a failure object is returned.

=head3 LU-csolve(Math::Libgsl::Matrix::Complex64 $LU where *.matrix.size1 == $LU.matrix.size2, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector::Complex64 $b where *.vector.size == $LU.matrix.size1 --> Math::Libgsl::Vector::Complex64)

This function solves the square system Ax = b using the LU decomposition of A into (LU, p) given by the output of LU-decomp.
In case of error a failure object is returned.

=head3 LU-csvx(Math::Libgsl::Matrix::Complex64 $LU where *.matrix.size1 == $LU.matrix.size2, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector::Complex64 $x where *.vector.size == $LU.matrix.size1 --> Int)

This function solves the square system Ax = b in-place using the precomputed LU decomposition of A into (LU, p). On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LU-crefine(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Matrix::Complex64 $LU where { $LU.matrix.size1 == $LU.matrix.size2 && $A.matrix.size1 == $LU.matrix.size2 }, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector::Complex64 $b where *.vector.size == $LU.matrix.size1, Math::Libgsl::Vector::Complex64 $x where *.vector.size == $LU.matrix.size1 --> Int)

This function applies an iterative improvement to x, the solution of Ax = b, from the precomputed LU decomposition of A into (LU, p).
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LU-cinvert(Math::Libgsl::Matrix::Complex64 $LU, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1 --> Math::Libgsl::Matrix::Complex64)

This function computes the inverse of a matrix A from its LU decomposition (LU, p), returning the matrix inverse.
In case of error a failure object is returned.

=head3 LU-cdet(Math::Libgsl::Matrix::Complex64 $LU, Int $signum where * ~~ -1|1 --> Complex)

This function computes the determinant of a matrix A from its LU decomposition, $LU, and the sign of the permutation, $signum.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LU-clndet(Math::Libgsl::Matrix::Complex64 $LU --> Num)

This function computes the determinant the logarithm of the absolute value of the determinant of a matrix A, ln |det(A)|, from its LU decomposition, $LU.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 LU-csgndet(Math::Libgsl::Matrix::Complex64 $LU, Int $signum where * ~~ -1|1 --> Complex)

This function computes the sign or phase factor of the determinant of a matrix A, det(A)/|det(A)| from its LU decomposition, $LU.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-cdecomp(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> Int)

This function factorizes the symmetric, positive-definite square matrix $A into the Cholesky decomposition A = LT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used (the upper triangular part is ignored). On output the diagonal and lower triangular part of the input matrix $A contain the matrix L, while the upper triangular part contains the original matrix. If the matrix is not positive-definite then the decomposition will fail, returning the error code GSL_EDOM.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-csolve(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector::Complex64 $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector::Complex64)

These functions solve the system Ax = b using the Cholesky decomposition of $A which must have been previously computed by cholesky-cdecomp.
The function returns the Math::Libgsl::Vector::Complex64 $x.
In case of error a failure object is returned.

=head3 cholesky-csvx(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector::Complex64 $x where *.vector.size == $A.matrix.size1 --> Int)

This function solves the system Ax = b in-place using the Cholesky decomposition of $A which must have been previously computed by cholesky-cdecomp.
On input $x should contain the right-hand side b, which is replaced by the solution on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 cholesky-cinvert(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> Int)

These functions compute the inverse of a matrix from its Cholesky decomposition $A, which must have been previously computed by cholesky-cdecomp.
On output, the inverse is stored in-place in $A.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 hermtd-cdecomp(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector::Complex64)

This function factorizes the hermitian matrix $A into the symmetric tridiagonal decomposition UTT(U). On output the diagonal and subdiagonal part of the input matrix $A contain the tridiagonal matrix T. The remaining lower triangular part of the input matrix contains the Householder vectors which, together with the Householder coefficients tau returned as a Math::Libgsl::Vector::Complex64 object, encode the orthogonal matrix U. The upper triangular part of $A and the imaginary parts of the diagonal are not referenced.
In case of error a failure object is returned.

=head3 hermtd-cunpack(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector::Complex64 $tau  where *.vector.size == $A.matrix.size1 - 1 --> List)

This function unpacks the encoded symmetric tridiagonal decomposition ($A, $tau) obtained from hermtd-cdecomp.
The function returns a List: the unitary matrix $U as a Math::Libgsl::Matrix::Complex64, the real vector of diagonal elements $diag as a Math::Libgsl::Vector, and the real vector of subdiagonal elements $subdiag as a Math::Libgsl::Vector.
In case of error a failure object is returned.

=head3 hermtd-cunpack-T(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> List)

This function unpacks the diagonal and subdiagonal of the encoded symmetric tridiagonal decomposition ($A, $tau) obtained from hermtd-cdecomp.
The function returns a List of two real Math::Libgsl::Vector: $diag and $subdiag.
In case of error a failure object is returned.

=head3 householder-ctransform(Math::Libgsl::Vector::Complex64 $w --> Complex)

This function prepares a Householder transformation H = I − τvT(v) which can be used to zero all the elements of the input vector $w except the first. On output the Householder vector v is stored in $w and the scalar τ is returned. The householder vector v is normalized so that v[0] = 1, however this 1 is not stored in the output vector. Instead, $w[0] is set to the first element of the transformed vector, so that if u = Hw, w[0] = u[0] on output and the remainder of u is zero.
This function returns a Complex: $τ.

=head3 householder-chm(Complex $tau, Math::Libgsl::Vector::Complex64 $v, Math::Libgsl::Matrix::Complex64 $A --> Int)

This function applies the Householder matrix H defined by the scalar $tau and the vector $v to the left-hand side of the matrix $A. On output the result HA is stored in $A.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 householder-cmh(Complex $tau, Math::Libgsl::Vector::Complex64 $v, Math::Libgsl::Matrix::Complex64 $A --> Int)

This function applies the Householder matrix H defined by the scalar $tau and the vector $v to the right-hand side of the matrix $A. On output the result AH is stored in $A.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 householder-chv(Complex $tau, Math::Libgsl::Vector::Complex64 $v, Math::Libgsl::Vector::Complex64 $w --> Int)

This function applies the Householder transformation H defined by the scalar $tau and the vector $v to the vector $w. On output the result Hw is stored in w.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tri-cinvert(Int $Uplo, Int $Diag, Math::Libgsl::Matrix::Complex64 $A --> Int)

This function is available only from the C library v2.6.
This function computes the in-place inverse of the triangular matrix $T, stored in the lower triangle when $Uplo = CblasLower and upper triangle when $Uplo = CblasUpper. The parameter $Diag = CblasUnit, CblasNonUnit specifies whether the matrix is unit triangular.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tri-cLHL(Math::Libgsl::Matrix::Complex64 $L --> Int)

This function is available only from the C library v2.6.
This function computes the product LT(L) in-place and stores it in the lower triangle of $L on output.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head3 tri-cUL(Math::Libgsl::Matrix::Complex64 $LU --> Int)

This function is available only from the C library v2.6.
This function compute the product $UL where U is upper triangular and L is unit lower triangular, stored in $LU, as computed by LU-cdecomp.
This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

=head1 C Library Documentation

For more details on libgsl see L<https://www.gnu.org/software/gsl/>.
The excellent C Library manual is available here L<https://www.gnu.org/software/gsl/doc/html/index.html>, or here L<https://www.gnu.org/software/gsl/doc/latex/gsl-ref.pdf> in PDF format.

=head1 Prerequisites

This module requires the libgsl library to be installed. Please follow the instructions below based on your platform:

=head2 Debian Linux and Ubuntu 20.04

=begin code
sudo apt install libgsl23 libgsl-dev libgslcblas0
=end code

That command will install libgslcblas0 as well, since it's used by the GSL.

=head2 Ubuntu 18.04

libgsl23 and libgslcblas0 have a missing symbol on Ubuntu 18.04.
I solved the issue installing the Debian Buster version of those three libraries:

=item L<http://http.us.debian.org/debian/pool/main/g/gsl/libgslcblas0_2.5+dfsg-6_amd64.deb>
=item L<http://http.us.debian.org/debian/pool/main/g/gsl/libgsl23_2.5+dfsg-6_amd64.deb>
=item L<http://http.us.debian.org/debian/pool/main/g/gsl/libgsl-dev_2.5+dfsg-6_amd64.deb>

=head1 Installation

To install it using zef (a module management tool):

=begin code
$ zef install Math::Libgsl::LinearAlgebra
=end code

=head1 AUTHOR

Fernando Santagata <nando.santagata@gmail.com>

=head1 COPYRIGHT AND LICENSE

Copyright 2020 Fernando Santagata

This library is free software; you can redistribute it and/or modify it under the Artistic License 2.0.

=end pod
