use v6;

unit class Math::Libgsl::LinearAlgebra::Complex64:ver<0.0.4>:auth<zef:FRITH>;

use Math::Libgsl::Raw::LinearAlgebra :ALL;
use NativeCall;
use Math::Libgsl::Vector::Complex64;
use Math::Libgsl::Matrix::Complex64;
use Math::Libgsl::Raw::Complex;
use Math::Libgsl::Permutation;
use Math::Libgsl::Exception;
use Math::Libgsl::Constants;

# LU Decomposition
sub LU-cdecomp(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> List) is export
{
  my Math::Libgsl::Permutation $p .= new: $A.matrix.size1;
  my int32 $s = 0;
  my $ret = gsl_linalg_complex_LU_decomp($A.matrix, $p.p, $s);
  fail X::Libgsl.new: errno => $ret, error => "Error in LU-decomp" if $ret ≠ GSL_SUCCESS;
  return $s, $p;
}

sub LU-csolve(Math::Libgsl::Matrix::Complex64 $LU where *.matrix.size1 == $LU.matrix.size2,
              Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1,
              Math::Libgsl::Vector::Complex64 $b where *.vector.size == $LU.matrix.size1
              --> Math::Libgsl::Vector::Complex64) is export
{
  my Math::Libgsl::Vector::Complex64 $x .= new: $LU.matrix.size2;
  my $ret = gsl_linalg_complex_LU_solve($LU.matrix, $p.p, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in LU-csolve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub LU-csvx(Math::Libgsl::Matrix::Complex64 $LU where *.matrix.size1 == $LU.matrix.size2,
            Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1,
            Math::Libgsl::Vector::Complex64 $x where *.vector.size == $LU.matrix.size1 --> Int) is export
{
  gsl_linalg_complex_LU_svx($LU.matrix, $p.p, $x.vector);
}

sub LU-crefine(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2,
               Math::Libgsl::Matrix::Complex64 $LU where { $LU.matrix.size1 == $LU.matrix.size2 && $A.matrix.size1 == $LU.matrix.size2 },
               Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1,
               Math::Libgsl::Vector::Complex64 $b where *.vector.size == $LU.matrix.size1,
               Math::Libgsl::Vector::Complex64 $x where *.vector.size == $LU.matrix.size1 --> Int) is export
{
  my Math::Libgsl::Vector::Complex64 $work .= new($LU.matrix.size2);
  gsl_linalg_complex_LU_refine($A.matrix, $LU.matrix, $p.p, $b.vector, $x.vector, $work.vector);
}

sub LU-cinvert(Math::Libgsl::Matrix::Complex64 $LU,
               Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1
               --> Math::Libgsl::Matrix::Complex64) is export
{
  my Math::Libgsl::Matrix::Complex64 $inverse .= new($LU.matrix.size1, $LU.matrix.size1);
  my $ret = gsl_linalg_complex_LU_invert($LU.matrix, $p.p, $inverse.matrix);
  fail X::Libgsl.new: errno => $ret, error => "Error in LU-invert" if $ret ≠ GSL_SUCCESS;
  $inverse;
}

sub LU-cdet(Math::Libgsl::Matrix::Complex64 $LU, Int $signum where * ~~ -1|1 --> Complex) is export
{
  my gsl_complex $res = alloc_gsl_complex;
  mgsl_linalg_complex_LU_det($LU.matrix, $signum, $res);
  my Complex $nc = $res.dat[0] + i * $res.dat[1];
  free_gsl_complex($res);
  $nc;
}

sub LU-clndet(Math::Libgsl::Matrix::Complex64 $LU --> Num) is export
{
  gsl_linalg_complex_LU_lndet($LU.matrix);
}

sub LU-csgndet(Math::Libgsl::Matrix::Complex64 $LU, Int $signum where * ~~ -1|1 --> Complex) is export
{
  my gsl_complex $res = alloc_gsl_complex;
  mgsl_linalg_complex_LU_sgndet($LU.matrix, $signum, $res);
  my Complex $nc = $res.dat[0] + i * $res.dat[1];
  free_gsl_complex($res);
  $nc;
}

# Cholesky Decomposition
sub cholesky-cdecomp(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> Int) is export
{
  gsl_linalg_complex_cholesky_decomp($A.matrix);
}

sub cholesky-csolve(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2,
                    Math::Libgsl::Vector::Complex64 $b where *.vector.size == $A.matrix.size1
                    --> Math::Libgsl::Vector::Complex64) is export
{
  my Math::Libgsl::Vector::Complex64 $x .= new: $A.matrix.size2;
  my $ret = gsl_linalg_complex_cholesky_solve($A.matrix, $b.vector, $x.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in cholesky-csolve" if $ret ≠ GSL_SUCCESS;
  return $x;
}

sub cholesky-csvx(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2,
                  Math::Libgsl::Vector::Complex64 $x where *.vector.size == $A.matrix.size1 --> Int) is export
{
  gsl_linalg_complex_cholesky_svx($A.matrix, $x.vector);
}

sub cholesky-cinvert(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> Int) is export
{
  gsl_linalg_complex_cholesky_invert($A.matrix);
}

# Tridiagonal Decomposition of Hermitian Matrices
sub hermtd-cdecomp(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2
                   --> Math::Libgsl::Vector::Complex64) is export
{
  my Math::Libgsl::Vector::Complex64 $tau .= new: $A.matrix.size1 - 1;
  my $ret = gsl_linalg_hermtd_decomp($A.matrix, $tau.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in hermtd-cdecomp" if $ret ≠ GSL_SUCCESS;
  return $tau;
}

sub hermtd-cunpack(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2,
                   Math::Libgsl::Vector::Complex64 $tau  where *.vector.size == $A.matrix.size1 - 1 --> List) is export
{
  my $M = $A.matrix.size1;
  my Math::Libgsl::Matrix::Complex64 $U .= new: $M, $M;
  my Math::Libgsl::Vector $diag         .= new: $A.matrix.size1;
  my Math::Libgsl::Vector $subdiag      .= new: $A.matrix.size1 - 1;
  my $ret = gsl_linalg_hermtd_unpack($A.matrix, $tau.vector, $U.matrix, $diag.vector, $subdiag.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in hermtd-cunpack" if $ret ≠ GSL_SUCCESS;
  return $U, $diag, $subdiag;
}

sub hermtd-cunpack-T(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> List) is export
{
  my $M = $A.matrix.size1;
  my Math::Libgsl::Vector $diag    .= new: $A.matrix.size1;
  my Math::Libgsl::Vector $subdiag .= new: $A.matrix.size1 - 1;
  my $ret = gsl_linalg_hermtd_unpack_T($A.matrix, $diag.vector, $subdiag.vector);
  fail X::Libgsl.new: errno => $ret, error => "Error in hermtd-cunpack-T" if $ret ≠ GSL_SUCCESS;
  return $diag, $subdiag;
}

# Householder Transformations
sub householder-ctransform(Math::Libgsl::Vector::Complex64 $w --> Complex) is export
{
  my gsl_complex $res = alloc_gsl_complex;
  mgsl_linalg_complex_householder_transform($w.vector, $res);
  my Complex $nc = $res.dat[0] + i * $res.dat[1];
  free_gsl_complex($res);
  $nc;
}

sub householder-chm(Complex $tau, Math::Libgsl::Vector::Complex64 $v, Math::Libgsl::Matrix::Complex64 $A --> Int) is export
{
  my gsl_complex $τ = alloc_gsl_complex;
  mgsl_complex_rect($tau.re, $tau.im, $τ);
  my $res = mgsl_linalg_complex_householder_hm($τ, $v.vector, $A.matrix);
  free_gsl_complex($res);
  $res;
}

sub householder-cmh(Complex $tau, Math::Libgsl::Vector::Complex64 $v, Math::Libgsl::Matrix::Complex64 $A --> Int) is export
{
  my gsl_complex $τ = alloc_gsl_complex;
  mgsl_complex_rect($tau.re, $tau.im, $τ);
  my $res = mgsl_linalg_complex_householder_mh($τ, $v.vector, $A.matrix);
  free_gsl_complex($res);
  $res;
}

sub householder-chv(Complex $tau, Math::Libgsl::Vector::Complex64 $v, Math::Libgsl::Vector::Complex64 $w --> Int) is export
{
  my gsl_complex $τ = alloc_gsl_complex;
  mgsl_complex_rect($tau.re, $tau.im, $τ);
  my $res = mgsl_linalg_complex_householder_hv($τ, $v.vector, $w.vector);
  free_gsl_complex($res);
  $res;
}

# Triangular Systems
sub tri-cinvert(Int $Uplo, Int $Diag, Math::Libgsl::Matrix::Complex64 $A --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-decomp: version < v2.6" if $gsl-version < v2.6;
  gsl_linalg_complex_tri_invert($Uplo, $Diag, $A.matrix);
}

sub tri-cLHL(Math::Libgsl::Matrix::Complex64 $L --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-decomp: version < v2.6" if $gsl-version < v2.6;
  gsl_linalg_complex_tri_LHL($L.matrix);
}

sub tri-cUL(Math::Libgsl::Matrix::Complex64 $LU --> Int) is export
{
  fail X::Libgsl.new: errno => GSL_FAILURE, error => "Error in ldlt-decomp: version < v2.6" if $gsl-version < v2.6;
  gsl_linalg_complex_tri_UL($LU.matrix);
}
