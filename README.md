[![Actions Status](https://github.com/frithnanth/raku-Math-Libgsl-LinearAlgebra/workflows/test/badge.svg)](https://github.com/frithnanth/raku-Math-Libgsl-LinearAlgebra/actions) [![Build Status](https://travis-ci.org/frithnanth/raku-Math-Libgsl-LinearAlgebra.svg?branch=master)](https://travis-ci.org/frithnanth/raku-Math-Libgsl-LinearAlgebra)

NAME
====

Math::Libgsl::LinearAlgebra - An interface to libgsl, the Gnu Scientific Library - Linear Algebra.

SYNOPSIS
========

```perl6
use Math::Libgsl::LinearAlgebra;
```

DESCRIPTION
===========

Math::Libgsl::LinearAlgebra is an interface to the linear algebra functions of libgsl, the GNU Scientific Library. This package provides both the low-level interface to the C library (Raw) and a more comfortable interface layer for the Raku programmer.

This module provides functions for Num and Complex data types.

Num
---

### LU-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> List)

This function factorizes the matrix A into the LU decomposition PA = LU. The factorization is done in place, so on output the diagonal and upper triangular (or trapezoidal) part of the input matrix A contain the matrix U. The lower triangular (or trapezoidal) part of the input matrix (excluding the diagonal) contains L. The diagonal elements of L are unity, and are not stored. The return value is a List: the sign of the permutation and a permutation object, which encodes the permutation matrix P. In case of error a failure object is returned.

### LU-solve(Math::Libgsl::Matrix $LU where *.matrix.size1 == $LU.matrix.size2, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $LU.matrix.size1 --> Math::Libgsl::Vector)

This function solves the square system Ax = b using the LU decomposition of A into (LU, p) given by the output of LU-decomp. In case of error a failure object is returned.

### LU-svx(Math::Libgsl::Matrix $LU where *.matrix.size1 == $LU.matrix.size2, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector $x where *.size == $LU.matrix.size1 --> Int)

This function solves the square system Ax = b in-place using the precomputed LU decomposition of A into (LU, p). On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LU-refine(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Matrix $LU where { $LU.matrix.size1 == $LU.matrix.size2 && $A.matrix.size1 == $LU.matrix.size2 }, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $LU.matrix.size1, Math::Libgsl::Vector $x where *.vector.size == $LU.matrix.size1 --> Int)

This function applies an iterative improvement to x, the solution of Ax = b, from the precomputed LU decomposition of A into (LU, p). This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LU-invert(Math::Libgsl::Matrix $LU, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1 --> Math::Libgsl::Matrix)

This function computes the inverse of a matrix A from its LU decomposition (LU, p), returning the matrix inverse. In case of error a failure object is returned.

### LU-det(Math::Libgsl::Matrix $LU, Int $signum where * ~~ -1|1 --> Num)

This function computes the determinant of a matrix A from its LU decomposition, $LU, and the sign of the permutation, $signum. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LU-lndet(Math::Libgsl::Matrix $LU --> Num)

This function computes the determinant the logarithm of the absolute value of the determinant of a matrix A, ln |det(A)|, from its LU decomposition, $LU. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LU-sgndet(Math::Libgsl::Matrix $LU, Int $signum where * ~~ -1|1 --> Int)

This function computes the sign or phase factor of the determinant of a matrix A, det(A)/|det(A)| from its LU decomposition, $LU. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### QR-decomp(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector)

This function factorizes the M-by-N matrix $A into the QR decomposition A = QR. On output the diagonal and upper triangular part of the input matrix contain the matrix R. The returned vector and the columns of the lower triangular part of the matrix A contain the Householder coefficients and Householder vectors which encode the orthogonal matrix Q. In case of error a failure object is returned.

### QR-solve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector)

This function solves the square system Ax = b using the QR decomposition of A held in ($QR, $tau) which must have been computed previously with QR-decomp(). In case of error a failure object is returned.

### QR-svx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $QR.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size1 --> Int)

This function solves the square system Ax = b in-place using the QR decomposition of A held in ($QR, $tau) which must have been computed previously by QR-decomp(). On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### QR-lssolve(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> List)

This function finds the least squares solution to the overdetermined system Ax = b where the matrix A has more rows than columns. The least squares solution minimizes the Euclidean norm of the residual, ||Ax − b||. The routine requires as input the QR decomposition of A into ($QR, $tau) given by QR-decomp(). The function returns a List of two Math::Libgsl::Vector objects: the solution x and the residual. In case of error a failure object is returned.

### QR-QTvec(Math::Libgsl::Matrix $QR, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Vector $v where *.vector.size == $QR.matrix.size1 --> Int)

These function applies the matrix T(Q) encoded in the decomposition (QR, tau) to the vector $v, storing the result T(Q) v in $v. The matrix multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix T(Q). This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### QR-Qvec(Math::Libgsl::Matrix $QR, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Vector $v where *.vector.size == $QR.matrix.size1 --> Int)

This function applies the matrix Q encoded in the decomposition ($QR, $tau) to the vector $v, storing the result Qv in $v. The matrix multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix Q. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### QR-QTmat(Math::Libgsl::Matrix $QR, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Matrix $B where *.matrix.size1 == $QR.matrix.size1 --> Int)

This function applies the matrix T(Q) encoded in the decomposition ($QR, $tau) to the M-by-K matrix $B, storing the result T(Q) B in $B. The matrix multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix T(Q). This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### QR-Rsolve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector)

This function solves the triangular system Rx = b and returns the Math::Libgsl::Vector object $x. In case of error a failure object is returned.

### QR-Rsvx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size2 --> Int)

This function solves the triangular system Rx = b for x in-place. On input $x should contain the right-hand side b and is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### QR-unpack(Math::Libgsl::Matrix $QR, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2) --> List)

This function unpacks the encoded QR decomposition ($QR, $tau) into the matrices Q and R. The function returns a List of two Math::Libgsl::Matrix objects: $Q which is M-by-M and $R which is M-by-N. In case of error a failure object is returned.

### QR-QRsolve(Math::Libgsl::Matrix $Q, Math::Libgsl::Matrix $R where { $R.matrix.size1 == $R.matrix.size2 && $Q.matrix.size1 == $R.matrix.size1 }, Math::Libgsl::Vector $b where *.vector.size == $R.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system Rx = T(Q) b for x. It can be used when the QR decomposition of a matrix is available in unpacked form as ($Q, $R). The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### QR-update(Math::Libgsl::Matrix $Q, Math::Libgsl::Matrix $R where { $Q.matrix.size1 == $R.matrix.size1 && $Q.matrix.size2 == $R.matrix.size1 }, Math::Libgsl::Vector $w where *.vector.size == $R.matrix.size1, Math::Libgsl::Vector $v where *.vector.size == $R.matrix.size2 --> Int)

This function performs a rank-1 update wT(v) of the QR decomposition ($Q, $R). The update is given by Q'R' = Q(R+wT(v)) where the output matrices $Q and $R are also orthogonal and right triangular. Note that $w is destroyed by the update. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### R-solve(Math::Libgsl::Matrix $R where *.matrix.size1 == $R.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $R.matrix.size1 --> Math::Libgsl::Vector)

This function solves the triangular system Rx = b for the N-by-N matrix $R. The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### R-svx(Math::Libgsl::Matrix $R where *.matrix.size1 == $R.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $R.matrix.size2 --> Int)

This function solves the triangular system Rx = b in-place. On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### QRPT-decomp(Math::Libgsl::Matrix $A --> List)

This function factorizes the M-by-N matrix $A into the QRT(P) decomposition A = QRT(P). On output the diagonal and upper triangular part of the input matrix contain the matrix R. The function's output is a List of three objects: the Math::Libgsl::Vector $tau, the Math::Libgsl::Permutation $p and the sign of the permutation Int $signum. In case of error a failure object is returned.

### QRPT-decomp2(Math::Libgsl::Matrix $A --> List)

This function factorizes the matrix $A into the decomposition A = QRT(P) without modifying $A itself. The function returns a List: the Math::Libgsl::Matrix $Q, the Math::Libgsl::Matrix $R, the Math::Libgsl::Permutation $p, and the sign of the permutation Int $signum. In case of error a failure object is returned.

### QRPT-solve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector)

This function solves the square system Ax = b using the QRT(P) decomposition of A held in ($QR, $tau, $p) which must have been computed previously by QRPT-decomp. The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### QRPT-svx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size1, Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size2 --> Int)

This function solves the square system Ax = b in-place using the QRT(P) decomposition of A held in ($QR, $tau, $p). On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### QRPT-lssolve(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> List)

This function finds the least squares solution to the overdetermined system Ax = b where the matrix A has more rows than columns and is assumed to have full rank. The least squares solution minimizes the Euclidean norm of the residual, ||b − Ax||. The routine requires as input the QR decomposition of A into ($QR, $tau, $p) given by QRPT-decomp. The function returns a List of two Math::Libgsl::Vector objects: the solution x and the residual. In case of error a failure object is returned.

### QRPT-lssolve2(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($QR.matrix.size1, $QR.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1, Int $rank where 0 < * ≤ $QR.matrix.size2 --> List)

This function finds the least squares solution to the overdetermined system Ax = b where the matrix A has more rows than columns and has rank given by the input rank. If the user does not know the rank of A, it may be estimated by calling QRPT-rank. The routine requires as input the QR decomposition of A into ($QR, $tau, $p) given by QRPT-decomp. The function returns a List of two Math::Libgsl::Vector objects: the solution x and the residual. In case of error a failure object is returned.

### QRPT-QRsolve(Math::Libgsl::Matrix $Q where *.matrix.size1 == $Q.matrix.size2, Math::Libgsl::Matrix $R where {$R.matrix.size1 == $R.matrix.size2 && $R.matrix.size1 == $Q.matrix.size1}, Math::Libgsl::Permutation $p where *.p.size == $Q.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $Q.matrix.size1 --> Math::Libgsl::Vector)

This function solves the square system RT(P) x = T(Q) b for x. It can be used when the QR decomposition of a matrix is available in unpacked form as ($Q, $R). The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### QRPT-update(Math::Libgsl::Matrix $Q, Math::Libgsl::Matrix $R where { $Q.matrix.size1 == $R.matrix.size1 && $Q.matrix.size2 == $R.matrix.size1 }, Math::Libgsl::Permutation $p where *.p.size == $R.matrix.size1, Math::Libgsl::Vector $w where *.vector.size == $R.matrix.size1, Math::Libgsl::Vector $v where *.vector.size == $R.matrix.size2 --> Int)

This function performs a rank-1 update wT(v) of the QRT(P) decomposition ($Q, $R, $p). The update is given by Q' R' = Q(R + wT(v) P) where the output matrices Q' and R' are also orthogonal and right triangular. Note that $w is destroyed by the update. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### QRPT-Rsolve(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $QR.matrix.size1 --> Math::Libgsl::Vector)

This function solves the triangular system RT(P) x = b for the N-by-N matrix R contained in $QR. The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### QRPT-Rsvx(Math::Libgsl::Matrix $QR where *.matrix.size1 == $QR.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $QR.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $QR.matrix.size2 --> Int)

This function solves the triangular system RT(P) x = b in-place for the N-by-N matrix R contained in $QR. On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### QRPT-rank(Math::Libgsl::Matrix $QR, Num() $tolerance --> Int)

This function returns the rank of the triangular matrix R contained in $QR.

### QRPT-rcond(Math::Libgsl::Matrix $QR where *.matrix.size1 ≥ $QR.matrix.size2 --> Num)

This function returns the reciprocal condition number (using the 1-norm) of the R factor, stored in the upper triangle of $QR. In case of error a failure object is returned.

### LQ-decomp(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector)

This function factorizes the M-by-N matrix $A into the LQ decomposition A = LQ. On output the diagonal and lower trapezoidal part of the input matrix contain the matrix L. This function returns the Math::Libgsl::Vector $tau. The vector $tau and the elements above the diagonal of the matrix $A contain the Householder coefficients and Householder vectors which encode the orthogonal matrix Q. In case of error a failure object is returned.

### LQ-solve-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2), Math::Libgsl::Vector $b where *.vector.size == $LQ.matrix.size2 --> Math::Libgsl::Vector)

This function finds the solution to the system Ax = b. The routine requires as input the LQ decomposition of A into ($LQ, $tau) given by LQ-decomp. The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### LQ-svx-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2), --> Math::Libgsl::Vector)

This function finds the solution to the system Ax = b. The routine requires as input the LQ decomposition of A into ($LQ, $tau) given by LQ-decomp. On input $x should contain the right-hand side b, which is replaced by the solution on output. In case of error a failure object is returned.

### LQ-lssolve-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 ≥ $LQ.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2), Math::Libgsl::Vector $b where *.vector.size == $LQ.matrix.size1 --> List)

This function finds the minimum norm least squares solution to the underdetermined system Ax = b, where the M-by-N matrix A has M ≤ N. The routine requires as input the LQ decomposition of A into ($LQ, $tau) given by LQ-decomp. The function returns a List of two Math::Libgsl::Vector objects: the solution x and the residual. In case of error a failure object is returned.

### LQ-Lsolve-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $LQ.matrix.size1 --> Math::Libgsl::Vector)

The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### LQ-Lsvx-T(Math::Libgsl::Matrix $LQ where *.matrix.size1 == $LQ.matrix.size2 --> Math::Libgsl::Vector)

On input $x should contain the right-hand side b, which is replaced by the solution on output. In case of error a failure object is returned.

### L-solve-T(Math::Libgsl::Matrix $L where *.matrix.size1 == $L.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $L.matrix.size2 --> Math::Libgsl::Vector)

The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### LQ-vecQ(Math::Libgsl::Matrix $LQ, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2), Math::Libgsl::Vector $v where *.vector.size == $LQ.matrix.size1 --> Int)

This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LQ-vecQT(Math::Libgsl::Matrix $LQ, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2), Math::Libgsl::Vector $v where *.vector.size == $LQ.matrix.size1 --> Int)

This function applies T(Q) to the vector v, storing the result T(Q) v in $v on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LQ-unpack(Math::Libgsl::Matrix $LQ, Math::Libgsl::Vector $tau where *.vector.size == min($LQ.matrix.size1, $LQ.matrix.size2) --> List)

This function unpacks the encoded LQ decomposition ($LQ, $tau). The function outputs a List: the Math::Libgsl::Matrix $Q and the Math::Libgsl::Matrix $L. In case of error a failure object is returned.

### LQ-update(Math::Libgsl::Matrix $Q, Math::Libgsl::Matrix $L where { $L.matrix.size2 == $Q.matrix.size1 && $L.matrix.size2 == $Q.matrix.size2 }, Math::Libgsl::Vector $v where *.vector.size == $L.matrix.size1, Math::Libgsl::Vector $w where *.vector.size == $L.matrix.size2 --> Int)

This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LQ-LQsolve(Math::Libgsl::Matrix $Q, Math::Libgsl::Matrix $L where { $L.matrix.size1 == $L.matrix.size2 && $Q.matrix.size1 == $L.matrix.size2 }, Math::Libgsl::Vector $b where *.vector.size == $L.matrix.size2 --> Math::Libgsl::Vector)

In case of error a failure object is returned.

### COD-decomp(Math::Libgsl::Matrix $A --> List)

### COD-decomp-e(Math::Libgsl::Matrix $A, Num() $tolerance --> List)

These functions factor the M-by-N matrix $A into the decomposition A = QRZT(P). On output the matrix R₁₁ is stored in the upper rank-by-rank block of $A. The matrices Q and Z are encoded in packed storage in $A on output. This function outputs a List: two Math::Libgsl::Vector objects which contain the Householder scalars corresponding to the matrices Q and Z respectively, a Math::Libgsl::Permutation object which contain the permutation matrix P, and the rank of $A. In case of error a failure object is returned.

### COD-lssolve(Math::Libgsl::Matrix $QRZT where *.matrix.size1 ≥ $QRZT.matrix.size2, Math::Libgsl::Vector $tau-Q where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QRZT.matrix.size2, Int $rank where * ≤ min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Vector $b where *.vector.size == $QRZT.matrix.size1 --> List)

This function finds the unique minimum norm least squares solution to the overdetermined system Ax = b where the matrix $A has more rows than columns. The least squares solution minimizes the Euclidean norm of the residual, ||b − Ax|| as well as the norm of the solution ||x||. The routine requires as input the QRZT decomposition of $A into ($QRZT, $tau_Q, $tau_Z, $p, $rank) given by COD-decomp. The function outputs a List: a Math::Libgsl::Vector object which is the solution x, and a Math::Libgsl::Vector object which stores the residual b − Ax. In case of error a failure object is returned.

### COD-lssolve2(Math::Libgsl::Matrix $QRZT where *.matrix.size1 ≥ $QRZT.matrix.size2, Math::Libgsl::Vector $tau-Q where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Permutation $p where *.p.size == $QRZT.matrix.size2, Int $rank where * ≤ min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Vector $b where *.vector.size == $QRZT.matrix.size1, Num() $lambda --> List)

This function finds the solution to the regularized least squares problem in Tikhonov standard form, minₓ||b − Ax||² + λ²||x||². The routine requires as input the QRZT decomposition of A into ($QRZT, $tau_Q, $tau_Z, $p, $rank) given by COD-decomp. The parameter λ is supplied in $lambda. The function outputs a List: a Math::Libgsl::Vector object which is the solution x, and a Math::Libgsl::Vector object which stores the residual b − Ax. In case of error a failure object is returned.

### COD-unpack(Math::Libgsl::Matrix $QRZT, Math::Libgsl::Vector $tau-Q where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Int $rank where * ≤ min($QRZT.matrix.size1, $QRZT.matrix.size2) --> List)

This function unpacks the encoded QRZT decomposition ($QRZT, $tau_Q, $tau_Z, $rank). The function returns a List of three Math::Libgsl::Matrix objects: $Q, $R, $Z. In case of error a failure object is returned.

### COD-matZ(Math::Libgsl::Matrix $QRZT, Math::Libgsl::Vector $tau-Z where *.vector.size == min($QRZT.matrix.size1, $QRZT.matrix.size2), Math::Libgsl::Matrix $A where *.matrix.size2 == $QRZT.matrix.size2, Int $rank --> Int)

This function multiplies the input matrix $A on the right by Z, A’ = AZ using the encoded QRZT decomposition ($QRZT, $tau_Z, $rank). This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### SV-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List)

This function factorizes the M-by-N matrix $A into the singular value decomposition A = UST(V) for M ≥ N. On output the matrix $A is replaced by U. The function returns a List: a Math::Libgsl::Matrix object which contains the elements of V in untransposed form, and a Math::Libgsl::Vector object which contains the diagonal elements of the singular value matrix S. The singular values are non-negative and form a non-increasing sequence from S₁ to Sₙ. In case of error a failure object is returned.

### SV-decomp-mod(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List)

This function computes the SVD using the modified Golub-Reinsch algorithm, which is faster for M ≫ N. The function returns a List: a Math::Libgsl::Matrix object which contains the elements of V in untransposed form, and a Math::Libgsl::Vector object which contains the diagonal elements of the singular value matrix S. The singular values are non-negative and form a non-increasing sequence from S₁ to Sₙ. In case of error a failure object is returned.

### SV-decomp-jacobi(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List)

This function computes the SVD of the M-by-N matrix A using one-sided Jacobi orthogonalization for M ≥ N. The function returns a List: a Math::Libgsl::Matrix object which contains the elements of V in untransposed form, and a Math::Libgsl::Vector object which contains the diagonal elements of the singular value matrix S. The singular values are non-negative and form a non-increasing sequence from S₁ to Sₙ. In case of error a failure object is returned.

### SV-solve(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2, Math::Libgsl::Matrix $V where { $V.matrix.size1 == $V.matrix.size2 && $V.matrix.size1 == $A.matrix.size2 }, Math::Libgsl::Vector $S where *.vector.size == $A.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system Ax = b using the singular value decomposition (U, S, V) of A which must have been computed previously with COD-decomp. Only non-zero singular values are used in computing the solution. The parts of the solution corresponding to singular values of zero are ignored. Other singular values can be edited out by setting them to zero before calling this function. The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### SV-leverage(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector)

This function computes the statistical leverage values hᵢ of a matrix A using its singular value decomposition (U, S, V) previously computed with COD-decomp. The function returns a Math::Libgsl::Vector object which stores the diagonal values of the matrix A(T(A)A)⁻¹T(A) and depend only on the matrix U which is the input to this function. In case of error a failure object is returned.

### cholesky-decomp1(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int)

This function factorizes the symmetric, positive-definite square matrix $A into the Cholesky decomposition A = LT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used (the upper triangular part is ignored). On output the diagonal and lower triangular part of the input matrix $A contain the matrix L, while the upper triangular part contains the original matrix. If the matrix is not positive-definite then the decomposition will fail, returning the error code GSL_EDOM. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-solve(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector)

These functions solve the system Ax = b using the Cholesky decomposition of $A which must have been previously computed by cholesky-decomp1. The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### cholesky-svx(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $A.matrix.size1 --> Int)

This function solves the system Ax = b in-place using the Cholesky decomposition of $A which must have been previously computed by cholesky-decomp1. On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-solve-mat(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Matrix $B where *.matrix.size1 == $A.matrix.size1 --> Math::Libgsl::Matrix)

In case of error a failure object is returned.

### cholesky-svx-mat(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Matrix $X where *.matrix.size1 == $A.matrix.size1 --> Int)

This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-invert(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int)

These functions compute the inverse of a matrix from its Cholesky decomposition $A, which must have been previously computed by cholesky-decomp1. On output, the inverse is stored in-place in $A. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-decomp2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector)

This function calculates a diagonal scaling transformation S for the symmetric, positive-definite square matrix $A, and then computes the Cholesky decomposition SAS = LT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used (the upper triangular part is ignored). On output the diagonal and lower triangular part of the input matrix $A contain the matrix L, while the upper triangular part of the input matrix is overwritten with T(L) (the diagonal terms being identical for both L and L T). If the matrix is not positive-definite then the decomposition will fail, returning the error code GSL_EDOM. The function returns a Math::Libgsl::Vector object which stores the diagonal scale factors. In case of error a failure object is returned.

### cholesky-solve2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $s where *.vector.size == $A.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system (SAS)(S⁻¹ x) = Sb using the Cholesky decomposition of SAS held in the matrix $A which must have been previously computed by cholesky-decomp2. The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### cholesky-svx2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $s where *.size == $A.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $A.matrix.size2 --> Int)

This function solves the system (SAS)(S⁻¹ x) = Sb using the Cholesky decomposition of SAS held in the matrix $A which must have been previously computed by cholesky-decomp2. On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-decomp-unit(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector)

The function returns the Math::Libgsl::Vector $d. In case of error a failure object is returned.

### cholesky-scale(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector)

This function calculates a diagonal scaling transformation of the symmetric, positive definite matrix $A, such that SAS has a condition number within a factor of N of the matrix of smallest possible condition number over all possible diagonal scalings. The function outputs a Math::Libgsl::Vector object which contains the scale factors. In case of error a failure object is returned.

### cholesky-scale-apply(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $s where *.vector.size == $A.matrix.size2 --> Int)

This function applies the scaling transformation S to the matrix $A. On output, $A is replaced by SAS. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-rcond(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num)

This function estimates the reciprocal condition number (using the 1-norm) of the symmetric positive definite matrix A, using its Cholesky decomposition provided in $A. The function returns a Math::Libgsl::Vector object which stores the reciprocal condition number estimate, defined as 1/(||A||₁ · ||A⁻¹||₁). In case of error a failure object is returned.

### pcholesky-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Permutation)

This function factors the symmetric, positive-definite square matrix $A into the Pivoted Cholesky decomposition PAT(P) = LDT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used to construct the factorization. On output the diagonal of the input matrix $A stores the diagonal elements of D, and the lower triangular portion of $A contains the matrix L. Since L has ones on its diagonal these do not need to be explicitely stored. The upper triangular portion of $A is unmodified. The function returns the Math::Libgsl::Permutation object which stores the permutation matrix P. In case of error a failure object is returned.

### pcholesky-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system Ax = b using the Pivoted Cholesky decomposition of A held in the matrix $LDLT and permutation $p which must have been previously computed by pcholesky-decomp. The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### pcholesky-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int)

This function solves the system Ax = b using the Pivoted Cholesky decomposition of A held in the matrix $LDLT and permutation $p which must have been previously computed by pcholesky-decomp. On input, $x contains the right hand side vector b which is replaced by the solution vector on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### pcholesky-decomp2(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> List)

This function computes the pivoted Cholesky factorization of the matrix SAS, where the input matrix $A is symmetric and positive definite, and the diagonal scaling matrix S is computed to reduce the condition number of $A as much as possible. On input, the values from the diagonal and lower-triangular part of the matrix $A are used to construct the factorization. On output the diagonal of the input matrix $A stores the diagonal elements of D, and the lower triangular portion of $A contains the matrix L. Since L has ones on its diagonal these do not need to be explicitely stored. The upper triangular portion of $A is unmodified. The function returns a List: the permutation matrix P is stored in a Math::Libgsl::Permutation object, the diagonal scaling transformation is stored in a Math::Libgsl::Vector object. In case of error a failure object is returned.

### pcholesky-solve2(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $s where *.vector.size == $LDLT.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system (SAS)(S⁻¹x) = Sb using the Pivoted Cholesky decomposition of SAS held in the matrix $LDLT, permutation $p, and vector $S, which must have been previously computed by pcholesky-decomp2. The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### pcholesky-svx2(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $s where *.vector.size == $LDLT.matrix.size1, Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int)

This function solves the system (SAS)(S⁻¹x) = Sb using the Pivoted Cholesky decomposition of SAS held in the matrix $LDLT, permutation $p, and vector $S, which must have been previously computed by pcholesky-decomp2. On input, $x contains the right hand side vector b which is replaced by the solution vector on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### pcholesky-invert(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1 --> Math::Libgsl::Matrix)

This function computes the inverse of the matrix A, using the Pivoted Cholesky decomposition stored in $LDLT and $p. The function returns the Math::Libgsl::Matrix A⁻¹. In case of error a failure object is returned.

### pcholesky-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1 --> Num)

This function estimates the reciprocal condition number (using the 1-norm) of the symmetric positive definite matrix A, using its pivoted Cholesky decomposition provided in $LDLT. The function returns the reciprocal condition number estimate, defined as 1/(||A||₁ · ||A⁻¹||₁). This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### mcholesky-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Bool :$perturbation = True --> List)

This function factors the symmetric, indefinite square matrix $A into the Modified Cholesky decomposition P(A + E)T(P) = LDT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used to construct the factorization. On output the diagonal of the input matrix $A stores the diagonal elements of D, and the lower triangular portion of $A contains the matrix L. Since L has ones on its diagonal these do not need to be explicitely stored. The upper triangular portion of $A is unmodified. The function returns a List: the permutation matrix P, stored in a Math::Libgsl::Permutation object and the diagonal perturbation matrix, stored in a Math::Libgsl::Vector object. In case of error a failure object is returned.

### mcholesky-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1 --> Math::Libgsl::Vector)

This function solves the perturbed system (A + E)x = b using the Cholesky decomposition of A + E held in the matrix $LDLT and permutation $p which must have been previously computed by mcholesky-decomp. The function returns the Math::Libgsl::Vector $x. In case of error a failure object is returned.

### mcholesky-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1, Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int)

This function solves the perturbed system (A + E)x = b using the Cholesky decomposition of A + E held in the matrix $LDLT and permutation $p which must have been previously computed by mcholesky-decomp. On input, $x contains the right hand side vector b which is replaced by the solution vector on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### mcholesky-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1 --> Num)

This function estimates the reciprocal condition number (using the 1-norm) of the perturbed matrix A+E, using its pivoted Cholesky decomposition provided in $LDLT. The function returns the reciprocal condition number estimate, defined as 1/(||A + E||₁ · ||(A + E)⁻¹||₁). In case of error a failure object is returned.

### mcholesky-invert(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Permutation $p where *.p.size == $LDLT.matrix.size1 --> Math::Libgsl::Matrix)

The function returns the inverted matrix. In case of error a failure object is returned.

### ldlt-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int)

This function factorizes the symmetric, non-singular square matrix $A into the decomposition A = LDT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used. The upper triangle of $A is used as temporary workspace. On output the diagonal of $A contains the matrix D and the lower triangle of $A contains the unit lower triangular matrix L. The matrix 1-norm, ||A||₁ is stored in the upper right corner on output This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error. This function is available only in the C library starting from v2.6.

### ldlt-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system Ax = b using the LDT(L) decomposition of A held in the matrix $LDLT which must have been previously computed by ldlt-decomp. In case of error a failure object is returned. This function is available only in the C library starting from v2.6.

### ldlt-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int)

This function solves the system Ax = b using the LDT(L) decomposition of A held in the matrix $LDLT which must have been previously computed by ldlt-decomp. On input $x should contain the right-hand side b, which is replaced by the solution on output. In case of error a failure object is returned. This function is available only in the C library starting from v2.6.

### ldlt-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2 --> Num)

This function estimates the reciprocal condition number (using the 1-norm) of the symmetric nonsingular matrix A, using its LDT(L) decomposition provided in $LDLT. The function returns the reciprocal condition number estimate, defined as 1/(||A + E||₁ · ||(A + E)⁻¹||₁). In case of error a failure object is returned. This function is available only in the C library starting from v2.6.

### symmtd-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector)

This function factorizes the symmetric square matrix $A into the symmetric tridiagonal decomposition QTT(Q). On output the diagonal and subdiagonal part of the input matrix $A contain the tridiagonal matrix T. The remaining lower triangular part of the input matrix contains the Householder vectors which, together with the Householder coefficients tau returned as a Math::Libgsl::Vector object, encode the orthogonal matrix Q. The upper triangular part of $A is not referenced. In case of error a failure object is returned.

### symmtd-unpack(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $A.matrix.size1 - 1 --> List)

This function unpacks the encoded symmetric tridiagonal decomposition ($A, $tau) obtained from symmtd-decomp. The function returns a List: the orthogonal matrix $Q as a Math::Libgsl::Matrix, the vector of diagonal elements $diag as a Math::Libgsl::Vector, and the vector of subdiagonal elements $subdiag as a Math::Libgsl::Vector. In case of error a failure object is returned.

### symmtd-unpack-T(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> List)

This function unpacks the diagonal and subdiagonal of the encoded symmetric tridiagonal decomposition ($A, $tau) obtained from symmtd-decomp. The function returns a List of two Math::Libgsl::Vector: $diag and $subdiag. In case of error a failure object is returned.

### hessenberg-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector)

This function computes the Hessenberg decomposition of the matrix $A by applying the similarity transformation H = T(U)AU. On output, H is stored in the upper portion of $A. The information required to construct the matrix U is stored in the lower triangular portion of $A. U is a product of N−2 Householder matrices. The Householder vectors are stored in the lower portion of $A (below the subdiagonal) and the Householder coefficients are returned as a Math::Libgsl::Vector. In case of error a failure object is returned.

### hessenberg-unpack(Math::Libgsl::Matrix $H where *.matrix.size1 == $H.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $H.matrix.size1 --> Math::Libgsl::Matrix)

This function constructs the orthogonal matrix U from the information stored in the Hessenberg matrix $H along with the vector $tau. $H and $tau are outputs from hessenberg-decomp. The function returns the Math::Libgsl::Matrix which stores $U. In case of error a failure object is returned.

### hessenberg-unpack-accum(Math::Libgsl::Matrix $H where *.matrix.size1 == $H.matrix.size2, Math::Libgsl::Vector $tau where *.vector.size == $H.matrix.size1 --> Math::Libgsl::Matrix)

This function is similar to hessenberg-unpack, except it accumulates the matrix U into V, so that V′ = VU. The function returns the Math::Libgsl::Matrix which stores $V. In case of error a failure object is returned.

### hessenberg-set-zero(Math::Libgsl::Matrix $H where *.matrix.size1 == $H.matrix.size2 --> Int)

This function sets the lower triangular portion of $H, below the subdiagonal, to zero. It is useful for clearing out the Householder vectors after calling hessenberg-decomp. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### hesstri-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Matrix $B where { $B.matrix.size1 == $A.matrix.size1 && $B.matrix.size2 == $A.matrix.size2}, Bool :$similarity = True --> List)

This function computes the Hessenberg-Triangular decomposition of the matrix pair ($A, $B). On output, H is stored in $A, and R is stored in $B. If the $similarity Bool parameter is True (default), then the similarity transformations are returned as a List of Math::Libgsl::Matrix objects. In case of error a failure object is returned.

### bidiag-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List)

This function factorizes the M-by-N matrix $A into bidiagonal form UBT(V). The diagonal and superdiagonal of the matrix B are stored in the diagonal and superdiagonal of $A. The orthogonal matrices U and V are stored as compressed Householder vectors in the remaining elements of $A. This function returns a List of Math::Libgsl::Vector objects: the two Householder coefficients $tau_U and $tau_V. In case of error a failure object is returned.

### bidiag-unpack(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2, Math::Libgsl::Vector $tau_U where *.vector.size == $A.matrix.size2, Math::Libgsl::Vector $tau_V where *.vector.size == $A.matrix.size2 - 1 --> List)

This function unpacks the bidiagonal decomposition of $A produced by bidiag-decomp, ($A, $tau_U, $tau_V) into the separate orthogonal matrices U, V and the diagonal vector diag and superdiagonal superdiag. Note that U is stored as a compact M-by-N orthogonal matrix satisfying T(U)U = I for efficiency. The function returns a List of four objects: the Math::Libgsl::Matrix $U and $V, and the Math::Libgsl::Vector $diag and $sdiag. In case of error a failure object is returned.

### bidiag-unpack2(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2, Math::Libgsl::Vector $tau_U where *.vector.size == $A.matrix.size2, Math::Libgsl::Vector $tau_V where *.vector.size == $A.matrix.size2 - 1 --> Math::Libgsl::Matrix)

This function unpacks the bidiagonal decomposition of $A produced by bidiag-decomp, ($A, $tau_U, $tau_V) into the separate orthogonal matrices U, V and the diagonal vector diag and superdiagonal superdiag. The matrix U is stored in-place in $A. The function returns the Math::Libgsl::Matrix object $V. In case of error a failure object is returned.

### bidiag-unpack-B(Math::Libgsl::Matrix $A where *.matrix.size1 ≥ $A.matrix.size2 --> List)

This function unpacks the diagonal and superdiagonal of the bidiagonal decomposition of $A from bidiag-decomp. The function returns a List of two objects: the Math::Libgsl::Vector $diag and $sdiag. In case of error a failure object is returned.

### givens(Num() $a, Num() $b --> List)

This function computes c = cos θ and s = sin θ so that the Givens matrix G(θ) acting on the vector (a, b) produces (r, 0), with r = √ a² + b². The function returns a List of Num: the c and s elements of the Givens matrix.

### givens-gv(Math::Libgsl::Vector $v, Int, $i, Int $j, Num() $c, Num() $s)

This function applies the Givens rotation defined by c = cos θ and s = sin θ to the i and j elements of v. On output, (v(i), v(j)) ← G(θ)(v(i), v(j)). This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### householder-transform(Math::Libgsl::Vector $w --> Num)

This function prepares a Householder transformation H = I − τvT(v) which can be used to zero all the elements of the input vector $w except the first. On output the Householder vector v is stored in $w and the scalar τ is returned. The householder vector v is normalized so that v[0] = 1, however this 1 is not stored in the output vector. Instead, $w[0] is set to the first element of the transformed vector, so that if u = Hw, w[0] = u[0] on output and the remainder of u is zero. This function returns a Num: $τ.

### householder-hm(Num() $tau, Math::Libgsl::Vector $v, Math::Libgsl::Matrix $A --> Int)

This function applies the Householder matrix H defined by the scalar $tau and the vector $v to the left-hand side of the matrix $A. On output the result HA is stored in $A. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### householder-mh(Num() $tau, Math::Libgsl::Vector $v, Math::Libgsl::Matrix $A --> Int)

This function applies the Householder matrix H defined by the scalar $tau and the vector $v to the right-hand side of the matrix $A. On output the result AH is stored in $A. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### householder-hv(Num() $tau, Math::Libgsl::Vector $v, Math::Libgsl::Vector $w --> Int)

This function applies the Householder transformation H defined by the scalar $tau and the vector $v to the vector $w. On output the result Hw is stored in w. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### HH-solve(Math::Libgsl::Matrix $A where *.matrix.size1 ≤ $A.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector)

This function solves the system Ax = b directly using Householder transformations. $b is not modified. The matrix $A is destroyed by the Householder transformations. The function returns a Math::Libgsl::Vector object: the solution $x. In case of error a failure object is returned.

### HH-svx(Math::Libgsl::Matrix $A where *.matrix.size1 ≤ $A.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $A.matrix.size2 --> Int)

This function solves the system Ax = b in-place using Householder transformations. On input $x should contain the right-hand side b, which is replaced by the solution on output. The matrix $A is destroyed by the Householder transformations. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tridiag-solve(Math::Libgsl::Vector $diag, Math::Libgsl::Vector $abovediag where *.vector.size == $diag.vector.size - 1, Math::Libgsl::Vector $belowdiag where *.vector.size == $diag.vector.size - 1, Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size --> Math::Libgsl::Vector)

This function solves the general N-by-N system Ax = b where A is tridiagonal (N ≥ 2). The super-diagonal and sub-diagonal vectors must be one element shorter than the diagonal vector diag. The function returns the Math::Libgsl::Vector solution $x. In case of error a failure object is returned.

### tridiag-symm-solve(Math::Libgsl::Vector $diag, Math::Libgsl::Vector $offdiag where *.vector.size == $diag.vector.size - 1, Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size --> Math::Libgsl::Vector)

This function solves the general N-by-N system Ax = b where A is symmetric tridiagonal (N ≥ 2). The off-diagonal vector $offdiag must be one element shorter than the diagonal vector $diag. The function returns the Math::Libgsl::Vector solution $x. In case of error a failure object is returned.

### tridiag-cyc-solve(Math::Libgsl::Vector $diag where *.vector.size ≥ 3, Math::Libgsl::Vector $abovediag where *.vector.size == $diag.vector.size, Math::Libgsl::Vector $belowdiag where *.vector.size == $diag.vector.size, Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size --> Math::Libgsl::Vector)

This function solves the general N-by-N system Ax = b where A is cyclic tridiagonal (N ≥ 3). The cyclic super-diagonal and sub-diagonal vectors $abovediag and $belowdiag must have the same number of elements as the diagonal vector $diag. The function returns the Math::Libgsl::Vector solution $x. In case of error a failure object is returned.

### tridiag-symm-cyc-solve(Math::Libgsl::Vector $diag where *.vector.size ≥ 3, Math::Libgsl::Vector $offdiag where *.vector.size == $diag.vector.size, Math::Libgsl::Vector $rhs where *.vector.size == $diag.vector.size --> Math::Libgsl::Vector)

This function solves the general N-by-N system Ax = b where A is symmetric cyclic tridiagonal (N ≥ 3). The cyclic off-diagonal vector $offdiag must have the same number of elements as the diagonal vector $diag. The function returns the Math::Libgsl::Vector solution $x. In case of error a failure object is returned.

### tri-upper-rcond(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num)

This function estimates the reciprocal condition number. In case of error a failure object is returned.

### tri-lower-rcond(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num)

This function estimates the reciprocal condition number. In case of error a failure object is returned.

### tri-upper-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int)

This function computes the in-place inverse of the triangular matrix $T, stored in the upper triangle. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tri-lower-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int)

This function computes the in-place inverse of the triangular matrix $T, stored in the lower triangle. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tri-upper-unit-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int)

This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tri-lower-unit-invert(Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int)

This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tri-invert(Int $Uplo, Int $Diag, Math::Libgsl::Matrix $T where *.matrix.size1 == $T.matrix.size2 --> Int)

This function is available only from the C library v2.6. This function computes the in-place inverse of the triangular matrix $T, stored in the lower triangle when $Uplo = CblasLower and upper triangle when $Uplo = CblasUpper. The parameter $Diag = CblasUnit, CblasNonUnit specifies whether the matrix is unit triangular. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tri-LTL(Math::Libgsl::Matrix $L where *.matrix.size1 == $L.matrix.size2 --> Int)

This function is available only from the C library v2.6. This function computes the product LT(L) in-place and stores it in the lower triangle of $L on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tri-UL(Math::Libgsl::Matrix $LU where *.matrix.size1 == $LU.matrix.size2 --> Int)

This function is available only from the C library v2.6. This function compute the product $UL where U is upper triangular and L is unit lower triangular, stored in $LU, as computed by LU-decomp. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tri-rcond(Int $Uplo, Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Num)

This function is available only from the C library v2.6. This function estimates the 1-norm reciprocal condition number of the triangular matrix $A, using the lower triangle when $Uplo is CblasLower and upper triangle when $Uplo is CblasUpper. The function returns the reciprocal condition number 1/(||A||₁||A⁻¹||₁). In case of error a failure object is returned.

### cholesky-band-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int)

This function is available only from the C library v2.6. This function factorizes the symmetric, positive-definite square matrix $A into the Cholesky decomposition A = LT(L). The input matrix $A is given in symmetric banded format, and has dimensions N-by-(p + 1), where p is the lower bandwidth of the matrix. On output, the entries of $A are replaced by the entries of the matrix L in the same format. In addition, the lower right element of $A is used to store the matrix 1-norm, used later by cholesky-band-rcond() to calculate the reciprocal condition number. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-band-solve(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $LLT.matrix.size1 --> Math::Libgsl::Vector)

This function is available only from the C library v2.6. This function solves the symmetric banded system Ax = b using the Cholesky decomposition of A held in the matrix $LLT which must have been previously computed by cholesky-band-decomp. The function returns the Math::Libgsl::Vector solution $x. In case of error a failure object is returned.

### cholesky-band-svx(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $LLT.matrix.size1 --> Int)

This function is available only from the C library v2.6. This function solves the symmetric banded system Ax = b using the Cholesky decomposition of A held in the matrix $LLT which must have been previously computed by cholesky-band-decomp. On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-band-invert(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2 --> Math::Libgsl::Matrix)

This function is available only from the C library v2.6. This function computes the inverse of a symmetric banded matrix from its Cholesky decomposition $LLT, which must have been previously computed by cholesky-band-decomp. The function returns the inverse matrix as a Math::Libgsl::Matrix solution object. In case of error a failure object is returned.

### cholesky-band-unpack(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2 --> Math::Libgsl::Matrix)

This function is available only from the C library v2.6. This function unpacks the lower triangular Cholesky factor from $LLT, which returns as a Math::Libgsl::Matrix object. In case of error a failure object is returned.

### cholesky-band-rcond(Math::Libgsl::Matrix $LLT where *.matrix.size1 == $LLT.matrix.size2 --> Num)

This function is available only from the C library v2.6. This function estimates the reciprocal condition number (using the 1-norm) of the symmetric banded positive definite matrix A, using its Cholesky decomposition provided in $LLT. The reciprocal condition number estimate, defined as 1/(||A||₁ · ||A⁻¹||₁), is returned.

### ldlt-band-decomp(Math::Libgsl::Matrix $A where *.matrix.size1 == $A.matrix.size2 --> Int)

This function is available only from the C library v2.6. This function factorizes the symmetric, non-singular square matrix $A into the decomposition A = LDT(L). The input matrix $A is given in symmetric banded format, and has dimensions N-by-(p + 1), where p is the lower bandwidth of the matrix. On output, the entries of $A are replaced by the entries of the matrices D and L in the same format. If the matrix is singular then the decomposition will fail, returning the error code GSL_EDOM . This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### ldlt-band-solve(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Vector $b where *.vector.size == $LDLT.matrix.size1 --> Math::Libgsl::Vector)

This function is available only from the C library v2.6. This function solves the symmetric banded system Ax = b using the LDT(L) decomposition of A held in the matrix $LDLT which must have been previously computed by ldlt-band-decomp. The function returns the Math::Libgsl::Vector solution $x. In case of error a failure object is returned.

### ldlt-band-svx(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2, Math::Libgsl::Vector $x where *.vector.size == $LDLT.matrix.size1 --> Int)

This function is available only from the C library v2.6. This function solves the symmetric banded system Ax = b using the LDT(L) decomposition of A held in the matrix $LDLT which must have been previously computed by ldlt-band-decomp. On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### ldlt-band-unpack(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2 --> List)

This function is available only from the C library v2.6. This function unpacks the unit lower triangular factor L from $LDLT. This function returns a Math::Libgsl::Matrix object, which holds the lower triangular portion of the matrix L, and a Math::Libgsl::Vector object, which contains the diagonal matrix D. In case of error a failure object is returned.

### ldlt-band-rcond(Math::Libgsl::Matrix $LDLT where *.matrix.size1 == $LDLT.matrix.size2 --> Num)

This function is available only from the C library v2.6. This function estimates the reciprocal condition number (using the 1-norm) of the symmetric banded nonsingular matrix A, using its LDT(L) decomposition provided in $LDLT. The reciprocal condition number estimate, defined as 1/(||A||₁ · ||A⁻¹||₁), is returned.

### balance-matrix(Math::Libgsl::Matrix $A --> Math::Libgsl::Vector)

This function replaces the matrix $A with its balanced counterpart and returns the diagonal elements of the similarity transformation as a Math::Libgsl::Vector object. In case of error a failure object is returned.

Complex
-------

### LU-cdecomp(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> List)

This function factorizes the matrix A into the LU decomposition PA = LU. The factorization is done in place, so on output the diagonal and upper triangular (or trapezoidal) part of the input matrix A contain the matrix U. The lower triangular (or trapezoidal) part of the input matrix (excluding the diagonal) contains L. The diagonal elements of L are unity, and are not stored. The return value is a List: the sign of the permutation and a permutation object, which encodes the permutation matrix P. In case of error a failure object is returned.

### LU-csolve(Math::Libgsl::Matrix::Complex64 $LU where *.matrix.size1 == $LU.matrix.size2, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector::Complex64 $b where *.vector.size == $LU.matrix.size1 --> Math::Libgsl::Vector::Complex64)

This function solves the square system Ax = b using the LU decomposition of A into (LU, p) given by the output of LU-decomp. In case of error a failure object is returned.

### LU-csvx(Math::Libgsl::Matrix::Complex64 $LU where *.matrix.size1 == $LU.matrix.size2, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector::Complex64 $x where *.vector.size == $LU.matrix.size1 --> Int)

This function solves the square system Ax = b in-place using the precomputed LU decomposition of A into (LU, p). On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LU-crefine(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Matrix::Complex64 $LU where { $LU.matrix.size1 == $LU.matrix.size2 && $A.matrix.size1 == $LU.matrix.size2 }, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1, Math::Libgsl::Vector::Complex64 $b where *.vector.size == $LU.matrix.size1, Math::Libgsl::Vector::Complex64 $x where *.vector.size == $LU.matrix.size1 --> Int)

This function applies an iterative improvement to x, the solution of Ax = b, from the precomputed LU decomposition of A into (LU, p). This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LU-cinvert(Math::Libgsl::Matrix::Complex64 $LU, Math::Libgsl::Permutation $p where *.size == $LU.matrix.size1 --> Math::Libgsl::Matrix::Complex64)

This function computes the inverse of a matrix A from its LU decomposition (LU, p), returning the matrix inverse. In case of error a failure object is returned.

### LU-cdet(Math::Libgsl::Matrix::Complex64 $LU, Int $signum where * ~~ -1|1 --> Complex)

This function computes the determinant of a matrix A from its LU decomposition, $LU, and the sign of the permutation, $signum. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LU-clndet(Math::Libgsl::Matrix::Complex64 $LU --> Num)

This function computes the determinant the logarithm of the absolute value of the determinant of a matrix A, ln |det(A)|, from its LU decomposition, $LU. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### LU-csgndet(Math::Libgsl::Matrix::Complex64 $LU, Int $signum where * ~~ -1|1 --> Complex)

This function computes the sign or phase factor of the determinant of a matrix A, det(A)/|det(A)| from its LU decomposition, $LU. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-cdecomp(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> Int)

This function factorizes the symmetric, positive-definite square matrix $A into the Cholesky decomposition A = LT(L). On input, the values from the diagonal and lower-triangular part of the matrix $A are used (the upper triangular part is ignored). On output the diagonal and lower triangular part of the input matrix $A contain the matrix L, while the upper triangular part contains the original matrix. If the matrix is not positive-definite then the decomposition will fail, returning the error code GSL_EDOM. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-csolve(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector::Complex64 $b where *.vector.size == $A.matrix.size1 --> Math::Libgsl::Vector::Complex64)

These functions solve the system Ax = b using the Cholesky decomposition of $A which must have been previously computed by cholesky-cdecomp. The function returns the Math::Libgsl::Vector::Complex64 $x. In case of error a failure object is returned.

### cholesky-csvx(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector::Complex64 $x where *.vector.size == $A.matrix.size1 --> Int)

This function solves the system Ax = b in-place using the Cholesky decomposition of $A which must have been previously computed by cholesky-cdecomp. On input $x should contain the right-hand side b, which is replaced by the solution on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### cholesky-cinvert(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> Int)

These functions compute the inverse of a matrix from its Cholesky decomposition $A, which must have been previously computed by cholesky-cdecomp. On output, the inverse is stored in-place in $A. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### hermtd-cdecomp(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> Math::Libgsl::Vector::Complex64)

This function factorizes the hermitian matrix $A into the symmetric tridiagonal decomposition UTT(U). On output the diagonal and subdiagonal part of the input matrix $A contain the tridiagonal matrix T. The remaining lower triangular part of the input matrix contains the Householder vectors which, together with the Householder coefficients tau returned as a Math::Libgsl::Vector::Complex64 object, encode the orthogonal matrix U. The upper triangular part of $A and the imaginary parts of the diagonal are not referenced. In case of error a failure object is returned.

### hermtd-cunpack(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2, Math::Libgsl::Vector::Complex64 $tau where *.vector.size == $A.matrix.size1 - 1 --> List)

This function unpacks the encoded symmetric tridiagonal decomposition ($A, $tau) obtained from hermtd-cdecomp. The function returns a List: the unitary matrix $U as a Math::Libgsl::Matrix::Complex64, the real vector of diagonal elements $diag as a Math::Libgsl::Vector, and the real vector of subdiagonal elements $subdiag as a Math::Libgsl::Vector. In case of error a failure object is returned.

### hermtd-cunpack-T(Math::Libgsl::Matrix::Complex64 $A where *.matrix.size1 == $A.matrix.size2 --> List)

This function unpacks the diagonal and subdiagonal of the encoded symmetric tridiagonal decomposition ($A, $tau) obtained from hermtd-cdecomp. The function returns a List of two real Math::Libgsl::Vector: $diag and $subdiag. In case of error a failure object is returned.

### householder-ctransform(Math::Libgsl::Vector::Complex64 $w --> Complex)

This function prepares a Householder transformation H = I − τvT(v) which can be used to zero all the elements of the input vector $w except the first. On output the Householder vector v is stored in $w and the scalar τ is returned. The householder vector v is normalized so that v[0] = 1, however this 1 is not stored in the output vector. Instead, $w[0] is set to the first element of the transformed vector, so that if u = Hw, w[0] = u[0] on output and the remainder of u is zero. This function returns a Complex: $τ.

### householder-chm(Complex $tau, Math::Libgsl::Vector::Complex64 $v, Math::Libgsl::Matrix::Complex64 $A --> Int)

This function applies the Householder matrix H defined by the scalar $tau and the vector $v to the left-hand side of the matrix $A. On output the result HA is stored in $A. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### householder-cmh(Complex $tau, Math::Libgsl::Vector::Complex64 $v, Math::Libgsl::Matrix::Complex64 $A --> Int)

This function applies the Householder matrix H defined by the scalar $tau and the vector $v to the right-hand side of the matrix $A. On output the result AH is stored in $A. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### householder-chv(Complex $tau, Math::Libgsl::Vector::Complex64 $v, Math::Libgsl::Vector::Complex64 $w --> Int)

This function applies the Householder transformation H defined by the scalar $tau and the vector $v to the vector $w. On output the result Hw is stored in w. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tri-cinvert(Int $Uplo, Int $Diag, Math::Libgsl::Matrix::Complex64 $A --> Int)

This function is available only from the C library v2.6. This function computes the in-place inverse of the triangular matrix $T, stored in the lower triangle when $Uplo = CblasLower and upper triangle when $Uplo = CblasUpper. The parameter $Diag = CblasUnit, CblasNonUnit specifies whether the matrix is unit triangular. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tri-cLHL(Math::Libgsl::Matrix::Complex64 $L --> Int)

This function is available only from the C library v2.6. This function computes the product LT(L) in-place and stores it in the lower triangle of $L on output. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

### tri-cUL(Math::Libgsl::Matrix::Complex64 $LU --> Int)

This function is available only from the C library v2.6. This function compute the product $UL where U is upper triangular and L is unit lower triangular, stored in $LU, as computed by LU-cdecomp. This function returns GSL_SUCCESS if successful, or one of the error codes listed in Math::Libgsl::Constants::gsl-error.

C Library Documentation
=======================

For more details on libgsl see [https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/). The excellent C Library manual is available here [https://www.gnu.org/software/gsl/doc/html/index.html](https://www.gnu.org/software/gsl/doc/html/index.html), or here [https://www.gnu.org/software/gsl/doc/latex/gsl-ref.pdf](https://www.gnu.org/software/gsl/doc/latex/gsl-ref.pdf) in PDF format.

Prerequisites
=============

This module requires the libgsl library to be installed. Please follow the instructions below based on your platform:

Debian Linux and Ubuntu 20.04
-----------------------------

    sudo apt install libgsl23 libgsl-dev libgslcblas0

That command will install libgslcblas0 as well, since it's used by the GSL.

Ubuntu 18.04
------------

libgsl23 and libgslcblas0 have a missing symbol on Ubuntu 18.04. I solved the issue installing the Debian Buster version of those three libraries:

  * [http://http.us.debian.org/debian/pool/main/g/gsl/libgslcblas0_2.5+dfsg-6_amd64.deb](http://http.us.debian.org/debian/pool/main/g/gsl/libgslcblas0_2.5+dfsg-6_amd64.deb)

  * [http://http.us.debian.org/debian/pool/main/g/gsl/libgsl23_2.5+dfsg-6_amd64.deb](http://http.us.debian.org/debian/pool/main/g/gsl/libgsl23_2.5+dfsg-6_amd64.deb)

  * [http://http.us.debian.org/debian/pool/main/g/gsl/libgsl-dev_2.5+dfsg-6_amd64.deb](http://http.us.debian.org/debian/pool/main/g/gsl/libgsl-dev_2.5+dfsg-6_amd64.deb)

Installation
============

To install it using zef (a module management tool):

    $ zef install Math::Libgsl::LinearAlgebra

AUTHOR
======

Fernando Santagata <nando.santagata@gmail.com>

COPYRIGHT AND LICENSE
=====================

Copyright 2020 Fernando Santagata

This library is free software; you can redistribute it and/or modify it under the Artistic License 2.0.

