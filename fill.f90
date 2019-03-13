Module extrapolate

implicit none

contains

SUBROUTINE fill(i1,i2,j1,j2,tx,critx,cor,mxs,za,nx,ny)

    ! ** This routine originates from the USA.
    ! ** It was delivered (by 3.5" diskette)
    ! ** to DNMI in nov. 1990.

    ! ***** Re-edited at DNMI in January '91 by H.Engedahl.
    ! ***** Again re-edited at DNMI in March '92 by H.Engedahl to fit the
    ! ***** ECOM3D output configuration.
    ! ***** DNMI/FoU 25.08.1993 A.Foss : Recoded for speedup.

    !// Solves Laplace's equation with Neumann boundary conditions
    !// (dA/dn = 0) in RECTANGULAR coordinates by an iterative method to
    !// fill in reasonable values at gridpoints containing values like
    !// "undef".

    !// NOTE: It is impossibel to make a really parallel (MPP) version of this
    !//       routine. One MPP node should work on a complete field.
    !//       Parallelization is done on a 'higher' level.  (A.Foss '98)

    !-----------------------------------------------------------------------
    !//     nx   = Array x dimension
    !//     ny   = Array y dimension
    !//   i1,i2  = Subarea in x direction to be filled
    !//   j1,j2  = Subarea in y direction to be filled
    !//     za   = Array to be filled (REAL)
    !//     tx   = All values in array A GREATER than Tx are filled (REAL)
    !//    critx = Criteria for relaxation, DEL**2 = CRIT
    !//            (Usually 4 orders of magnitude DOWN from data in A)
    !//     cor  = Coef. of overrelaxation, between +1.2 and +2.0
    !//     mxs  = Max. allowed no. of scans in relaxation procedure.
    !//   rmask  = Work array
    !//   error  = Array containing the errors in the relaxation procedure.
    !//   nvalue = No. of gridpoints with value (possibly 0) ... output
    !-----------------------------------------------------------------------
    !
    ! Converted to use as Python module with F2PY
    ! by Trond Kristiansen, 31.05.2012. Notice that the indexes (j,i) are
    ! flipped compared to original Fortran version.
    !
    ! f2py --verbose -c -m extrapolate fill.f90

    IMPLICIT NONE

    INTEGER  :: nx, ny, i1, i2, j1, j2, mxs
    REAL     :: tx,critx,cor

    REAL, dimension(ny,nx):: rmask, error, za

    INTEGER :: nvalue
    INTEGER              :: n, j, i, i1p1, i2m1, j1p1, j2m1, nnn, nbad
    REAL                 :: suma, asuma, crit, crtest

!f2py intent(in,out,overwrite)  za
!f2py intent(in) nx, ny, i1, i2, j1, j2, tx, critx, cor, mxs, nvalue
   
    n = 0
    suma = 0.

    DO j=j1,j2
    DO i=i1,i2
      IF (za(j,i) < tx) THEN
        suma = suma + za(j,i)
        n = n+1
      END IF
    END DO
    END DO

    nvalue = n

    IF (n < 1) THEN
  !    PRINT *
  !    PRINT *,"******************  WARNING  *******************"
  !    PRINT *,"SUBROUTINE FILL : NO USEFUL DATA IN THE FIELD"
  !    PRINT *,"ALL DATA IN THE INTERIOR (SUBAREA) WERE > ",tx
  !    PRINT *,"******************  WARNING  *******************"
  !    PRINT *
      RETURN
    END IF

    suma = suma/n
    asuma = 0.

    DO j=j1,j2
    DO i=i1,i2
      IF (za(j,i) < tx) THEN
        asuma = asuma + ABS(za(j,i)-suma)
        rmask(j,i) = 0.
      ELSE
        za(j,i) = suma
        rmask(j,i) = 1.
      END IF
    END DO
    END DO

    asuma = asuma/n

    crit = critx*asuma  ! af:  Using asuma not suma

    !// The value of suma ( i.e. the MEAN value of the array A ) is filled in
    !// the array A at all points with an initial value >= t.
    !// za(j,i) = suma may be regarded as the "first guess" in the iterative
    !// method.

    i1p1 = i1 + 1
    i2m1 = i2 - 1
    j1p1 = j1 + 1
    j2m1 = j2 - 1

    DO j=j1p1,j2m1
    DO i=i1p1,i2m1
      rmask(j,i) = cor*rmask(j,i)
    END DO
    END DO

    DO nnn=1,mxs

      DO j=j1p1,j2m1
        DO i=i1p1,i2m1
          error(j,i) = (za(j,i+1)+za(j+1,i)+za(j-1,i))*0.25-za(j,i)
        END DO
        DO i=i1p1,i2m1
          error(j,i) = error(j,i) + za(j,i-1)*0.25
          za(j,i) = za(j,i) + error(j,i)*rmask(j,i)
        END DO
      END DO

      !// Test convergence now and then (slow test loop)
      IF (nnn < mxs-5 .AND. MOD(nnn,10) == 0) THEN
        crtest = crit*cor
        nbad = 0
        j = j1
        DO WHILE (nbad == 0 .AND. j < j2m1)
          j=j+1
          DO i=i1p1,i2m1
            IF (ABS(error(j,i))*rmask(j,i) > crtest) nbad = 1
          END DO
        END DO
        IF (nbad == 0) GOTO 300
      END IF

      DO j=j1p1,j2m1
        za(j,i1) = za(j,i1) + (za(j,i1p1)-za(j,i1))*rmask(j,i1)
        za(j,i2) = za(j,i2) + (za(j,i2m1)-za(j,i2))*rmask(j,i2)
      END DO
      DO i=i1,i2
         za(j1,i) = za(j1,i) + (za(j1p1,i)-za(j1,i))*rmask(j1,i)
         za(j2,i) = za(j2,i) + (za(j2m1,i)-za(j2,i))*rmask(j2,i)
      END DO

    ENDDO

    300 CONTINUE

    RETURN

END SUBROUTINE

END MODULE
