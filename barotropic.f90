
Module velocity
        implicit none

        contains

            subroutine ubar(dat,outdat,z_w,z_wu,Nroms,II,JJ,xi_rho,eta_rho)

            ! ----------------------------------
            ! Program : ubar
            !
            ! Trond Kristiansen, March 04 2009, bug fix 09.11.2009
            ! Rutgers University, NJ.
            ! -------------------------------------------------------------------------------------------------------
            !
            ! USAGE: Compile this routine using Intel Fortran compiler and create
            ! a python module using the command:
            ! f2py --verbose --fcompiler=intelem -c -m barotropic barotropic.f90
            ! or
            ! f2py --verbose c -m barotropic barotropic.f90
            ! The resulting module is imported to python using:
            ! import barotropic
            ! To call the function from python use:
            ! barotropic.ubar(dat,outdat,zr,zw,Nroms,II,JJ)
            !
            ! where: dat is the data such as temperature (3D structure (z,y,x))
            !
            !        outdat is a 3D output array with the correct size (Nroms,JJ,II)
            !        zr is the depth matrix for the output grid (Nroms,JJ,II)
            !        zs is the 1D SODA depth-matrix (e.g. zs=[5,10,20,30])
            !        Nroms is the total depth levels in output grid
            !        JJ is the total grid points in eta direction
            !        II is the total grid points in xi direction
            ! -------------------------------------------------------------------------------------------------------

            integer eta_rho, xi_rho, II, JJ, ic, jc, kc, Nroms
            REAL(4), dimension(Nroms,JJ,II) :: dat
            REAL(4), dimension(JJ,II) :: outdat
            REAL(4), dimension(Nroms+1,eta_rho,xi_rho) ::  z_w
            REAL(4), dimension(Nroms+1,JJ,II) ::  z_wu

!f2py intent(in,overwrite) dat,z_w,z_wu,Nroms, JJ, II, xi_rho, eta_rho
!f2py intent(in,out,overwrite) outdat
!f2py intent(hide) ic,jc,kc

            print*,'--->Started ubar calculations'
            ! average z_w to Arakawa-C u,v-points (z_wu, z_wv)
            do jc=1,JJ
              do ic=2,II+1
                  do kc=1,Nroms+1
                    z_wu(kc,jc,ic-1) = 0.5*(z_w(kc,jc,ic-1)+z_w(kc,jc,ic))
                  end do
               end do
            end do

            do jc=1,JJ
              do ic=1,II
                 outdat(jc,ic)=0.0
                  do kc=1,Nroms
                        outdat(jc,ic) = outdat(jc,ic) + dat(kc,jc,ic)*abs(z_wu(kc+1,jc,ic) - z_wu(kc,jc,ic))
                  end do
                  if (abs(z_wu(1,jc,ic)) > 0.0) then
                    outdat(jc,ic) = outdat(jc,ic)/abs(z_wu(1,jc,ic))
                  else
                    outdat(jc,ic) = 0.0
                  end if
              end do
            end do

            end subroutine ubar

            subroutine vbar(dat,outdat,z_w,z_wv,Nroms,II,JJ,xi_rho,eta_rho)

            ! ----------------------------------
            ! Program : vbar
            !

            !
            ! Trond Kristiansen, March 04 2009, bug fix 09.11.2009
            ! Rutgers University, NJ.
            ! -------------------------------------------------------------------------------------------------------
            !
            ! USAGE: Compile this routine using Intel Fortran compiler and create
            ! a python module using the command:
            ! f2py --verbose --fcompiler=intelem -c -m barotropic barotropic.f90
            ! or
            ! f2py --verbose --fcompiler=intelem  -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m barotropic barotropic.f90
            ! The resulting module is imported to python using:
            ! import barotropic
            ! To call the function from python use:
            ! barotropic.ubar(dat,bathymetry,outdat,zr,zw,Nroms,II,JJ)
            !
            ! where: dat is the data such as temperature (3D structure (z,y,x))
            !
            !        outdat is a 3D output array with the correct size (Nroms,JJ,II)
            !        zr is the depth matrix for the output grid (Nroms,JJ,II)
            !        zs is the 1D SODA depth-matrix (e.g. zs=[5,10,20,30])
            !        Nroms is the total depth levels in output grid
            !        JJ is the total grid points in eta direction
            !        II is the total grid points in xi direction
            ! -------------------------------------------------------------------------------------------------------

            integer eta_rho, xi_rho, II, JJ, ic, jc, kc, Nroms
            REAL(4), dimension(Nroms,JJ,II) :: dat
            REAL(4), dimension(JJ,II) :: outdat
            REAL(4), dimension(Nroms+1,eta_rho,xi_rho) ::    z_w
            REAL(4), dimension(Nroms+1,JJ,II) ::  z_wv

!f2py intent(in,overwrite) dat, z_w, z_wv, Nroms, JJ, II, xi_rho, eta_rho
!f2py intent(in,out,overwrite) outdat
!f2py intent(hide) ic,jc,kc

            print*,'--->Started vbar calculations'
            do jc=2,JJ+1
              do ic=1,II
                  do kc=1,Nroms+1
                    z_wv(kc,jc-1,ic) = 0.5*(z_w(kc,jc-1,ic)+z_w(kc,jc,ic))
                  end do
               end do
            end do

            do jc=1,JJ
              do ic=1,II

                 outdat(jc,ic)=0.0

                  do kc=1,Nroms
                        outdat(jc,ic) = outdat(jc,ic) + dat(kc,jc,ic)*abs(z_wv(kc+1,jc,ic) - z_wv(kc,jc,ic))
                  end do

                  if (abs(z_wv(1,jc,ic)) > 0.0) then
                    outdat(jc,ic) = outdat(jc,ic)/abs(z_wv(1,jc,ic))
                  else
                    outdat(jc,ic) = 0.0
                  end if

              end do
            end do

            end subroutine vbar

     end module velocity