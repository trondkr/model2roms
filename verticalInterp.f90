       
Module interpolation
        implicit none
                
        contains
  
            subroutine doVertInter(dat,outdat,zr,zs,Nroms,Nsoda,II,JJ)
            
            ! ----------------------------------
            ! Program : DOVERTINTER
            !
            ! This routine interpolates from z-levels to sigma levels using linear interpolation.
            ! This is a simplified approach and should be changed to bilinear interpolation.
            ! The input values to this routine reflects the python index form (from 0 to N), and
            ! for that reason the counters in the Fortran routine adds 1 to each counter loop since Fortran
            ! counts from 1 to N+1.This means that the do loop in python would be one less than in fortran.
            !
            ! Trond Kristiansen, December 2008
            ! Rutgers University, NJ.
            ! ----------------------------------
            !
            ! USAGE: Compile this routine using Intel Fortran compiler and create
            ! a python module using the command:
            ! f2py --verbose --fcompiler=intel -c -m vertInterp verticalInterp.f90
            !
            ! The resulting module is imported to python using:
            ! import vertInterp as interp2D
            ! To call the function from python use:
            ! interp2D.doHorInterpolation(..args..)
            !
            ! -----------------------------------
            
            double precision rz2, rz1
            integer II, JJ, ic, jc, kc, kT, Nsoda, Nroms
            double precision, dimension(Nsoda,JJ,II) :: dat
            double precision, dimension(Nroms,JJ,II) :: zr, outdat
            double precision, dimension(Nsoda) ::  zs
       
!f2py intent(in) dat, zr, zs, Nroms, Nsoda, JJ, II
!f2py intent(in,out) outdat
!f2py intent(hide) ic,jc,kc,kT,rz1,rz2
            
            do jc=1,JJ
              do ic=1,II
                  do kc=1,Nroms
            
                      ! If deepest depth in new grid is deeper than
                      ! any depth in the old grid, use the values found
                      ! at the deepest depth in the old grid in the new grid.
                      IF (zr(kc,jc,ic) .LT. zs(Nsoda)) THEN
                          outdat(kc,jc,ic)=dat(Nsoda,jc,ic)
                          !print*,zr(kc,jc,ic),zs(Nsoda),'deepest'
                      ! If the shallowest depth in the new grid is shallower than
                      ! any depth in the old grid, use the shallowest values.
                      ELSEIF (zr(kc,jc,ic) .GT. zs(1)) THEN
                          outdat(kc,jc,ic)=dat(1,jc,ic)
                          !print*,zr(kc,jc,ic),zs(1),'shallowest'
                      ! If none of the above, do linear interpolation
                      ELSE
                          DO kT=1,Nsoda+1
                             
                              IF ( (zr(kc,jc,ic) .LT. zs(kT)) .AND.            &
                              (zr(kc,jc,ic) .GE. zs(kT+1)) ) THEN
                                 
                                  rz1 = abs(zr(kc,jc,ic)-zs(kT))/                 &
                                   abs(zs(kT+1)-zs(kT))
                                  
                                  rz2 = 1.0-rz1 
        
                                  outdat(kc,jc,ic) = rz1*dat(kT,jc,ic)            &
                                   + rz2*dat(kT+1,jc,ic)
                                  
                                  EXIT
                              end if
                          end do
                      end if
           
                  end do
              end do
            end do
        
      !  print*,outdat(:,300,300)
      !  print*,shape(outdat)
        
      end subroutine doVertInter
    end module interpolation
    
