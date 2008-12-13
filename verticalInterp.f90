    
!verticalInterp.f90   
Module interpolation
        implicit none
                
        contains
  
            subroutine doVertInter(dat,outdat,zr,zs,Nroms,Nsoda,II,JJ)
           
            double precision rz2, rz1
            integer II, JJ, ic, jc, kc, kT, Nsoda, Nroms
            double precision, dimension(Nsoda,JJ,II) :: dat
            double precision, dimension(Nroms,JJ,II) :: zr, outdat
            double precision, dimension(Nsoda) ::  zs
       
!cf2py intent(in) dat, zr, zs, Nroms, Nsoda
!cf2py intent(out) outdat
!cf2py intent(hide) ic,jc,kc,kT,rz1,rz2,JJ,II
           
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
                          DO kT=1,Nsoda
                             
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
        
        print*,outdat(:,200,200)
        print*,shape(outdat)
        
      end subroutine doVertInter
    end module interpolation
    