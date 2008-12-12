    
   
Module interpolation
        
        public :: doVertInter
        
        contains
  
            subroutine doVertInter(dat,zr,zs,Nroms,Nsoda,II,JJ,KK)
           
            double precision rz2, rz1
            integer II, JJ, KK, ic, jc, kc, kT, Nsoda, Nroms
            double precision, dimension(KK,JJ,II) :: tmp, zr, dat    
            double precision, dimension(Nsoda) ::  zs
            print*,'Hello from fortran'
!cf2py intent(in) :: dat, zr, zs, Nroms, Nsoda 
!cf2py intent(out) :: tmp
!cf2py intent(hide) :: ic,jc,kc,kT,rz1,rz2,JJ,KK,II
     
            do jc=1,JJ-1
              do ic=1,II-1
                  do kc=1,KK-1
            
                       
                      IF (zr(kc,jc,ic) .LT. abs(zs(0))) THEN
                          tmp(kc,jc,ic)=dat(Nsoda,jc,ic)
                          
                      ELSEIF (abs(zr(kc,jc,ic)) .GT. abs(zs(Nsoda-1))) THEN
                          tmp(kc,jc,ic)=dat(Nsoda-1,jc,ic)
                         
                      ELSE
                          DO kT=1,Nsoda-1
                            
                              IF ( (zr(kc,jc,ic).LT.zs(kT+1)).AND.(zr(kc,jc,ic)    &
                               .GE. zs(kT)) ) THEN
                              
                                 
                                  rz2 = (zr(kc,jc,ic)-zs(kT))/                  &
                                   (zs(kT+1)-zs(kT)) 
                                  rz1 = 1.0-rz2 
                                  tmp(kc,jc,ic) = rz1*dat(kT,jc,ic)              &
                                   + rz2*dat(kT+1,jc,ic)
                                  exit
                              end if
                          end do
                      end if
                  end do
              end do
            end do
        
      end subroutine doVertInter
    end module interpolation
    