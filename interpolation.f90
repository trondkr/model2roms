       
Module interpolation
        implicit none
                
        contains
  
            subroutine doVertInter(outdat,dat,bathymetry,zr,zs,Nroms,Nsoda,II,JJ,xi_rho,eta_rho)
            
            ! ----------------------------------
            ! Program : doVertInter
            !
            ! This routine interpolates from z-levels to sigma levels using linear interpolation.
            !
            ! The index values in python goes from 0 toN while in Fortran they run from 1 to N+1. This is important to
            ! remember if one wants to compare input index wtih output index in fortran and python.
            !
            ! This routine assumes that the two depth matrixes zr (ROMS) and zs (SODA) are arranged from shallowest
            ! (index=1) to deepest (index=N+1). The depth matrizes must also be negative (if positive, reverse all
            ! comparison signs (e.g. LT, GT) in the program or multiply with minus 1). The input data are arranged with
            ! deepest values at highest index (N+1 e.g. dat(N+1)==bottom, dat(1)==surface). This is done so because
            ! it is the way SODA data are organized (bottom at highest index). However, ROMS output files are organized vice versa, so
            ! to accomodate that the output values are stored according to the ROMS structure. Highest index (N+1) equals surface,
            ! while lowest index equals bottom (index=1)(see how outdat(kc,jc,ic) is used opposite of the loop over kc).
            !
            ! Trond Kristiansen, December 2008, January, and March 2009
            ! Rutgers University, NJ.
            ! -------------------------------------------------------------------------------------------------------
            !
            ! USAGE: Compile this routine using Intel Fortran compiler and create
            ! a python module using the command:
            !
            ! For Python 3 also add the LDFLAG:
            ! https://github.com/conda-forge/f90wrap-feedstock/issues/1
            ! export LDFLAGS="-undefined dynamic_lookup -bundle"
            ! f2py --verbose --fcompiler=intelem -c -m interpolation interpolation.f90
            !
            ! The resulting module is imported to python using:
            ! import vertInterp as interp2D
            ! To call the function from python use:
            ! interp2D.doHorInterpolation(dat,bathymetry,outdat,zr,zs,Nroms,Nsoda,II,JJ)
            !
            ! where: dat is the data such as temperature (3D structure (z,y,x))
            !        bathymetry is the 2D bottom matrix from the output grid (in ROMS this is usually 'h')
            !        outdat is a 3D output array with the correct size (Nroms,JJ,II)
            !        zr is the depth matrix for the output grid (Nroms,JJ,II)
            !        zs is the 1D SODA depth-matrix (e.g. zs=[5,10,20,30])
            !        Nroms is the total depth levels in output grid
            !        JJ is the total grid points in eta direction
            !        II is the total grid points in xi direction
            ! -------------------------------------------------------------------------------------------------------
            
            REAL(4) rz2, rz1, fill
            integer eta_rho, xi_rho, II, JJ, ic, jc, kc, kT, kkT, Nsoda, Nroms
            REAL(4), dimension(Nsoda,JJ,II) :: dat
            REAL(4), dimension(eta_rho,xi_rho) :: bathymetry
            REAL(4), dimension(Nroms,JJ,II) :: outdat
            REAL(4), dimension(Nsoda) ::  zs
            REAL(4), dimension(Nroms,eta_rho,xi_rho) :: zr

!f2py intent(in,out,overwrite) outdat       
!f2py intent(in,overwrite) dat, bathymetry, zr, zs
!f2py intent(in,overwrite) Nroms, Nsoda, JJ, II, xi_rho, eta_rho
!f2py intent(hide) ic,jc,kc,kT,rz1,rz2, kkT
            fill=-10000
            do jc=1,JJ
              do ic=1,II
                  do kc=1,Nroms
                      ! CASE 1: ROMS deeper than SODA. This part searches for deepest good value if ROMS depth is deeper
                      ! than SODA. This means that if no value, or only fill_value, is available from SODA where ROMS is
                      ! deepest, the closest value from SODA is found by looping upward in the water column.
                      IF (zr(kc,jc,ic) .LT. zs(Nsoda)) THEN
                          outdat(kc,jc,ic)=dat(Nsoda,jc,ic)
                          if (MAXVAL(dat(:,jc,ic)) .GT. fill) then
                            if (dat(Nsoda,jc,ic) .LT. fill) then
                              !print*,'Inside dovert and finding deepest depth with good values. current',dat(Nsoda,jc,ic)
                              DO kT=1,Nsoda
                                if (dat(Nsoda-kT,jc,ic) .GT. fill) then
                                    print*,'working on depth',kT,'with value',dat(kT,jc,ic)
                                    outdat(kc,jc,ic)=dat(Nsoda-kT,jc,ic)
                                    print*,'Able to find good value at depth ', Nsoda-kT
                                    exit
                                end if
                              end do
                             end if
                            end if
                          !print*,zr(kc,jc,ic),zs(Nsoda),dat(Nsoda,jc,ic),jc,ic,'case 1'
                    
                      ! CASE 2: ROMS shallower than SODA
                      ELSE IF (zr(kc,jc,ic) .GT. zs(1)) THEN
                          outdat(kc,jc,ic)=dat(1,jc,ic)
                     
                      ELSE
                          ! DO LOOP BETWEEN SURFACE AND BOTTOM
                          DO kT=1,Nsoda
                              ! CASE 3: ROMS deeper than SODA for one layer, but shallower than next SODA layer (bottom in between)
                              ! Deeper than some SODA depth layer, but shallower than next layer which is below bottom
                              IF (zr(kc,jc,ic) .LE. zs(kT) .AND.               &
                                -(bathymetry(jc,ic)) .GT. zs(kT+1)) THEN
                                outdat(kc,jc,ic)=dat(kT,jc,ic)
                                
                                ! We do not want to give the deepest depth a fill_value, so we find
                                ! the closest value to deepest depth.
                                if (MAXVAL(dat(:,jc,ic)) .GT. fill) then
                               
                                    if (dat(kT,jc,ic) .LE. fill) then
                                       !print*,'case3:Need to find better value for depth ',kT,'which has value ',dat(kT,jc,ic)
                                        DO kkT=1,Nsoda
                                            if (dat(kT-kkT,jc,ic) .GT. fill) then
                                                 outdat(kc,jc,ic)=dat(kT-kkT,jc,ic)
                            
                                                exit
                                            end if
                                        end do
                                    end if
                                end if
                                
                                ! CASE 4: Special case where ROMS layers are much deeper than in SODA
                                ELSE IF (zr(kc,jc,ic) .LE. zs(kT) .AND. dat(kT,jc,ic) .GT. fill &
                                .AND. dat(kT+1,jc,ic) .LE. fill) THEN
                                outdat(kc,jc,ic)=dat(kT,jc,ic)
                              
                              
                              ! CASE 5: ROMS layer in between two SODA layers
                              ! This is the typical case for most layers
                              ELSE IF ( (zr(kc,jc,ic) .LE. zs(kT)) .AND.       &
                              (zr(kc,jc,ic) .GE. zs(kT+1)) .AND.               &
                              (-bathymetry(jc,ic) .LE. zs(kT+1)) ) THEN
                              
                                 rz2 = abs((zr(kc,jc,ic)-zs(kT+1))/            &
                                 (abs(zs(kT+1))-abs(zs(kT))))
                                 
                                 rz1 = 1.0-rz2
        
                                 outdat(kc,jc,ic) = (rz1*dat(kT+1,jc,ic) &
                                 + rz2*dat(kT,jc,ic))
                                
                                if (MAXVAL(dat(:,jc,ic)) .GT. fill) then
                               
                                    if (dat(kT,jc,ic) .LE. fill .OR. dat(kT+1,jc,ic) .LE. fill) then
                                       !print*,'case4:Need to find better value for depth ',kT,kT+1,'which has &
                                   !    values ',dat(kT,jc,ic),dat(kT+1,jc,ic)
                                        DO kkT=1,Nsoda
                                            if (dat(kT-kkT,jc,ic) .GT. fill .and. dat(kT-kkT+1,jc,ic) .GT. fill  ) then
                                                 !print*,'CASE 4: Found good value at depth',kT-kkT,kt-kkT+1
                                                 !print*,'with values',dat(kT-kkT,jc,ic), dat(kt-kkT+1,jc,ic)
                                                 outdat(kc,jc,ic) = (rz1*dat(kT+1-kkT,jc,ic) &
                                                 + rz2*dat(kT-kkT,jc,ic))
                            
                                                exit
                                            end if
                                        END DO
                                    end if
                                 end if
                                
                                 exit
                                 
                              END IF
                          ! DO LOOP BETWEEN SURFACE AND BOTTOM: CASE 3,4,5
                          END DO
                      ! TEST ALL CASES IF LOOP: CASE 1,2,3,4,5
                      END IF
                     
                  end do
              end do
            end do
            
        
      end subroutine doVertInter
      
      subroutine rho2u(udata,rhodata,II,JJ,KK)
       
            ! ----------------------------------
            ! Program : rho2u
            !
            ! This routine interpolates RHO points to U points using simple linear interpolation
            ! The input matrix (rhodata) is a matrix of size (JJ,II). The output matrix is the
            ! interpolated RHO values at U points with dimensions (JJ,II-1).
            ! Trond Kristiansen, January 2009
            ! Rutgers University, NJ.
            ! -------------------------------------------------------------------------------------------------------
          
           
            integer KK, II, JJ, kc, ic, jc, fill
            REAL(4), dimension(KK,JJ,II) :: rhodata
            REAL(4), dimension(KK,JJ,II-1) :: udata
       
!f2py intent(in,out,overwrite) udata
!f2py intent(in,overwrite) rhodata, KK, JJ, II
!f2py intent(hide) ic,jc,kc, fill

            fill=10000
            print*,'---> Started horisontal rho2u interpolation'
            do kc=1,KK
                do jc=1,JJ
                    do ic=2,II-1
                        
                        udata(kc,jc,1)=rhodata(kc,jc,1)
                        ! Now make sure that if we have two stations where one has values
                        ! and the other not, we only use the good value
                        ! case 1: one value is good (jc+1) other bad (jc-1)
                        if (abs(rhodata(kc,jc,ic-1)) > fill .AND. abs(rhodata(kc,jc,ic+1)) < fill) then
                            udata(kc,jc,ic)=(rhodata(kc,jc,ic+1))
                        ! case 2: one value is good (jc-1) other bad (jc+1)
                        else if (abs(rhodata(kc,jc,ic-1)) < fill .AND. abs(rhodata(kc,jc,ic+1)) > fill) then
                            udata(kc,jc,ic)=(rhodata(kc,jc,ic-1))
                        ! Both values are bad:
                        else if (abs(rhodata(kc,jc,ic-1)) > fill .AND. abs(rhodata(kc,jc,ic+1)) > fill) then
                            udata(kc,jc,ic)=0.0
                        ! Both values are good and we do linear interpolation
                        else 
                            udata(kc,jc,ic)=(rhodata(kc,jc,ic-1)+rhodata(kc,jc,ic+1))*0.5
                        end if
                    end do
                end do
            end do
            print*,'-----> Ended horisontal rho2u interpolation'
        end subroutine rho2u
            
      subroutine rho2v(vdata,rhodata,II,JJ,KK)
       
            ! ----------------------------------
            ! Program : rho2v
            !
            ! This routine interpolates RHO points to V points using simple linear interpolation
            ! The input matrix (rhodata) is a matrix of size (JJ,II). The output matrix is the
            ! interpolated RHO values at U points with dimensions (JJ-1,II).
            ! Trond Kristiansen, January, February, and March2009
            ! Rutgers University, NJ.
            ! -------------------------------------------------------------------------------------------------------
          
           integer KK, II, JJ, kc, ic, jc, fill
           REAL(4), dimension(KK,JJ,II) :: rhodata
           REAL(4), dimension(KK,JJ-1,II) :: vdata

!f2py intent(in,out,overwrite) vdata
!f2py intent(in,overwrite) rhodata, KK, JJ, II
!f2py intent(hide) ic,jc,kc, fill
            
            fill=10000
            print*,'---> Started horisontal rho2v interpolation'
            do kc=1,KK
                do jc=2,JJ-1
                    do ic=1,II
                       
                        
                        vdata(kc,1,ic)=rhodata(kc,1,ic)
                            
                        if (abs(rhodata(kc,jc-1,ic)) > fill .AND. abs(rhodata(kc,jc+1,ic)) < fill) then
                            vdata(kc,jc,ic)=(rhodata(kc,jc+1,ic))
                        else if (abs(rhodata(kc,jc-1,ic)) < fill .AND. abs(rhodata(kc,jc+1,ic)) > fill) then
                            vdata(kc,jc,ic)=(rhodata(kc,jc-1,ic))
                        else if (abs(rhodata(kc,jc-1,ic)) > fill .AND. abs(rhodata(kc,jc+1,ic)) > fill) then
                            vdata(kc,jc,ic)=0.0
                        else
                            vdata(kc,jc,ic)=(rhodata(kc,jc-1,ic)+rhodata(kc,jc+1,ic))*0.5
                        end if
                    end do
                end do
            end do
            print*,'-----> Ended horisontal rho2v interpolation'
        end subroutine rho2v
            
   
        subroutine rotate(urot,vrot,u_rho,v_rho,angle,II,JJ,KK)
            ! ----------------------------------
            ! Program : rotate
            !
            ! This routine rotates u and v velocities in the North-South grid to an
            ! the output North-South grid with angle "angle"
            ! Trond Kristiansen, January 2009
            ! Rutgers University, NJ.
            ! -------------------------------------------------------------------------------------------------------
            
           REAL(4), dimension(KK,JJ,II) :: urot, vrot
           REAL(4), dimension(KK,JJ,II) :: u_rho, v_rho
           REAL(4), dimension(JJ,II)  :: angle
           integer KK, II, JJ, kc, ic, jc
          
!f2py intent(in,out,overwrite) urot, vrot
!f2py intent(in,overwrite)  u_rho, v_rho, angle, KK, JJ, II
!f2py intent(hide) ic,jc,kc
    
           print*,'---> Started rotation of velocities'
           do kc=1,KK
             do jc=1,JJ
                do ic=1,II
                    
                    ! https://www.myroms.org/forum/viewtopic.php?f=3&t=295

                    urot(kc,jc,ic)=u_rho(kc,jc,ic)*COS(angle(jc,ic)) + v_rho(kc,jc,ic)*SIN(angle(jc,ic))
                    vrot(kc,jc,ic)=-u_rho(kc,jc,ic)*SIN(angle(jc,ic)) + v_rho(kc,jc,ic)*COS(angle(jc,ic))
                     
                    !print*, vrot(kc,jc,ic), urot(kc,jc,ic), kc,jc,ic !, sin(angle(jc,ic)), cos(angle(jc,ic))
                    !print*, v_rho(kc,jc,ic), u_rho(kc,jc,ic), ic,jc
                end do
             end do
            end do
            print*,'-----> Ended rotation of velocities'
        end subroutine rotate
    
     end module interpolation