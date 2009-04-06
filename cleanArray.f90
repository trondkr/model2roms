        
Module cleanArray

        implicit none
                
        contains
  
 
 
 subroutine sweep(FILL,KK,datain,dataout,mask,II,JJ)
       
            ! ------------------------------------------------------------------
            ! Program : sweep
            !
            ! This routine finds points in the interpolated fields that was masked by
            ! the original land-array, but which is no longer masked in the new land-ocean array.
            ! Sometimes this leads to ocean posint in the new grid that are still
            ! ddefined as land because of the prevsious masking. This routine sweeps through
            ! the array and replaces fill values over ocean with the nearest points values.
            !
            ! Trond Kristiansen, March 2009
            ! Rutgers University, NJ.
            !
            ! Compile with:
            ! f2py --verbose -DF2PY_REPORT_ON_ARRAY_COPY=1 --fcompiler=intel -c -m cl cleanArray.f90
            !
            ! USAGE:
            ! import cl
            ! Zin = cl.cleanArray.sweep(Zin,Zout, maskarray,i,j,k)
            !
            !  (I-1,J+1)       (I,J+1)        (I+1,J+1)                                                       
            !   O----------------O----------------O                                                                        
            !   |                |                |                                                
            !   |                |                |                                                  
            !  (I-1,J)        X (I,J)          (I+1,J)
            !   |                |                |   
            !   |                |                |                                                  
            !   O----------------O----------------O                                                  
            !  (I-1,J-1)       (I,J-1)        (I+1,J-1)                                       
            !
    
            ! ------------------------------------------------------------------
          
           
            integer II, JJ, KK, ic, jc, jcc, icc, fill, kc, ff, total

            double precision, dimension(JJ,II) :: dataout,datain
            double precision, dimension(JJ,II) :: mask
            double precision, dimension(2,3)  :: foundWeight, foundData
            double precision counter, point, hori, vert
            
!f2py intent(in,out,overwrite) dataout
!f2py intent(in,overwrite) datain, mask
!f2py intent(in) II, JJ, FILL, KK
!f2py intent(hide) ic, jc,kc, ff, icc, jcc, counter, point, total, hori, verti
!f2py intent(hide) foundWeight, foundData

            
            print*,'---> Started cleaning of array'
            !print*,'---> Using a radius of ', KK,'points to fill in gaps'
            total   = 0    
            
            do jc=1,JJ
                do ic=1,II
                
                    counter = 0
                    point   = 0.0
                    
                    if (abs(datain(jc,ic)) > FILL .AND. mask(jc,ic)==1) THEN
                        ! Found one position that needs to be filled
                       ! print*,'Found bad values at:', jc,ic
                       ! print*, 'values:',abs(datain(jc,ic)), mask(jc,ic)
                        
                        total=total+1
                        foundWeight(:,:)=-99
                        foundData(:,:)=-99
                        ! Horisontally find values around problem value
                        do ff=-1,1,1
                        do kc=-KK,KK
                        
                            hori=ic+kc
                            
                            if (hori .LE. 1 ) then
                                hori=1
                            else if (hori .GE. II) then
                                hori=II
                            end if
                            
                            if (abs(datain(jc+ff,hori)) < FILL .AND. mask(jc+ff,hori)==1 .AND. foundWeight(1,ff+2) .EQ. 0) THEN
                                counter = counter + 1
                                foundData(1,ff+2)=datain(jc+ff,hori)
                                foundWeight(1,ff+2)=kc
                                !print*,'found value hor',datain(jc+ff, hori)
                            end if
                       
                         ! Vertically find values around problem value
                      
                            vert=jc+kc
                            
                            if (vert .LE. 1 ) then
                                vert=1
                            else if (vert .GE. JJ) then
                                vert=JJ
                            end if
                            
                            if (abs(datain(vert,ic+ff)) < FILL .AND. mask(vert,ic+ff)==1 .AND. foundWeight(2,ff+2) .EQ. 0) THEN
                                counter = counter + 1
                                foundData(2,ff+2)=datain(vert,ic+ff)
                                foundWeight(2,ff+2)=kc
                                !print*,'found value vert',datain(vert,ic+ff)
                                !print*,'weight',kc, jc, ic
                            end if
                        end do
                        end do   
                       
        
                        ! Assign the final value to data array using weights
                        if (counter > 0) then
                            do ff=1,1,3
                            do kc 1,1,2
                                
                                dataout(jc,ic) = point
                        else if (counter==0) then
                            print*,'Failed to fill value at ', jc, ic
                            !dataout(jc,ic)=0.0
                        end if
                       ! print*,'New value at point', point, dataout(jc,ic)
                       ! print*,'New point consist of ', counter,' values close to masked point'
                       ! print*,''
                        
                   
                    end if
                end do
                ic=1
            end do
     
            !print*,'Cleaned up',total,'grid points that had bad values'
        end subroutine sweep
        
    end module cleanArray
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            