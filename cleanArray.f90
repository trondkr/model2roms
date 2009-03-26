        
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
          
           
            integer II, JJ, KK, ic, jc, jcc, icc, fill, kc, total
            integer go1,go2,go3,go4,go5,go6,go7,go8,go9,go10,go11,go12,go13,go14,go15,go16
            double precision, dimension(JJ,II) :: dataout,datain
            double precision, dimension(JJ,II) :: mask
            double precision counter, point, neg, pos
            
!f2py intent(in,out) dataout
!f2py intent(in) datain, mask
!f2py intent(in) II, JJ, FILL, KK
!f2py intent(hide) ic, jc, icc, jcc, counter, point, total, neg, pos

            !FILL=10000
            !KK=20
            !print*,'---> Started cleaning of array'
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
                        
                        go1=0
                        go2=0
                        go3=0
                        go4=0
                        go5=0
                        go6=0
                        go7=0
                        go8=0
                        go9=0
                        go10=0
                        go11=0
                        go12=0
                        go13=0
                        go14=0
                        go15=0
                        go16=0
                        
                        ! ROW 1 - keep JC, JC+1, JC-1 constant, move I +- kc positions
                        do kc=1,KK
                            neg=ic-kc
                            pos=ic+kc
                            
                            if (neg .LE. 1) then
                                neg=1
                            end if
                            if (pos .GE. II) then
                                pos=II
                            end if
                            
                            if (abs(datain(jc+1,neg)) < FILL .AND. go1==0) THEN
                                point = point + datain(jc+1,neg)
                                counter = counter + 1
                                go1=1
                            end if
                            if (abs(datain(jc+1,ic)) < FILL .AND. go2==0) THEN
                                point = point + datain(jc+1,ic)
                                counter = counter + 1
                                go2=1
                            end if
                            if (abs(datain(jc+1,pos)) < FILL .AND. go3==0) THEN
                                point = point + datain(jc+1,pos)
                                counter = counter + 1
                                go3=1
                            end if

                            ! ROW 2
                            if (abs(datain(jc,ic-kc)) < FILL .AND. go4==0) THEN
                                point = point + datain(jc,neg)
                                counter = counter + 1
                                go4=1
                            end if
                            if (abs(datain(jc,ic+kc)) < FILL .AND. go5==0) THEN
                                point = point + datain(jc,pos)
                                counter = counter + 1
                                go5=1
                            end if
                        
                            ! ROW 3
                            if (abs(datain(jc-1,ic-1)) < FILL .AND. go6==0) THEN
                                point = point + datain(jc-1,neg)
                                counter = counter + 1
                                go6=1
                            end if
                            if (abs(datain(jc-1,ic)) < FILL .AND. go7==0) THEN
                                point = point + datain(jc-1,ic)
                                counter = counter + 1
                                go7=1
                            end if
                            if (abs(datain(jc-1,ic+1)) < FILL .AND. go8==0) THEN
                                point = point + datain(jc-1,pos)
                                counter = counter + 1
                                go8=1
                            end if
                        
                            ! SMOOTH VERTICALLY
                            ! ROW 1 - keep IC, IC+1, IC-1 constant, move J +- kc positions
                     
                            neg=jc-kc
                            pos=jc+kc
                            
                            if (neg .LE. 1) then
                                neg=1
                            end if
                            if (pos .GE. JJ) then
                                pos=JJ
                            end if
                            
                            if (abs(datain(pos,ic-1)) < FILL .AND. go9==0) THEN
                                point = point + datain(pos,ic-1)
                                counter = counter + 1
                                go9=1
                            end if
                            if (abs(datain(jc,ic-1)) < FILL .AND. go10==0) THEN
                                point = point + datain(jc,ic-1)
                                counter = counter + 1
                                go10=1
                            end if
                            if (abs(datain(neg,ic-1)) < FILL .AND. go11==0) THEN
                                point = point + datain(neg,ic-1)
                                counter = counter + 1
                                go11=1
                            end if
                        
                            ! ROW 2 (midpoints)
                            if (abs(datain(pos,ic)) < FILL .AND. go12==0) THEN
                                point = point + datain(pos,ic)
                                counter = counter + 1
                                go12=1
                            end if
                            if (abs(datain(neg,ic)) < FILL .AND. go13==0) THEN
                                point = point + datain(neg,ic)
                                counter = counter + 1
                                go13=1
                            end if
                        
                            ! ROW 3
                            if (abs(datain(pos,ic+1)) < FILL .AND. go14==0) THEN
                                point = point + datain(pos,ic+1)
                                counter = counter + 1
                                go14=1
                            end if
                            if (abs(datain(jc,ic+1)) < FILL .AND. go15==0) THEN
                                point = point + datain(jc,ic+1)
                                counter = counter + 1
                                go15=1
                            end if
                            if (abs(datain(neg,ic+1)) < FILL .AND. go16==0) THEN
                                point = point + datain(neg,ic+1)
                                counter = counter + 1
                                go16=1
                            end if
                        
                        end do
                        ! Assign the final value to data array
                        if (counter > 0) then
                            point = point/(counter*1.0)
                            dataout(jc,ic) = point
                        else if (counter==0 .OR. point > FILL) then
                            dataout(jc,ic)=0.0
                        end if
                       ! print*,'New value at point', point, dataout(jc,ic)
                       ! print*,'New point consist of ', counter,' values close to masked point'
                       ! print*,''
                        
                        
                    end if
                end do
            end do
     
            !print*,'Cleaned up',total,'grid points that had bad values'
        end subroutine sweep
        
    end module cleanArray
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            