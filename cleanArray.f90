        
Module cleanArray

        implicit none
                
        contains
  
 
 
 subroutine sweep(bathymetry,sodadepth,maxpoints,FILL,KK,VV,datain,dataout,mask,II,JJ)
       
            ! ------------------------------------------------------------------
            ! Program : sweep
            !
            ! This routine finds points in the interpolated fields that was masked by
            ! the original land-array, but which is no longer masked in the new land-ocean array.
            ! Sometimes this leads to ocean posint in the new grid that are still
            ! defined as land because of the previous masking. This routine sweeps through
            ! the array and replaces fill values over ocean with the nearest points values.
            !
            ! The routine finds the points in KKxVV distance from missing value points.
            ! This means that a maximum of KK points in the horizontal direction and a maximum
            ! of VV points in the vertical are searched for values to be used to fill in
            ! the missing value. A maximum of "maxpoints" within the distance KK and VV are
            ! then used to make sure you use mostly values close to missing value. However, with the
            ! ability to have high values of KK and VV you make sure that a value is
            ! found to fill in if no-values exists nearby. The values are weighted
            ! with the distance from the missing value point. If no value is found within KKxVV distance,
            ! the search radius is increased (INCREASED=1) and for each iteration no values are found,
            ! the KK and VV distance increases with 10 points. As soon as a good value is found,
            ! the search radius decreases back to the original size determined in grd.py. This
            ! increase/decrease is necessary for some points in difficult topographical areas (like Bay of Fundy).
            ! The program increases the searh radius 3 times before giving up to find values.
            !
            ! Trond Kristiansen, March and April 2009
            ! Rutgers University, NJ.
            !
            ! Compile with:
            ! f2py-64 --verbose -DF2PY_REPORT_ON_ARRAY_COPY=1 --fcompiler=intelem -c -m cl cleanArray.f90
            ! or
            ! f2py --verbose --fcompiler=intelem -c -m cl cleanArray.f90
            !
            ! USAGE:
            ! import cl
            ! Zin = cl.cleanArray.sweep(Zin,Zout, maskarray,i,j,k)
            !
            !                 (I,J+1)                                                               
            !   |----------------O----------------|                                                                        
            !   |                |                |                                                
            !   |                |                |                                                  
            !  (I-1,J)        X (I,J)          (I+1,J)
            !   |                |                |   
            !   |                |                |                                                  
            !   |----------------O----------------|                                                  
            !                 (I,J-1)                                           
            !
            ! From each of these 6 points (5 unique) we move up and down VV points,
            ! and KK points horisontally. A maximum of "maxpoints" are stored in each direction (foundData array)
            ! together with the distance from the missing value point (foundWeight array) 
    
            ! ------------------------------------------------------------------
          
           
            integer maxpoints, II, JJ, KK, VV, ic, jc, jcc, icc, fill, kc, ff, total
            integer pointsHpos, pointsHneg, pointsVpos, pointsVneg, points, INCREASED, orgVV, orgKK, TRIES

            double precision, dimension(JJ,II) :: dataout,datain
            double precision, dimension(JJ,II) :: mask
            double precision, dimension(4,3,maxpoints)  :: foundWeight, foundData
             
            double precision counter, hori, vert, sumweights, distance, sodadepth, weight
            double precision, dimension(JJ,II) :: bathymetry
            
!f2py intent(in,out,overwrite) dataout
!f2py intent(in,overwrite) datain, mask, bathymetry
!f2py intent(in) maxpoints, II, JJ, FILL, KK, VV, sodadepth
!f2py intent(hide) ic, jc,kc, ff, icc, jcc, counter, total, hori, verti
!f2py intent(hide) foundWeight, foundData, sumweights, distance, weight
!f2py intent(hide) pointsHpos, pointsHneg, pointsVpos, pointsVneg, points, INCREASED
!f2py intent(hide) orgVV, orgKK, TRIES
            
           ! print*,'---> Started cleaning of array with box size of',KK,VV
            total   = 0    
            orgVV=VV
            orgKK=KK
            
            do jc=1,JJ
                do ic=1,II
                    TRIES=1
                    counter  = 0
                    pointsHpos  = 1
                    pointsVpos  = 1
                    pointsHneg  = 1
                    pointsVneg  = 1
                    
                    if (abs(datain(jc,ic)) > FILL .AND. mask(jc,ic)==1 .AND. sodadepth >= -(bathymetry(jc,ic))) THEN
                        ! Have to check only grid cells that are deeper than the sodadepth. If not, we will fill in all
                        ! values toward the coast for each iteration (sodadepth >= -(bathymetry(jc,ic)))).
                        total=total+1
                        foundWeight(:,:,:)=0
                        foundData(:,:,:)=0
                        
50                      continue

                        do ff=-1,1,1
                        do kc=1,KK
                        
                            hori=ic+kc
                            
                            if (hori .LE. 1 ) then
                                hori=1
                            else if (hori .GE. II) then
                                hori=II
                            end if
                         
                            vert=jc+ff
                            
                            if (vert .LE. 1 ) then
                                vert=1
                            else if (vert .GE. JJ) then
                                vert=JJ
                            end if

                            if (abs(datain(vert,hori)) < FILL .AND. mask(vert,hori)==1 .AND. &
                            & foundWeight(1,ff+2,pointsHneg) .EQ. 0 .AND. pointsHneg .LE. maxpoints) THEN
                                counter = counter + 1
                                foundData(1,ff+2,pointsHneg)=datain(vert,hori)
                                foundWeight(1,ff+2,pointsHneg)=sqrt(((hori) - ic)**2.0 + ((vert) - jc)**2.0)
                                pointsHneg=pointsHneg+1
                            end if
                            
                            hori=ic+kc*(-1)
                            if (hori .LE. 1 ) then
                                hori=1
                            else if (hori .GE. II) then
                                hori=II
                            end if
                            
                            if (abs(datain(vert,hori)) < FILL .AND. mask(vert,hori)==1 .AND. &
                            & foundWeight(3,ff+2,pointsHpos) .EQ. 0 .AND. pointsHpos .LE. maxpoints) THEN
                                counter = counter + 1
                                foundData(3,ff+2,pointsHpos)=datain(vert,hori)
                                foundWeight(3,ff+2,pointsHpos)=sqrt(((hori) - ic)**2.0 + ((vert) - jc)**2.0)
                                pointsHpos=pointsHpos+1
                            end if
                        end do
                        end do
                        
                        do ff=-1,1,1
                        do kc=1,VV
                             ! Vertically find values around problem value
                      
                            vert=jc+kc
                            
                            if (vert .LE. 1 ) then
                                vert=1
                            else if (vert .GE. JJ) then
                                vert=JJ
                            end if
                            
                            hori=ic+ff
                            
                             if (hori .LE. 1 ) then
                                hori=1
                            else if (hori .GE. II) then
                                hori=II
                            end if
                            
                            if (abs(datain(vert,hori)) < FILL .AND. mask(vert,hori)==1 &
                            & .AND. foundWeight(2,ff+2, pointsVneg) .EQ. 0 .AND. pointsVneg .LE. maxpoints) THEN
                                counter = counter + 1
                                foundData(2,ff+2, pointsVneg)=datain(vert,hori)
                                foundWeight(2,ff+2, pointsVneg)=sqrt(((hori) - ic)**2.0 + ((vert) - jc)**2.0)
                                pointsVneg=pointsVneg+1
                            end if
                       
                       
                            vert=jc+kc*(-1)
                            
                            if (vert .LE. 1 ) then
                                vert=1
                            else if (vert .GE. JJ) then
                                vert=JJ
                            end if
                            
                           
                             if (abs(datain(vert,hori)) < FILL .AND. mask(vert,hori)==1 &
                            & .AND. foundWeight(4,ff+2, pointsVpos) .EQ. 0 .AND. pointsVpos .LE. maxpoints) THEN
                                counter = counter + 1
                                foundData(4,ff+2, pointsVpos)=datain(vert,hori)
                                foundWeight(4,ff+2, pointsVpos)=sqrt(((ic+ff) - ic)**2.0 + ((vert) - jc)**2.0)
                                pointsVpos=pointsVpos+1
                            end if
                        end do
                        end do   
                      
                        ! Assign the final value to data array using weights
                        if (counter > 0) then
                            
                            distance = sum(abs(foundWeight))
                            sumweights=0.0
                            !print*,'total number of points',counter
                            
                            dataout(jc,ic)=0.0
                            do kc=1,4
                            do ff=1,3
                            do points=1,maxpoints
                                if (foundWeight(kc,ff,points) .NE. 0) then
                                    
                                    if (counter > 1) then
                                        weight=((1.0-(abs(foundWeight(kc,ff,points)))/distance)/(counter-1.0)*1.0)
                                        sumweights=sumweights + ((1.0-(abs(foundWeight(kc,ff,points)))/distance)/(counter-1.0)*1.0)
                                    else
                                        weight=1.0
                                        sumweights=sumweights + 1.0
                                    end if
                                    dataout(jc,ic) = dataout(jc,ic) + foundData(kc,ff,points)*weight
                                    
                                    !print*,foundData(kc,ff,points), distance, abs(foundWeight(kc,ff,points))/distance
                                    !print*,'total',sumweights, distance, 1-(abs(foundWeight(kc,ff,points))
                                end if
                            end do
                            end do
                            end do
                        end if

                        if (counter==0) then
                            !print*,'no data', dataout(jc,ic), mask(jc,ic), jc,ic
                            ! Here we increase the search area as no points within KKxVV
                            ! both vertically and horisontally was found.
                            KK=KK+10
                            VV=VV+10
                            INCREASED=1
                            TRIES=TRIES+1
                            !print*,'Increased distances for values', TRIES, KK, VV
                            if (TRIES .LT. 3) then
                                goto 50
                            else
                                goto 100
                            end if
                        else if (counter >0 .AND. INCREASED .EQ. 1) then
                            ! Good values found, now reset the search area to
                            ! what is defined in grd.py
                            INCREASED=0
                            KK=orgKK
                            VV=orgVV
                            !print*,'reduced distances', jc, ic, KK, VV
                        end if
                       
                    end if
                    
100                 continue                 
                end do
            end do
           
            !print*,'Number of points filled', total
            ! ----
           
        end subroutine sweep
        
    end module cleanArray
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            