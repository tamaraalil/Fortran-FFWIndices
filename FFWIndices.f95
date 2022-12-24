!
!   PROGRAM: FFWIncidices.f95
!   AUTHOR: Tamara Alilovic
!   SUBROUTINES FOR ffwi.f95
!

module FFWIndices
    implicit none
contains

    ! Get input and output file names
    subroutine filenames (inputFilename, outputFilename)
        character(len=30), intent(inout) :: inputFilename, outputFilename
        write (*,*) 'Enter input file name.'
        read(*,*) inputFilename
        write (*,*) 'Enter output file name.'
        read(*,*) outputFilename
    end subroutine filenames

    ! Fine Fuel Moisture code
    subroutine fineFuelMoisture (initRain, prevFFMC, currRain, startMC, FFM, wind, humidity, temp)
        real :: finalMC, Z, X, FFMC, FFMCAfterRain, correctionTerm, dryFFEMC, wetFFEMC
        real, intent(inout) :: initRain, prevFFMC, currRain, startMC, FFM, wind, humidity, temp
        if (initRain <= 0.5) then
            initRain = 0.0
            FFMCAfterRain = prevFFMC
        end if
        if (initRain > 0.5) then
                currRain = initRain
                if (currRain <= 1.45) then
                    FFMC = 123.85-(55.6*ALOG(currRain+1.016))
                else if ((currRain-5.75) <= 0) then
                    FFMC = 57.87-(18.2*ALOG(currRain-1.016))
                else if ((currRain-5.75) > 0) then
                    FFMC = 40.69-(8.25*ALOG(currRain-1.905))
                end if
                correctionTerm = 8.73*EXP(-0.1117*prevFFMC)
                FFMCAfterRain = (prevFFMC/100.)*FFMC+(1.0-correctionTerm)
                if (FFMCAfterRain < 0) then
                    FFMCAfterRain = 0.0
                end if
        end if 
        startMC = 101.-FFMCAfterRain
        dryFFEMC = 0.942*(humidity**0.679)+(11.*EXP((humidity-100.)/10.))+0.18*(21.1-temp)*(1.-1./EXP(0.115*humidity))
        if ((startMC - dryFFEMC) < 0) then
                wetFFEMC = 0.618*(humidity**0.753)+(10.*EXP((humidity-100.)/10.))+0.18*(21.1-temp)*(1.-1./EXP(0.115*humidity))
                if (startMC < wetFFEMC) then
                    finalMC = wetFFEMC-(wetFFEMC-startMC)/1.9953
                else 
                    finalMC = startMC
                end if
        end if
        if ((startMC - dryFFEMC) == 0) then
                finalMC = startMC 
        end if
        if ((startMC - dryFFEMC) > 0) then
                Z = 0.424*(1.-(humidity/100.)**1.7)+(0.0694*(wind**0.5))*(1.-(humidity/100.)**8)
                X = Z*(0.463*(EXP(0.0365*temp)))
                finalMC = dryFFEMC+(startMC-dryFFEMC)/10.**X
        end if
        FFM=101.-finalMC
        if (FFM > 101) then
                FFM=101.
        end if
        if (FFM < 0) then
                FFM = 0.0
        end if
    end subroutine fineFuelMoisture

    ! DUFF Moisture code
    subroutine duffMoisture(temp, prevDMC, currRain, humidity, J, ELength, DMCfinal, effectiveRain, initRain)
        real :: rainDMC, moistureContent, WMR, drying, slope
        real, intent(inout) :: temp, currRain, humidity, DMCfinal, effectiveRain, initRain, prevDMC
        integer, intent(in) :: J
        real, dimension(12), intent(inout) :: ELength
        if ((temp + 1.1) >= 0) then
            drying=1.894*(temp+1.1)*(100.-humidity)*(ELength(J)*0.0001)
        else
                temp = -1.1
                drying=1.894*(temp+1.1)*(100.-humidity)*(ELength(J)*0.0001)
        end if
        drying=1.894*(temp+1.1)*(100.-humidity)*(ELength(J)*0.0001)
        if (initRain <= 1.5) then
                rainDMC = prevDMC
        end if
        if (initRain > 1.5) then
                currRain = initRain
                effectiveRain = 0.92*currRain-1.27
                moistureContent = 20.0+280./EXP(0.023*prevDMC)
                if (prevDMC <= 33) then
                    slope = 100./(0.5+0.3*prevDMC)
                else 
                    if ((prevDMC - 65) <= 0) then
                            slope = 14.-1.3*ALOG(prevDMC)
                    end if 
                    if ((prevDMC - 65) > 0) then
                            slope = 6.2*ALOG(prevDMC)-17.2
                    end if 
                end if 
                WMR = moistureContent+(1000.*effectiveRain)/(48.77+slope*effectiveRain)
                rainDMC = 43.43*(5.6348-ALOG(WMR-20.))
        end if 
        if (rainDMC < 0) then
                rainDMC = 0.0
        end if
        DMCfinal = rainDMC + drying
    end subroutine duffMoisture

    ! Drought code
    subroutine drought(temp, J, currRain, effectiveRain, prevDC, initRain, FLength, DCFinal)
        real, intent(inout) :: temp, currRain, effectiveRain, prevDC, initRain, DCFinal
        real :: rainDC, SMI, PE
        integer, intent(in) :: J
        real, dimension(12), intent(inout) :: FLength
        if ((temp + 2.8) < 0) then
            temp = -2.8
        end if
        PE=(.36*(temp+2.8)+FLength(J))/2.
        if (initRain > 2.8) then
                currRain = initRain
                effectiveRain = 0.83*currRain-1.27
                SMI = 800.*EXP(-prevDC/400.)
                rainDC = prevDC-400.*ALOG(1.+((3.937*effectiveRain)/SMI))
                if (rainDC <= 0) then
                    rainDC = 0.0
                    DCFinal = rainDC + PE
                else
                    DCFinal = rainDC + PE
                end if
        else
                rainDC = prevDC
                DCFinal = rainDC+PE
        end if
        if (DCFinal < 0) then
                DCFinal = 0.0
        end if
    end subroutine drought

    ! Initial spread index, buildup index, fire weather index
    subroutine spreadBuildupFWI(FFM, spreadIndex, buildupIndex, DCFinal, FWI, wind, DMCfinal)
        real, intent(inout) :: FFM, spreadIndex, buildupIndex, DCFinal, FWI, wind, DMCfinal
        real :: BB, CC, todayFFMC, SF, DMC, logFWI
        todayFFMC = 101.-FFM
        SF = 19.1152*EXP(-0.1386*todayFFMC)*(1.+todayFFMC**4.65/7950000.)
        spreadIndex = SF*EXP(0.05039*wind)
        buildupIndex = (0.8*DCFinal*DMCfinal)/(DMCfinal+0.4*DCFinal)
        if (buildupIndex < DMCfinal) then
            DMC = (DMCfinal-buildupIndex)/DMCfinal
            CC = 0.92+(0.0114*DMCfinal)**1.7
            buildupIndex = DMCfinal-(CC*DMC)
            if (buildupIndex < 0) then
                    buildupIndex = 0
            end if
        end if
        if (buildupIndex <= 80) then
            BB = 0.1*spreadIndex*(0.626*buildupIndex**0.809+2.)
        else 
            BB = 0.1*spreadIndex*(1000./(25.+108.64/EXP(0.023*buildupIndex)))
        end if
        if ((BB - 1.0) > 0) then
            logFWI=2.72*(0.43*ALOG(BB))**0.647
            FWI=EXP(logFWI)
        else
            FWI=BB 
        end if
    end subroutine spreadBuildupFWI

end module FFWIndices