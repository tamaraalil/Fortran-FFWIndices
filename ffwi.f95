!     PROGRAM: ffwi.f95
!     AUTHOR: Tamara Alilovic
!     FORTRAN IV PROGRAM TO CALCULATE CANADIAN FOREST
!     FIRE WEATHER INDEX FOR A DEC PDP 11 AT P.F.E.S.
!     READS DATA AND PRINTS OUT IN METRIC UNITS.

program main
      use FFWIndices
      implicit none

!     Variable declarations
      real :: prevFFMC, prevDMC, prevDC, currTemp, humidity, wind, rain
      real :: initRain, currRain, startMC, FFM, temp
      real :: effectiveRain, DMCfinal, DCFinal
      real :: spreadIndex, buildupIndex, FWI
      real, dimension(12) :: ELength, FLength
      integer, dimension(12) :: monthLength
      integer :: J, I, initHumidity, initWind, IFFM, IDMC, intDC, intSpreadIndex, intBuildupIndex, intFWI
      integer :: IDays, L, nDays, month, numDays
      character(len=30) :: inputFilename, outputFilename


!     Prompt user for input and output file names
      call filenames (inputFilename, outputFilename)
      open(9, file = inputFilename, status = 'old', action = 'read')
      open(1, file = outputFilename, status = 'unknown', action = 'write')


      write(1,1004)
1004  format(2X,'PROGRAM NO.: F-4')
100   format(I2,F4.1,F4.1)
101   format(F4.1,I4,I4,F4.1)
102   format(F4.1,F4.1,F5.1,I2,I2)


!     Reads lengths of months, and day length factors
      do J=1,12
            read(9,100) monthLength(J), ELength(J), FLength(J)
      end do


!     Reads initial values of FFMC, DMC, DC, starting month and number
!     of days of data starting month.
      read(9,102) prevFFMC, prevDMC, prevDC, month, nDays
      do J = month, 12
            numDays = monthLength(J)
      
1002        format(10(/),1X,'  DATE  TEMP  RH   WIND  RAIN   FFMC   DMC   DC     ISI   BUI   FWI'/)
1003        format(1X,'  mm/dd   C   %   km/hr   mm   '/)
            if (J == month) then
                  IDays = monthLength(J) - nDays + 1
            else 
                  IDays = 1
            end if

!           Reads daily weather data
            L = 0
            do I = IDays, numDays
                  L = L + 1
                  read(9,101,END = 2000) temp, initHumidity, initWind, initRain
                  if (L == 1) then
                        write(1,1002)
                        write(1,1003)
                  end if
                  currTemp = temp
                  humidity = initHumidity
                  wind = initWind
                  rain = initRain


!                 call subroutines to perform calcuations
                  call fineFuelMoisture(initRain, prevFFMC, currRain, startMC, FFM, wind, humidity, temp)

                  call duffMoisture(temp, prevDMC, currRain, humidity, J, ELength, DMCfinal, effectiveRain, initRain)

                  call drought(temp, J, currRain, effectiveRain, prevDC, initRain, FLength, DCFinal)

                  call spreadBuildupFWI(FFM, spreadIndex, buildupIndex, DCFinal, FWI, wind, DMCfinal)

!                 Convert calculations to integers
                  intDC = int(DCFinal + 0.5)
                  IFFM = int(FFM + 0.5)
                  IDMC = int(DMCfinal + 0.5)
                  intSpreadIndex = int(spreadIndex + 0.5)
                  intBuildupIndex = int(buildupIndex + 0.5)
                  intFWI = int(FWI + 0.5)

!                 Print results
                  write(1,1001) J,I,currTemp,initHumidity,initWind,rain,IFFM,IDMC,intDC,intSpreadIndex,intBuildupIndex,intFWI
1001              format(1X,2I3,F6.1,I4,I6,F7.1,6I6)
                  prevFFMC = FFM
                  prevDMC = DMCfinal
                  prevDC = DCFinal
                  end do
            end do
2000  stop
end program main