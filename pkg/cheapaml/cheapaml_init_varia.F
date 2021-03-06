#include "CHEAPAML_OPTIONS.h"

CBOP
C     !ROUTINE: CHEAPAML_INIT_VARIA
C     !INTERFACE:
      SUBROUTINE CHEAPAML_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE CHEAPAML_INIT_VARIA
C     | o Set cheapaml initial temp field
C     *==========================================================*

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "CHEAPAML.h"

C     !INPUT PARAMETERS:
C     myThid :: my Thread Id number
      INTEGER myThid
CEOP

C     !FUNCTIONS
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C     !LOCAL VARIABLES:
C     bi,bj  :: tile indices
C     i,j    :: grid-point indices
C     msgBuf :: Informational/error message buffer
      INTEGER bi, bj
      INTEGER i, j
      INTEGER iG,jG
      _RL local, localt
      _RL ssqa
      _RL recipNym1
      INTEGER iL, ioUnit
      CHARACTER*(MAX_LEN_MBUF) msgBuf
C     INTEGER prec
C     CHARACTER*(MAX_LEN_FNAM) fn

      ioUnit = standardMessageUnit
      recipNym1 = Ny - 1
      IF ( Ny.GT.1 ) recipNym1 = 1. _d 0 / recipNym1

C--   Initialise CheapAML variables in common block:
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          Tr            (i,j,bi,bj) = 0. _d 0
          qr            (i,j,bi,bj) = 0. _d 0
          Tair          (i,j,bi,bj) = 0. _d 0
          gTairm        (i,j,bi,bj) = 0. _d 0
          qair          (i,j,bi,bj) = 0. _d 0
          gqairm        (i,j,bi,bj) = 0. _d 0
          uWind         (i,j,bi,bj) = 0. _d 0
          vWind         (i,j,bi,bj) = 0. _d 0
          wWind         (i,j,bi,bj) = 0. _d 0
          solar         (i,j,bi,bj) = 0. _d 0
          ustress       (i,j,bi,bj) = 0. _d 0
          vstress       (i,j,bi,bj) = 0. _d 0
          wavesh        (i,j,bi,bj) = 0. _d 0
          wavesp        (i,j,bi,bj) = 0. _d 0
          cheapPrecip   (i,j,bi,bj) = 0. _d 0
          CheapHgrid    (i,j,bi,bj) = 0. _d 0
c         cheapPrGrid   (i,j,bi,bj) = 0. _d 0
          Cheapclouds   (i,j,bi,bj) = 0. _d 0
          Cheapdlongwave(i,j,bi,bj) = 0. _d 0
          Cheaptracer   (i,j,bi,bj) = 0. _d 0
          CheaptracerR  (i,j,bi,bj) = 0. _d 0
          gCheaptracerm (i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      IF ( startTime.EQ.baseTime .AND. nIter0.EQ.0
     &                           .AND. pickupSuff.EQ.' ' ) THEN

       IF ( AirTempFile .NE. ' ' ) THEN
         iL = ILNBLNK(AirTempFile)
         WRITE(msgBuf,'(4A)') 'CHEAPAML_INIT_VARIA: ',
     &      'Tair initialized from ->', AirTempFile(1:iL), '<-'
         CALL PRINT_MESSAGE(msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
         CALL READ_FLD_XY_RL( AirTempFile,' ',Tair,0,myThid )
       ELSE
         WRITE(msgBuf,'(4A)') 'CHEAPAML_INIT_VARIA: ',
     &      'Tair initialized using standard profile'
         CALL PRINT_MESSAGE(msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1,sNy
            DO i=1,sNx
             jG = myYGlobalLo-1+(bj-1)*sNy+j
             iG = myXGlobalLo-1+(bi-1)*sNx+i
             localt = 25. _d 0 - (jG-1)*recipNym1*10. _d 0
             localt = 20. _d 0
     &         + 10. _d 0*EXP( -( (jG-30)**2+(iG-30)**2 )/100. _d 0 )
             Tair(i,j,bi,bj) = localt
            ENDDO
           ENDDO
          ENDDO
         ENDDO
       ENDIF
       _EXCH_XY_RL( Tair, myThid )

C do specific humidity
       IF ( AirQFile .NE. ' ') THEN
         iL = ILNBLNK(AirQFile)
         WRITE(msgBuf,'(4A)') 'CHEAPAML_INIT_VARIA: ',
     &      'Qair initialized from ->', AirQFile(1:iL), '<-'
         CALL PRINT_MESSAGE(msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
         CALL READ_FLD_XY_RL( AirQFile,' ',qair,0,myThid )
       ELSE
C     default to 80% relative humidity
         WRITE(msgBuf,'(4A)') 'CHEAPAML_INIT_VARIA: ',
     &      'Qair initialized using standard profile'
         CALL PRINT_MESSAGE(msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1,sNy
            DO i=1,sNx
             local= Tair(i,j,bi,bj)+celsius2K
             ssqa = ssq0*EXP( lath*(ssq1-ssq2/local) ) / p0
             qair(i,j,bi,bj)=0.8 _d 0*ssqa
            ENDDO
           ENDDO
          ENDDO
         ENDDO
       ENDIF
       _EXCH_XY_RL( qair, myThid )

C do passive tracer
       IF ( TracerFile .NE. ' ') THEN
         iL = ILNBLNK(TracerFile)
         WRITE(msgBuf,'(4A)') 'CHEAPAML_INIT_VARIA: ',
     &      'Tracer initialized from ->', TracerFile(1:iL), '<-'
         CALL PRINT_MESSAGE(msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
         CALL READ_FLD_XY_RL( TracerFile,' ',Cheaptracer,0,myThid )
       ELSE
C default value at 290 (!)
         WRITE(msgBuf,'(4A)') 'CHEAPAML_INIT_VARIA: ',
     &      'Tracer initialized using standard profile'
         CALL PRINT_MESSAGE(msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1,sNy
            DO i=1,sNx
             Cheaptracer(i,j,bi,bj)=290.0 _d 0
            ENDDO
           ENDDO
          ENDDO
         ENDDO
       ENDIF
       _EXCH_XY_RL( Cheaptracer, myThid )

      ELSE
C Restart from cheapaml_pickups
       CALL CHEAPAML_READ_PICKUP( nIter0, myThid )
C End start-from-iter-zero if/else block
      ENDIF

C     construct cheaplayer thickness
       IF ( cheap_hFile .NE. ' ') THEN
         iL = ILNBLNK(cheap_hFile)
         WRITE(msgBuf,'(4A)') 'CHEAPAML_INIT_VARIA: ',
     &      'BL thickness taken from ->', cheap_hFile(1:iL), '<-'
         CALL PRINT_MESSAGE(msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
         CALL READ_FLD_XY_RL( cheap_hFile,' ',cheapHgrid,0,myThid )
       ELSE
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
               cheapHgrid(i,j,bi,bj) = cheapaml_h
             ENDDO
            ENDDO
          ENDDO
         ENDDO
       ENDIF
       _EXCH_XY_RL( CheapHgrid, myThid )

c!BD       IF ( cheap_prFile .NE. ' ') THEN
c!BD        write(*,*)'Conv precip taken from  ->',cheap_prFile
c!BD        CALL READ_FLD_XY_RL( cheap_prFile,' ',cheapPrGrid,0,myThid )
c!BD       ELSE
c!BD         DO bj = myByLo(myThid), myByHi(myThid)
c!BD           DO bi = myBxLo(myThid), myBxHi(myThid)
c!BD             DO j=1-OLy,sNy+OLy
c!BD               DO i=1-OLx,sNx+OLx
c!BD               cheapPrGrid(i,j,bi,bj) = 0.0 _d 0
c!BD               ENDDO
c!BD             ENDDO
c!BD           ENDDO
c!BD         ENDDO
c!BD       ENDIF
c!BD        _EXCH_XY_RL( cheapPrGrid, myThid )

C fill in outer edges
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=1-OLy,sNy+OLy
           jG = myYGlobalLo-1+(bj-1)*sNy+j
           DO i=1-OLx,sNx+OLx
             iG=myXGlobalLo-1+(bi-1)*sNx+i
             IF ( .NOT.cheapamlXperiodic .AND. iG.LT.1 ) THEN
                 Tair(i,j,bi,bj)=Tair(1,j,bi,bj)
                 qair(i,j,bi,bj)=qair(1,j,bi,bj)
                 Cheaptracer(i,j,bi,bj)=Cheaptracer(1,j,bi,bj)
                 CheapHgrid(i,j,bi,bj)=CheapHgrid(1,j,bi,bj)
             ELSEIF ( .NOT.cheapamlXperiodic .AND. iG.GT.Nx ) THEN
                 Tair(i,j,bi,bj)=Tair(sNx,j,bi,bj)
                 qair(i,j,bi,bj)=qair(sNx,j,bi,bj)
                 Cheaptracer(i,j,bi,bj)=Cheaptracer(sNx,j,bi,bj)
                 CheapHgrid(i,j,bi,bj)=CheapHgrid(sNx,j,bi,bj)
             ELSEIF ( .NOT.cheapamlYperiodic .AND.  jG.LT.1 ) THEN
                 Tair(i,j,bi,bj)=Tair(i,1,bi,bj)
                 qair(i,j,bi,bj)=qair(i,1,bi,bj)
                 Cheaptracer(i,j,bi,bj)=Cheaptracer(i,1,bi,bj)
                 CheapHgrid(i,j,bi,bj)=CheapHgrid(i,1,bi,bj)
             ELSEIF ( .NOT.cheapamlYperiodic .AND.  jG.GT.Ny ) THEN
                 Tair(i,j,bi,bj)=Tair(i,sNy,bi,bj)
                 qair(i,j,bi,bj)=qair(i,sNy,bi,bj)
                 Cheaptracer(i,j,bi,bj)=Cheaptracer(i,sNy,bi,bj)
                 CheapHgrid(i,j,bi,bj)=CheapHgrid(i,sNy,bi,bj)
             ENDIF
           ENDDO
          ENDDO
        ENDDO
       ENDDO

      IF ( debugLevel.GE.debLevB .AND. nIter0.EQ.0 ) THEN
        CALL WRITE_FLD_XY_RL('CheapHgrid', ' ', CheapHgrid, 0, myThid )
      ENDIF

      RETURN
      END
