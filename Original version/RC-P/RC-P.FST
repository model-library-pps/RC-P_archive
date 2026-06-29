*-------------------------------------------------------------------------*
* Copyright 2013. Wageningen University, Plant Production Systems group,  *
* P.O. Box 430, 6700 AK Wageningen, The Netherlands.                      *
* You may not use this work except in compliance with the Licence.        *
* You may obtain a copy of the Licence at:                                *
*                                                                         *
* http://models.pps.wur.nl/content/licence-agreement                      *
*                                                                         *
* Unless required by applicable law or agreed to in writing, software     *
* distributed under the Licence is distributed on an "AS IS" basis,       *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
*-------------------------------------------------------------------------*

DEFINE_CALL TOTRC (INPUT,INPUT, INPUT, OUTPUT, OUTPUT, OUTPUT)

TITLE RC

* FST version of RC-for model to calculate longterm recovery of applied
* fertilizer phosphorus
*
************************************************************************
***   1. Initial conditions and run control
************************************************************************

INITIAL


*     INITIAL AMOUNT OF LABILE AND OF STABLE PHOSPHORUS (KG P HA-1)
      LPI=MINP/FRL1
      SPI= LPI*TC2/TC1 + MINP*TC2

*     AMOUNT OF LABILE AND OF STABLE PHOSPHORUS, UNFERTILIZED AND
*     FERTILIZED (KG P HA-1)
      LPUI= LPI
      LPFI= LPI
      SPUI= SPI
      SPFI= SPI

INCON  ZERO= 0.

* Time step is one year
TIMER STTIME = 0.; FINTIM = 12.; DELT = 1.; PRDEL = 1.; IPFORM = 5
TRANSLATION_GENERAL DRIVER='EUDRIV'
PRINT RUNNR, FERT,UPTLPF,LPF,SPF,UPTLPU,LPU,SPU,RECPU,RECPF,DIF

DYNAMIC


************************************************************************
***   2. State variables
************************************************************************

*     AMOUNT OF LABILE AND OF STABLE PHOSPHORUS (KG P HA-1),WITH AND
*     WITHOUT FERTILIZATION
      LPF= INTGRL(LPFI, DLPF)
      SPF= INTGRL(SPFI, DSPF)

      LPU= INTGRL(LPUI, DLPU)
      SPU= INTGRL(SPUI, DSPU)

************************************************************************
***   3. Rate variables
************************************************************************

*     FERTILISER APPLICATION
      FERT = AFGEN(FERTTB, TIME)

*     PHOSPHORUS UPTAKE FROM LABILE POOL,WITH AND WITHOUT FERTILIZER
*     (KG P HA-1 YR-1)
      UPTLPU= FRL1 * AMIN1(PRT1,LPU) + ...
      FRL2 * AMIN1(PRT2-PRT1,MAX(0.,LPU-PRT1)) + ...
      FRL3 * AMIN1(PRT3-PRT2,MAX(0.,LPU-PRT2)) + ...
      FRL4 * AMIN1(PRT4-PRT3,MAX(0.,LPU-PRT3)) + ...
      FRL5 * AMIN1(PRT5-PRT4,MAX(0.,LPU-PRT4)) + ...
      FRL6 * MAX(0.,LPU-PRT5)

      UPTLPF= FRL1 * AMIN1(PRT1,LPF) + ...
      FRL2 * AMIN1(PRT2-PRT1,AMAX1(0.,LPF-PRT1)) + ...
      FRL3 * AMIN1(PRT3-PRT2,AMAX1(0.,LPF-PRT2)) + ...
      FRL4 * AMIN1(PRT4-PRT3,AMAX1(0.,LPF-PRT3)) + ...
      FRL5 * AMIN1(PRT5-PRT4,AMAX1(0.,LPF-PRT4)) + ...
      FRL6 * AMAX1(0.,LPF-PRT5)

*     TOTAL AMOUNT OF PHOSPHORUS RECOVERED WITH AND WITHOUT
*     FERTILIZATION (KG P HA-1 )
      CALL TOTRC (UPTLPF, UPTLPU, TIME, RECPF, RECPU, DIF)

*     RATE OF PRODUCTION OF LABILE AND OF STABLE PHOSPHORUS (KG P HA-1
*     YR-1), WITHOUT FERTILIZATION
      DLPU= SPU/TC2 - LPU/TC1 - UPTLPU
      DSPU= LPU/TC1 - SPU/TC2 + MINP

*     RATE OF PRODUCTION OF LABILE AND OF STABLE PHOSPHORUS (KG P HA-1
*     YR-1), WITH FERTILIZATION
      DLPF= - LPF/TC1 + SPF/TC2 - UPTLPF + FERT*FRF1
      DSPF= - SPF/TC2 + LPF/TC1 + MINP + FERT*FRF2


  PARAM RUN    = 0.
      RUNNR  = RUN


************************************************************************
***   4. Parameters
************************************************************************

*      TIME CONSTANT OF PRECIPITATION OF LABILE PHOSPHORUS (=1) AND OF
*      DISSOLUTION OF STABLE PHOSPHORUS (=2) (YR)
PARAM TC1= 5.; TC2= 30.; MINP= 3.

*     FRACTION OF APPLIED PHOSPHORUS FERTILIZER PARTITIONED TO LABILE
*     PHOSPHORUS (=1) AND TO STABLE PHOSPHORUS (=2) (KG P/KG P)
PARAM FRF2=0.2;  FRF1= 0.8

*     UPTAKE FRACTIONS OF LABILE PHOSPHORUS (KG P/KG P)
PARAM FRL1= 0.210; FRL2= 0.079; FRL3= 0.040; FRL4= 0.018; FRL5= 0.0; FRL6= 0.0

*     MAXIMUM AMOUNT OF LABILE PHOSPHORUS FOR WHICH UPTAKE
*     FRACTIONS ( FRL1-6) APPLY (KG P HA-1)
PARAM PRT1= 57.4; PRT2=113.4; PRT3=225.4; PRT4=10000.; PRT5=10000.


************************************************************************
***   5. Data
************************************************************************

*     FERTILIZER APPLICATION AT BEGINNING OF YEAR 1 TO 12 (KG P HA-1 YR-1)
FUNCTION FERTTB = 0., 35., 1., 35., 2., 35., 3., 35., 4., 35., 5., 0., 6., 0.,...
          7., 0., 8., 0., 9., 0., 10., 0., 11., 0., 12., 0.

**************************** RERUNS ************************************

* rerun
END
PARAM RUN = 1.; TC2= 20.
END
PARAM RUN= 2.,  TC2= 30., MINP= 0.5

STOP

************************** SUBROUTINES ****************************************

* -----------------------------------------------------------------------------
* Subroutine TOTRC
* Purpose: calculation of total P uptake  and the difference between fertilised and unfertilised P uptake
* ------------------------------------------------------------------------------
         SUBROUTINE TOTRC (UPTLPF, UPTLPU, TIME, RECPF, RECPU, DIF)

         IF (TIME .EQ. 0.) RECPF= 0.
         IF (TIME .EQ. 0.) RECPU= 0.
         RECPF= RECPF + UPTLPF
         RECPU= RECPU + UPTLPU

C       DIFFERENCE IN RECOVERY OF PHOSPHORUS, DUE TO FERTILISER APPLICATION
C       (KG P HA-1)
         DIF= RECPF - RECPU

         RETURN
         END
