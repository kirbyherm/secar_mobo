*******************************************************************************
*                                                                             *
*                                                                             *
*                        ELEMENTARY OPERATIONS PACKAGE                        *
*                                                                             *
*                           PART OF THE COSY SYSTEM                           *
*                                                                             *
*                                 VERSION 10.0                                *
*                                                                             *
*                           UPDATED IN MARCH, 2017                            *
*                                                                             *
*                                 INCLUDING                                   *
*                         DIFFERENTIAL ALGEBRA PACKAGE                        *
*                        UNDER LICENSE FROM MARTIN BERZ                       *
*                                                                             *
*     COPYRIGHT (C) MICHIGAN STATE UNIVERSITY BOARD OF TRUSTEES 1995 - 2017   *
*     SUBJECT TO LICENSING AGREEMENT - NOT TO BE DISTRIBUTED                  *
*                                                                             *
*     DISTRIBUTED BY M. BERZ AND K. MAKINO                                    *
*     DEPARTMENT OF PHYSICS AND ASTRONOMY                                     *
*     MICHIGAN STATE UNIVERSITY                                               *
*     EAST LANSING, MI 48824, USA                                             *
*     BERZ@MSU.EDU                                                            *
*     517-884-5583 (PHONE)                                                    *
*                                                                             *
*     VERSION FOR SETUP IN LINE THAT IS NOT COMMENTED OUT;                    *
*     TO CREATE DIFFERENT VERSIONS, USE THE PROGRAM 'VERSION'                 *
                                                                         *NORM*
*MPI                                                                     *MPI *
                                                                         *RND *
*                                                                             *
*******************************************************************************
*
*%%
      BLOCKDATA FOXBLD
*     ****************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*     THIS BLOCKDATA CONTAINS AND INITIALIZES ALL THE PARAMETERS AND
*     COMMON BLOCKS USED BY THE DIFFERENTIAL ALGEBRA PACKAGE AND
*     THE MEMORY MANAGEMENT
*
*     PARAMETERS:
*
*     LVAR: MAXIMUM NUMBER OF VARIABLES;       CAN BE CHANGED QUITE ARBITRARILY
*     LMEM: LENGTH OF MAIN STORAGE STACK;      CAN BE CHANGED QUITE ARBITRARILY
*     LDIM: MAXIMUM NUMBER OF ENTRIES OF NDIM; CAN BE CHANGED QUITE ARBITRARILY
*
*     LEA: MAXIMUM NUMBER OF MONOMIALS;        CAN BE INCREASED FOR LARGE NO,NV
*     LIA: DIMENSION OF IA1,IA2;               CAN BE INCREASED FOR LARGE NO,NV
*     LNO: MAXIMUM ORDER;                      CAN BE INCREASED TO ABOUT 1000
*     LNV: MAXIMUM NUMBER OF VARIABLES;        CAN BE INCREASED TO ABOUT 1000
*
*     ALL THE CHANGES IN THE VALUES OF PARAMETERS HAVE TO BE MADE BY GLOBAL
*     SUBSTITUTIONS IN ALL SUBROUTINES.
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
*     NTYP: INTEGER DESCRIBING THE TYPE
*     NBEG: BEGINNING ADDRESS IN CC, NC
*     NEND: ENDING ADDRESS IN CC, NC (MOMENTARY)
*     NMAX: MAXIMUM ENDING ADDRESS RESERVED
*
*     CC:   CONTAINS IMEM DOUBLE PRECISION NUMBERS
*     NC:   CONTAINS IMEM POSITION INTEGERS
*     NDIM: DIMENSIONS OF ARRAY VARIABLES
*     IDIM: MOMENTARY NUMBER OF ALLOCATED ARRAY LENGTH
*     IVAR: MOMENTARY NUMBER OF ALLOCATED VARIABLES
*     IMEM: MOMENTARY NUMBER OF ALLOCATED MEMORY POSITIONS
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
*     IE1:      CHARACTERISTIC INTEGER 1
*     IE2:      CHARACTERISTIC INTEGER 2
*     IEO:      ARRAY OF ORDERS
*     IA1:      REVERSE TO IE1 (CF DAINI)
*     IA2:      REVERSE TO IE2 (CF DAINI)
*
*     CDA:      DA/CD STORAGE FOR MULTIPLICATION ETC. (CF DAPAC/CDPAC)
*     NCFLT:    DA FILTERING TEMPLATE STORAGE (CF DAFSET)
*     LFLT:     DA FILTERING MODE
*     NFLT:     NUMBER OF TERMS IN DA FILTERING TEMPLATE
*
*     IEW:      WEIGHTED DA ORDER (CF DAINI)
*     IED:      WEIGHTED DA CODING NUMBER (CF DAINI,IJJ)
*     LEW:      FLAG OF WEIGHTED DA ORDER (CF DAINI)
*     LEWI:     FLAG OF INPUT OF WEIGHTED DA ORDER (CF DAINI,DANOTW)
*     IESP:     SPLITTER FOR DA CODING (CF DADEC,DAENC,IJJ)
*
*     NOMAX:    MAXIMUM REQUESTED ORDER  (CF DAINI)
*     NVMAX:    MAXIMUM REQUESTED NUMBER OF VARIABLES (CF DAINI)
*     NMMAX:    MAXIMUM NUMBER OF MONOMIALS FOR NOMAX, NVMAX (CF DAINI)
*     NOCUT:    MOMENTARY TRUNCATION ORDER
*     EPS:      TRUNCATION ACCURACY (CAN BE SET BY USER) (CF DAEPS)
*               BESIDES THE DATA STATEMENT INITIALIZATION, INRNDO INITIALIZES
*               IT TO EPSM*1.D-4, WHICH WOULD BE ABOUT 2.D-20. DAINI AND BASIC
*               TM ROUTINES CALL INRNDO TO DETERMINE EPS.
*     EPSMAC:   SMALLEST NUMBER REPRESENTABLE ON MACHINE (PESSIMISTIC ESTIMATE)
*
*     TMT:      TM TALLY COUNTER
*     TMS:      TM SWEEPING COUNTER
*     EPSM:     FLOATING POINT NUMBER ERROR (INRNDO)
*     TOLTMR:   TM REMAINDER BOUND MAXIMUM TOLERANCE
*     ITM:      TM COMPUTATION MODE
*     ITMPR:    FLAG OF TM OUTPUT OPTION
*
*-----ROUNDING --------------------------------------------------------------
      INTEGER INTINI
      DOUBLE PRECISION RINUP(2),RINDN(2),FINSR
      COMMON /INCOM/ RINUP,RINDN,FINSR,INTINI
*----------------------------------------------------------------------------
*
*     INTINI:   FLAG OF ROUNDING PARAMETER COMPUTATION
*
*-----DATA FOR DATRN, TMTRNF ------------------------------------------------
      PARAMETER(LBC=2000)
      INTEGER IBINO(0:LNO,0:LNO),ICBIN(LBC),ITRD1(0:LNO),ITRD2(0:LNO),
     *        ITREC(0:LNO),ITREA(0:LNO),JJTR(LNV),ICALLT
      DOUBLE PRECISION TRAR(0:LNO),TRCR(0:LNO),TRINUP,CBIN(LBC)
      COMMON /TRNCOM/ TRAR,TRCR,TRINUP,CBIN,IBINO,ICBIN,
     *       ITRD1,ITRD2,ITREC,ITREA,JJTR,ICALLT
*----------------------------------------------------------------------------
*
*     IBINO:    BINOMIAL COEFFICIENTS
*     CBIN:     SUPPLEMENT TO BINOMIAL COEFFICIENTS IF THE VALUE IS BEYOND
*               THE LIMIT OF INTEGER NUMBERS
*     ICBIN:    COUNTING OF EXTRA OPERATIONS FOR CBIN
*     ICALLT:   FLAG OF INITIAL DATRN, TMTRN CALLS
*
*-----NON-ABORTING ARITHMETIC -----------------------------------------------
      INTEGER LARI
      COMMON /ARIST/ LARI
*----------------------------------------------------------------------------
*
*     LARI:     FLAG OF NON-ABORTING ARITHMETIC
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
*     LGUI:     FLAG OF GUI MODE  1: FULL GUI ON   0: OFF --> ASCII GUI MODE
*     NGVAR:    NUMBER OF VARIABLE NUMBERS CURRENTLY IN IGVAR FOR EACH WINDOW
*               NGVAR(IWND,K)
*     IGVAR:    VARIABLE NUMBERS BOUND TO READ COMMANDS FOR EACH WINDOW
*               IGVAR(IWND,IGV,K)
*               IWND: THE GUI WINDOW NUMBER
*               IGV:  THE NUMBER OF THE READ COMMAND
*               K:    WINDOW MODE    K=1: HIDDEN    K=2: CURRENT
*
*%%
*
*     INITIALIZING THE VARIABLES IN THE COMMON BLOCKS
*     ***********************************************
*
      DATA EPS    / 1.D-20 /
      DATA EPSMAC / 1.D-99 /
      DATA TOLTMR / 1.D10  /
*
      DATA LEW    / 0 /
      DATA LEWI   / 0 /
      DATA LFLT   / 0 /
      DATA ITM    / 0 /
      DATA ITMPR  / 0 /
      DATA INTINI / -1 /
      DATA FINSR  / 10.D0 /
      DATA ICALLT / 0 /
      DATA LARI   / 1 /
      DATA LGUI   / 0 /
      DATA NGVAR  / LGW*0, LGW*0 /
*
      END
*%%
*
      SUBROUTINE FOXPAR(KMEM,KVAR,KDIM)
*     *********************************
*
*     THIS SUBROUTINE VERIFIES THAT THE PARAMETERS IN A "MEMORY MANAGEMENT"
*     COMMON BLOCK OUTSIDE THIS PACKAGE AGREE WITH THOSE HERE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF((KMEM.NE.LMEM).OR.(KVAR.NE.LVAR).OR.(KDIM.NE.LDIM)) THEN
         PRINT*,'@@@ ERROR, WRONG MEMORY MANAGEMENT '//
     *          'PARAMETERS IN DAFOX '
         CALL FOXSTP(1)
      ENDIF
*
      RETURN
      END
*%%
      SUBROUTINE FOXALL(IC,L,NLEN)
*     ****************************
*
*     THIS SUBROUTINE CREATES A NEW VARIABLE AND ALLOCATES SPACE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      INTEGER IC(*)
*
      CALL VARCHK(IVAR+L)
      CALL MEMCHK(IMEM+L*NLEN)
      IF(NLEN.LT.1) THEN
         PRINT*,'$$$ ERROR IN FOXALL, VARIABLE LENGTH IS TOO SMALL'
         CALL FOXDEB
      ENDIF
*
      DO 10 I=1,L
*
      IVAR = IVAR+1
      IC(I) = IVAR
*
      NTYP(IVAR) = NRE
      NBEG(IVAR) = IMEM+1
      NEND(IVAR) = NBEG(IVAR)
      NMAX(IVAR) = IMEM+NLEN
      IMEM = IMEM+NLEN
*
  10  CONTINUE
*
      RETURN
      END
*%%
      SUBROUTINE FOXDAL(IDAL,L)
*     *************************
*
*     THIS SUBROUTINE DEALLOCATES THE L VARIABLES IDAL
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      INTEGER IDAL(*)
*
      DO 10 I=L,1,-1
*
      IF(IDAL(I).EQ.IVAR) THEN
         IVAR = IVAR-1
         IMEM = NMAX(IVAR)
      ENDIF
*
  10  CONTINUE
*
      RETURN
      END
*%%
      SUBROUTINE FOXAAL(IC,L,NLEN)
*     ****************************
*
*     THIS SUBROUTINE CREATES A NEW ONE-DIMENSIONAL ARRAY WITH LENGTH L
*     AND ALLOCATES SPACE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      INTEGER LD(1)
*
      LD(1) = L
      CALL FOXAAM(IC,LD,1,NLEN)
*
      RETURN
      END
*%%
      SUBROUTINE FOXAAM(IC,LD,N,NLEN)
*     *******************************
*
*     THIS SUBROUTINE CREATES AN N-DIMENSIONAL ARRAY. THE LENGHT OF EACH
*     DIMENSION IS GIVEN IN LD AND NLEN IS THE SIZE OF EACH ENTRY
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      INTEGER LD(*)
*
      IF(NLEN.LT.1) THEN
         PRINT*,'$$$ ERROR IN FOXAAM, VARIABLE LENGTH IS TOO SMALL'
         CALL FOXDEB
      ENDIF
*
      IF(IDIM+N.GE.LDIM) THEN
         PRINT*,'!!! RUNTIME MEMORY EXHAUSTION, INCREASE LDIM'
         CALL FOXSTL
      ENDIF
*
      JNUM  = 1
      IDIM  = IDIM+1
      IDIMO = IDIM
      DO 100 J=1,N
      IDIM = IDIM+1
      NDIM(IDIM) = LD(J)
      JNUM = JNUM*NDIM(IDIM)
 100  CONTINUE
*
      NDIM(IDIMO) = JNUM
      IVAR = IVAR+1
      NTYP(IVAR) = -N
      NBEG(IVAR) = 0
      NEND(IVAR) = NLEN
      NMAX(IVAR) = IDIMO
      IC = IVAR
*
      CALL VARCHK(IVAR+JNUM)
      CALL MEMCHK(IMEM+JNUM*NLEN)
*
      DO 110 J=1,JNUM
      IVAR = IVAR+1
      NTYP(IVAR) = NRE
      NBEG(IVAR) = IMEM+1
      NEND(IVAR) = IMEM+1
      NMAX(IVAR) = IMEM+NLEN
      CC(IMEM+1) = 0.D0
      IMEM = IMEM+NLEN
 110  CONTINUE
*
      RETURN
      END
*%%
      SUBROUTINE FOXADA(IDAL,L)
*     *************************
*
*     THIS SUBROUTINE DEALLOCATES THE ONE DIMENSIONAL ARRAY WITH LENGTH L
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CALL FOXMAD(IDAL,L,1)
*
      RETURN
      END
*%%
      SUBROUTINE FOXMAD(IDAL,L,N)
*     ***************************
*
*     THIS SUBROUTINE DEALLOCATES AN N-DIMENSIONAL ARRAY WITH A TOTAL
*     NUMBER OF L ENTRIES
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(IDAL+L.NE.IVAR) THEN
         PRINT*,'$$$ ERROR IN FOXMAD, '//
     *          'VARIABLE DEALLOCATION OUT OF ORDER'
         CALL FOXDEB
      ENDIF
*
      DO 10 I=L,1,-1
      IMEM = NBEG(IVAR)-1
      IVAR = IVAR-1
 10   CONTINUE
*
      IVAR = IVAR-1
      IDIM = IDIM-N-1
*
      RETURN
      END
*%%
      INTEGER FUNCTION LOOKUP(IARR,JDIM,N)
*     ************************************
*
*     THIS FUNCTION PERFORMS AN ARRAY LOOKUP. IT RETURNS THE ADDRESS
*     OF THE COMPONENT OF IARR DESCRIBED BY THE N COMPONENTS OF JDIM
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      INTEGER JDIM(*)
*
      INDA  = -NTYP(IARR)
      IIDIM =  NMAX(IARR)
*
      LOOKUP  = 0
*
      IF(INDA.LE.0) THEN
         PRINT*,'$$$ ERROR IN LOOKUP, '//
     *          'INDEXED VARIABLE IS NOT DECLARED ARRAY'
         CALL FOXSTP(1)
      ENDIF
*
      DO 100 J=INDA,2,-1
      IJ = JDIM(J)
      NJ = NDIM(IIDIM+J)
      KJ = NDIM(IIDIM+J-1)
      IF(IJ.LT.1.OR.IJ.GT.NJ) THEN
         PRINT*,'$$$ ERROR IN LOOKUP, ARRAY INDEX ',J,' IS OUT OF BOUND'
         PRINT*,'    INDEX, BOUND = ',IJ,',',NJ
         CALL FOXSTP(1)
      ENDIF
      LOOKUP = LOOKUP+IJ-1
      LOOKUP = LOOKUP*KJ
 100  CONTINUE
*
      IJ = JDIM(1)
      NJ = NDIM(IIDIM+1)
      IF(IJ.LT.1.OR.IJ.GT.NJ) THEN
         PRINT*,'$$$ ERROR IN LOOKUP, ARRAY INDEX ',J,' IS OUT OF BOUND'
         PRINT*,'    INDEX, BOUND = ',IJ,',',NJ
         CALL FOXSTP(1)
      ENDIF
      LOOKUP = LOOKUP+IJ+IARR
*
      RETURN
      END
*%%
      SUBROUTINE FOXTYP(INA,INC)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = NTYP(INA)
      RETURN
      END
*
      SUBROUTINE FOXNTY(INB)
*     **********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      PARAMETER(LTYP=11,LSYNR=30,LOPS=20,LFUN=100,LSUB=200)
      CHARACTER CTYID(LTYP)*2
      COMMON /TYCID/ CTYID
*
      PRINT*,'$$$ ERROR, VARIABLE ',INB,' HAS WRONG TYPE: '
     *        //CTYID(NTYP(INB))
      CALL FOXDEB
      RETURN
      END
*
      SUBROUTINE FOXLEN(INA,INC)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      IF(NTYP(INA).EQ.NRE) THEN
         CC(NBEG(INC)) = 1
      ELSE
         CC(NBEG(INC)) = NEND(INA)-NBEG(INA)+1
      ENDIF
      RETURN
      END
*
      SUBROUTINE FOXMEM(INA,INC)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = NBEG(INA)
      RETURN
      END
*
      SUBROUTINE FOXPOI(INA,INC)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = INA
      RETURN
      END
*
      SUBROUTINE MEMCHK(JMEM)
*     ***********************
*
*     CHECKS FOR MAIN MEMORY OVERFLOW
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      COMMON /MSTAT/ MMEM,MVAR
*
      MMEM = MAX(MMEM,JMEM)
*
      IF(JMEM.LE.LMEM) RETURN
*
      WRITE(6,'(1X,A)') '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LMEM'
      CALL FOXSTL
      END
*
      SUBROUTINE VARCHK(JVAR)
*     ***********************
*
*     CHECKS FOR NUMBER OF VARIABLE OVERFLOW
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      COMMON /MSTAT/ MMEM,MVAR
*
      MVAR = MAX(MVAR,JVAR)
*
      IF(JVAR.LE.LVAR) RETURN
*
      WRITE(6,'(1X,A)') '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LVAR'
      CALL FOXSTL
      END
*
      SUBROUTINE OPENF(INUNIT,INAME,ISTAT)
*     ************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER NAME*200,STAT*100
*
      IF(NTYP(INUNIT).NE.NRE) CALL FOXNTY(INUNIT)
      IF(NTYP(INAME).NE.NST) CALL FOXNTY(INAME)
      IF(NTYP(ISTAT).NE.NST) CALL FOXNTY(ISTAT)
*
      IUNIT = NINT(CC(NBEG(INUNIT)))
*
      CALL STSTR(INAME,NAME,LNAME)
      CALL STSTR(ISTAT,STAT,IL)
      LSTAT = 0
      DO 10 I=1,IL
      IF(STAT(I:I).NE.' ') THEN
         LSTAT = LSTAT+1
         STAT(LSTAT:LSTAT) = STAT(I:I)
      ENDIF
 10   CONTINUE
      CALL CAPSTR(STAT,1,LSTAT)
*
      IF(LSTAT.EQ.7.AND.STAT(1:7).EQ.'SCRATCH') THEN
         OPEN(IUNIT,STATUS='SCRATCH')
      ELSE
         OPEN(IUNIT,FILE=NAME(1:LNAME),STATUS=STAT(1:LSTAT))
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE OPENFB(INUNIT,INAME,ISTAT)
*     *************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER NAME*200,STAT*100
*
      IF(NTYP(INUNIT).NE.NRE) CALL FOXNTY(INUNIT)
      IF(NTYP(INAME).NE.NST) CALL FOXNTY(INAME)
      IF(NTYP(ISTAT).NE.NST) CALL FOXNTY(ISTAT)
*
      IUNIT = NINT(CC(NBEG(INUNIT)))
*
      CALL STSTR(INAME,NAME,LNAME)
      CALL STSTR(ISTAT,STAT,IL)
      LSTAT = 0
      DO 10 I=1,IL
      IF(STAT(I:I).NE.' ') THEN
         LSTAT = LSTAT+1
         STAT(LSTAT:LSTAT) = STAT(I:I)
      ENDIF
 10   CONTINUE
      CALL CAPSTR(STAT,1,LSTAT)
*
      IF(LSTAT.EQ.7.AND.STAT(1:7).EQ.'SCRATCH') THEN
         OPEN(IUNIT,STATUS='SCRATCH',FORM='UNFORMATTED')
      ELSE
         OPEN(IUNIT,FILE=NAME(1:LNAME),STATUS=STAT(1:LSTAT),
     *           FORM='UNFORMATTED')
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE CLOSEF(IUNIT)
*     ************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(IUNIT).NE.NRE) CALL FOXNTY(IUNIT)
*
      IU = NINT(CC(NBEG(IUNIT)))
      CLOSE(IU)
*
      RETURN
      END
*
      SUBROUTINE REWF(IUNIT)
*     **********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(IUNIT).NE.NRE) CALL FOXNTY(IUNIT)
*
      IU = NINT(CC(NBEG(IUNIT)))
      REWIND(IU)
*
      RETURN
      END
*
      SUBROUTINE BACKF(IUNIT)
*     ***********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(IUNIT).NE.NRE) CALL FOXNTY(IUNIT)
*
      IU = NINT(CC(NBEG(IUNIT)))
      BACKSPACE(IU)
*
      RETURN
      END
*
      SUBROUTINE READS(IUNIT,INA)
*     ***************************
*
*     THIS SUBROUTINE READS A STRING FROM UNIT IUNIT TO INA WITHOUT THE USUAL
*     CONVERSION TO REAL PERFORMED IN READ.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER A*4096
*
      IF(NTYP(IUNIT).NE.NRE) CALL FOXNTY(IUNIT)
      IU = NINT(CC(NBEG(IUNIT)))
      IA = NBEG(INA)
*
      READ(IU,'(A)',END=20,ERR=20) A
      LA = ILAST(A,1,4096)
      IF(LA.EQ.0) THEN
*        DO NOT RETURN EMPTY STRINGS
         LA = 1
         NC(IA) = ICHAR(' ')
      ELSE
         LA = MIN(LA,NMAX(INA)-IA+1)
         DO 10 I=1,LA
         NC(IA+I-1) = ICHAR(A(I:I))
  10     CONTINUE
      ENDIF
      NTYP(INA) = NST
      NEND(INA) = IA+LA-1
*
      RETURN
*
*     ERROR AND END OF FILE HANDLING, RETURN AN EMPTY STRING
  20  NTYP(INA) = NST
      NEND(INA) = IA-1
      RETURN
*
      END
*
      SUBROUTINE READB(IUNIT,INA)
*     ***************************
*
*     THIS SUBROUTINE READS A COSY OBJECT FROM UNIT IUNIT TO INA IN BINARY FORM.
*
*     NOTE: AN ARITHMETIC FAILURE PROCESSING SUCH AS ABORTING IS NOT APPLIED.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----ROUNDING --------------------------------------------------------------
      INTEGER INTINI
      DOUBLE PRECISION RINUP(2),RINDN(2),FINSR
      COMMON /INCOM/ RINUP,RINDN,FINSR,INTINI
*----------------------------------------------------------------------------
*
      DIMENSION IEWIN(LNV)
      DATA IDBIN / 200907 /
*
      IF(NTYP(IUNIT).NE.NRE) CALL FOXNTY(IUNIT)
      IU = NINT(CC(NBEG(IUNIT)))
*
      READ(IU,ERR=90) IDBINI
      IF(IDBINI.NE.IDBIN) THEN
         PRINT*,'$$$ ERROR IN READB, COSY BINARY DATA IS INCOMPATIBLE.'
         PRINT*,'    CHECK THE CORRECTNESS OF THE BINARY DATA FILE,'
         PRINT*,'    OR UPDATE THE BINARY DATA FILE.'
         CALL FOXDEB
      ENDIF
*
      READ(IU,ERR=90) NTY,LENA
      IF(NMAX(INA)-NBEG(INA)+1.LT.LENA) CALL FOXERV('READB')
C
C     THE FOLLOWING CHECK MAY BE USEFUL, BUT PROBABLY AN OVERKILL.
C
C     IF(NTY.LE.0.OR.NTY.GE.12) THEN
C        PRINT*,'$$$ ERROR IN READB, ILLEGAL DATA TYPE ',NTY
C        CALL FOXDEB
C     ENDIF
C
      NEND(INA) = NBEG(INA)+LENA-1
      NTYP(INA) = NTY
*
      DO 10 I=NBEG(INA),NEND(INA)
      READ(IU,ERR=90) NC(I),CC(I)
 10   CONTINUE
*
C     IF(NTY.EQ.NDA.OR.NTY.EQ.NCD.OR.NTY.EQ.NTM) THEN
      IF(NTY.EQ.NDA.OR.NTY.EQ.NCD) THEN
         READ(IU,ERR=90) NOMAXI,NVMAXI,LEWIN
         READ(IU,ERR=90) (IEWIN(I),I=1,NVMAXI)
*
         IF(NOMAXI.NE.NOMAX.OR.NVMAXI.NE.NVMAX) THEN
            PRINT*,'$$$ ERROR IN READB, DAINI SETUP IS INCOMPATIBLE.'
            PRINT*,'    DAINI: ORDER, DIMENSION = ',NOMAX,NVMAX
            PRINT*,'    BINARY DATA: ',NOMAXI,NVMAXI
            CALL FOXDEB
         ELSEIF(LEWIN.NE.LEW) THEN
            PRINT*,'$$$ ERROR IN READB, DA WEIGHT IS INCOMPATIBLE.'
            PRINT*,'    DA WEIGHT: (0: OFF, 1: ON) = ',LEW
            PRINT*,'    BINARY DATA: ',LEWIN
            CALL FOXDEB
         ELSEIF(LEW.EQ.1) THEN
            DO 20 I=1,NVMAX
            IF(IEWIN(I).NE.IEW(I)) THEN
               PRINT*,'$$$ ERROR IN READB, DA WEIGHT MISMATCHES.'
               PRINT*,'    DA WEIGHT (',I,') = ',IEW(I)
               PRINT*,'    WEIGHT IN BINARY DATA: ',IEWIN(I)
               CALL FOXDEB
            ENDIF
 20         CONTINUE
         ENDIF
C
C        THE FOLLOWING PROCESS IS TO MAKE IT COMPATIBLE WITH TMREAB.
C        HOWEVER, THE EQUIVALENT ARI PROCESSING IS NOT DONE HERE.
C        THUS, RATHER MAKE IT CONSISTENT BETWEEN DA, CD AND TM.
C
C        IF(NTY.EQ.NTM) THEN
C           NOTM = NC(NEND(INA)-2)-1
C           IF(NOTM.GT.NOCUT) THEN
C              IF(INTINI.EQ.-1) CALL INRNDO(1.D0,1.D0,.FALSE.,.FALSE.,1)
C              INTINO = 1
C              IMD = NINT(CC(NEND(INA)))
C              IF(IMD.GT.0.AND.INTINI.NE.1) THEN
C                 INTINO = INTINI
C                 INTINI = 1
C              ENDIF
C              CALL TMNOT(INA,NOCUT,INA)
C              IF(INTINO.NE.1) INTINI = INTINO
C           ENDIF
C        ENDIF
*
      ENDIF
*
      RETURN
*
 90   CONTINUE
      PRINT*,
     *   '$$$ ERROR IN READB, BINARY DATA CANNOT BE READ SUCCESSFULLY.'
      CALL FOXDEB
*
      END
*
      SUBROUTINE WRITEB(IUNIT,INA)
*     ****************************
*
*     THIS SUBROUTINE WRITES INA TO UNIT IUNIT IN BINARY FORM.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DATA IDBIN / 200907 /
*
      IF(NTYP(IUNIT).NE.NRE) CALL FOXNTY(IUNIT)
      IU = NINT(CC(NBEG(IUNIT)))
*
      NTY = NTYP(INA)
      LENA = NEND(INA)-NBEG(INA)+1
*
      WRITE(IU) IDBIN
      WRITE(IU) NTY,LENA
*
      DO 10 I=NBEG(INA),NEND(INA)
      WRITE(IU) NC(I),CC(I)
 10   CONTINUE
*
C     IF(NTY.EQ.NDA.OR.NTY.EQ.NCD.OR.NTY.EQ.NTM) THEN
      IF(NTY.EQ.NDA.OR.NTY.EQ.NCD) THEN
         WRITE(IU) NOMAX,NVMAX,LEW
         WRITE(IU) (IEW(I),I=1,NVMAX)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE READM(INA,INFO,INCL,IMCC,IMNC,IDAPR)
*     ***********************************************
*
*     THIS SUBROUTINE READS A COSY OBJECT INA FROM ARRAYS IMCC AND IMNC
*     WITH LENGTH INCL. THIS IS THE COUNTERPART OF WRITEM.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DATA IDMEM / 201309 /
*
      IF(NTYP(INFO).NE.NVE) CALL FOXNTY(INFO)
      IF(NTYP(INCL).NE.NRE) CALL FOXNTY(INCL)
      IJ = NBEG(INFO)
*
      IF(NTYP(IMCC).GT.0.OR.NTYP(IMNC).GT.0) THEN
         PRINT*,'$$$ ERROR IN READM, 4TH OR 5TH ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(NEND(INFO)-IJ+1.LT.3) THEN
         PRINT*,'$$$ ERROR IN READM, 2ND ARGUMENT (VE) TOO SHORT.'
         CALL FOXDEB
      ENDIF
*
      NTY    = NINT(CC(IJ))
      LENA   = NINT(CC(IJ+1))
      IDMEMI = NINT(CC(IJ+2))
      NCL = NINT(CC(NBEG(INCL)))
*
      NEND(INA) = NBEG(INA)+LENA-1
      NTYP(INA) = NTY
*
      IF(IDMEMI.NE.IDMEM) THEN
         PRINT*,'$$$ ERROR IN READM, THE VERSION ID',IDMEMI,
     *          ' INCOMPATIBLE.'
         PRINT*,'    UPDATE TO THE CURRENT VERSION ',IDMEM,' BY WRITEM.'
         CALL FOXDEB
      ELSEIF(NTY.LT.1.OR.NTY.GT.11) THEN
         PRINT*,'$$$ ERROR IN READM, ILLEGAL DATA TYPE ',NTY
         CALL FOXDEB
      ELSEIF(NCL.LT.LENA) THEN
         PRINT*,'$$$ ERROR IN READM, '//
     *          'ARRAY LENGTH (3RD ARGUMENT) TOO SHORT.'
         CALL FOXDEB
      ELSEIF(NMAX(INA).LT.NEND(INA)) THEN
         CALL FOXERV('READM')
C     ELSEIF(NTY.EQ.NDA.OR.NTY.EQ.NCD.OR.NTY.EQ.NTM) THEN
      ELSEIF(NTY.EQ.NDA.OR.NTY.EQ.NCD) THEN
         CALL DACHK(IDAPR,'READM')
      ENDIF
*
      IA = NBEG(INA)
      DO 10 I=1,LENA
      CC(IA) = CC(NBEG(IMCC+I))
      NC(IA) = NINT(CC(NBEG(IMNC+I)))
      IA = IA+1
 10   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE WRITEM(INA,INFO,INCL,IMCC,IMNC,IDAPR)
*     ************************************************
*
*     THIS SUBROUTINE WRITES INA TO ARRAYS IMCC AND IMNC WITH LENGTH INCL.
*     THE VARIABLE INFORMATION INFO (VE) CONSISTS OF THE DATA TYPE,
*     THE LENGTH IN THE COSY MEMORY CC() AND NC(), THE WRITEM VERSION ID.
*     IF INA IS DA, TM, CD, THE DA INFORMATION IDAPR CONSISTING OF
*     NOMAX, NVMAX, IEW(J),J=1,NVMAX IS RETURNED (VE); ELSE 0 (RE) IS RETURNED.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DATA IDMEM / 201309 /
*
      IF(NTYP(INCL).NE.NRE) CALL FOXNTY(INCL)
      NCL = NINT(CC(NBEG(INCL)))
      NTY = NTYP(INA)
      LENA = NEND(INA)-NBEG(INA)+1
*
      IF(NTYP(IMCC).GT.0.OR.NTYP(IMNC).GT.0) THEN
         PRINT*,'$$$ ERROR IN WRITEM, 4TH OR 5TH ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(NBEG(IMCC+1).EQ.NBEG(IMNC+1)) THEN
         PRINT*,'$$$ ERROR IN WRITEM, THE SAME ARRAY IS USED FOR THE '
     *        //'4TH AND 5TH ARGUMENTS.'
         CALL FOXDEB
      ELSEIF(NCL.LT.LENA) THEN
         PRINT*,'$$$ ERROR IN WRITEM, '//
     *          'ARRAY LENGTH (3RD ARGUMENT) TOO SHORT.'
         CALL FOXDEB
      ENDIF
*
      NTYP(INFO) = NVE
      IJ = NBEG(INFO)
      CC(IJ) = NTY
      IJ = IJ+1
      CC(IJ) = LENA
      IJ = IJ+1
      CC(IJ) = IDMEM
      NEND(INFO) = IJ
      IF(IJ.GT.NMAX(INFO)) CALL FOXERV('WRITEM')
*
      IA = NBEG(INA)
      DO 10 I=1,LENA
      NTYP(IMCC+I) = NRE
      NTYP(IMNC+I) = NRE
      CC(NBEG(IMCC+I)) = CC(IA)
      CC(NBEG(IMNC+I)) = NC(IA)
      IA = IA+1
 10   CONTINUE
*
      ID = NBEG(IDAPR)
C     IF(NTY.EQ.NDA.OR.NTY.EQ.NCD.OR.NTY.EQ.NTM) THEN
      IF(NTY.EQ.NDA.OR.NTY.EQ.NCD) THEN
         NTYP(IDAPR) = NVE
         CC(ID) = NOMAX
         ID = ID+1
         CC(ID) = NVMAX
         IF(LEW.NE.0) THEN
            DO 20 J=1,NVMAX
            ID = ID+1
            CC(ID) = IEW(J)
 20         CONTINUE
         ENDIF
         IF(ID.GT.NMAX(IDAPR)) CALL FOXERV('WRITEM')
         NEND(IDAPR) = ID
      ELSE
         NTYP(IDAPR) = NRE
         CC(ID) = 0.D0
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DACHK(IDAPR,NAME)
*     ****************************
*
*     THIS SUBROUTINE CHECKS THE DA SETUP SUPPLIED IN IDAPR
*     BY THE ROUTINE NAME WITH THE CURRENT DAINI SETUP.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION IEWIN(LNV)
      CHARACTER NAME*(*),STR*20
*
      IF(NTYP(IDAPR).NE.NVE) CALL FOXNTY(IDAPR)
*
      IDB = NBEG(IDAPR)
      IDE = NEND(IDAPR)
      ID = IDB
      NOMAXI = NINT(CC(ID))
      ID = ID+1
      NVMAXI = NINT(CC(ID))
*
      LN = LEN(NAME)
      STR(1:LN) = NAME
*
      IF(NOMAXI.NE.NOMAX.OR.NVMAXI.NE.NVMAX) THEN
         WRITE(6,'(A/A,2I4,A,2I4)') 
     *      ' $$$ ERROR IN '//STR(1:LN)//', DAINI SETUP INCOMPATIBLE.',
     *      '     DAINI: ORDER, DIMENSION =',NOMAX,NVMAX,
     *      ',   '//STR(1:LN)//':',NOMAXI,NVMAXI
         CALL FOXDEB
      ENDIF
*
      DO 10 I=1,NVMAX
      IEWIN(I) = 1
 10   CONTINUE
*
      LD2 = IDE-IDB-1
      IF(LD2.GE.1) THEN
         DO 20 I=1,MIN(NVMAX,LD2)
         ID = ID+1
         IEWIN(I) = NINT(CC(ID))
 20      CONTINUE
      ENDIF
*
      DO 30 I=1,NVMAX
      IF(IEWIN(I).NE.IEW(I)) THEN
         WRITE(6,'(A/A,I3,A,I4,A,I4)') 
     *      ' $$$ ERROR IN '//STR(1:LN)//', DA WEIGHT MISMATCHES.',
     *      '     DAINI DA WEIGHT (',I,'):',IEW(I),
     *      ',   '//STR(1:LN)//' WEIGHT:',IEWIN(I)
         CALL FOXDEB
      ENDIF
 30   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE MEMMAX(INC)
*     **********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      COMMON /MSTAT/ MMEM,MVAR
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = MMEM
      RETURN
      END
*
      SUBROUTINE MEMALL(INC)
*     **********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = IMEM
      RETURN
      END
*
      SUBROUTINE MEMFRE(INC)
*     **********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = LMEM-IMEM
      RETURN
      END
*
      SUBROUTINE VARMAX(INC)
*     **********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      COMMON /MSTAT/ MMEM,MVAR
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = MVAR
      RETURN
      END
*
      SUBROUTINE VARALL(INC)
*     **********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = IVAR
      RETURN
      END
*
      SUBROUTINE CPUSEC(INC)
*     **********************
*
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      REAL TIM(2)
      NTYP(INC) = NRE
      CC(NBEG(INC)) = 0.D0
      CC(NBEG(INC)) = ETIME(TIM)
*     * CPUTIME COUNTING SUBROUTINE/FUNCTION DIFFERS BY SYSTEM.
*     * REFER TO THE SYSTEM MANUAL.
C     COSY-PLOOP ADD
C
C     RTC: REAL TIME CLOCK FOR IBM XLF FORTRAN COMPILER
C     USE RTC INSTEAD OF ETIME FOR SEABORG AT NERSC
C
C     CC(NBEG(INC)) = RTC()
C
      RETURN
      END
*
      SUBROUTINE SLEEPM(IMSEC)
*     ************************
*
*     THIS SUBROUTINE SUSPENDS EXECUTION OF THE PROGRAM FOR IMSEC (M SEC)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(IMSEC).NE.NRE) CALL FOXNTY(IMSEC)
*
*     SLEEPQQ is an Intel Fortran specific, platform independent call.
*
      CALL SLEEPQQ(NINT(CC(NBEG(IMSEC))))
*
*     SLEEP is not compiler specific and platform independent,
*     but only does full seconds, not milliseconds
*
C     CALL SLEEP(NINT(1.D-3*CC(NBEG(IMSEC))))
*
      RETURN
      END
*
      SUBROUTINE PWTIME(INC)
*     **********************
*
*     THIS SUBROUTINE RETURNS THE ELAPSED WALL-CLOCK TIME (SEC) ON MPI.
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION MPI_WTIME
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = 0.D0
*MPI  CC(NBEG(INC)) = MPI_WTIME()                                        *MPI
      CALL CPUSEC(INC)                                                   *NORM
      RETURN
      END
*
      SUBROUTINE PNPRO(INC)
*     *********************
*
*     THIS SUBROUTINE RETURNS THE NUMBER OF PROCESSES (TASKS) IN A GROUP ON MPI.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      INTEGER MTASK,ISROOT
      COMMON /PLCOM/ MTASK,ISROOT
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = DBFLOAT(MTASK)
      RETURN
      END
*
      SUBROUTINE PROOT(INC)
*     *********************
*
*     THIS SUBROUTINE RETURNS THE PROCESS NUMBER OF THE ROOT PROCESS ON MPI.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      INTEGER MTASK,ISROOT
      COMMON /PLCOM/ MTASK,ISROOT
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = DBFLOAT(ISROOT)
      RETURN
      END
*
      SUBROUTINE QUIT(I)
*     ******************
*
*     THIS SUBROUTINE TERMINATES EXECUTION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      IF(LGUI.EQ.0.OR.NINT(CC(NBEG(I))).NE.0)
     *   WRITE(6,'(1X,A8,I10)') ' QUIT AT',NINT(CC(NBEG(I)))
      CALL FOXSTP(0)
      END
*
      SUBROUTINE CRASH
*     ****************
*
*     THIS SUBROUTINE CRASHES THE SYSTEM IN ORDER TO PRODUCE A TRACEBACK
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*MPI  CALL FOXSTP(1)                                                     *MPI
*
      X = -1.D0
      X = X*X*X
      IF(X.LT.0.D0) X = X*X*X
      PRINT*,SQRT(X)
      CALL FOXSTP(1)
      END
*
      SUBROUTINE OS(ISYS)
*     *******************
*
*     THIS SUBROUTINE CALLS A SYSTEM CALL
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER SNAME*1000
*
      IF(NTYP(ISYS).NE.NST) CALL FOXNTY(ISYS)
*
      IS = 0
      DO 10 I=NBEG(ISYS),NEND(ISYS)
      IS = IS+1
 10   SNAME(IS:IS) = CHAR(NC(I))
*
      CALL SYSTEM(SNAME(1:IS))
*
      RETURN
      END
*
      SUBROUTINE ARGGET(INA,INC)
*     **************************
*
*     THIS SUBROUTINE INTERFACES FROM COSY PROGRAM TO FORTRAN INTRINSIC GETARG.
*     IT RETURNS THE COMMAND LINE ARGUMENT IN INC SPECIFIED BY THE POSITION
*     NUMBER OF THE ARGUMENT IN INA. 0 IS FOR THE COMMAND ITSELF.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER A*4096
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      CALL GETARG( NINT(CC(NBEG(INA))) ,A)
      IL = ILAST(A,1,NMAX(INC)-NBEG(INC)+1)
*
      IC = NBEG(INC)-1
      DO 10 I=1,IL
 10   NC(IC+I) = ICHAR(A(I:I))
      NEND(INC) = IC+IL
      NTYP(INC) = NST
*
      RETURN
      END
*
      SUBROUTINE FOXERV(NAME)
*     ***********************
*
*     THIS SUBROUTINE REPORTS THE VARIABLE EXHAUSTION OF INA AND STOPS
*
      CHARACTER NAME*(*),STR*20
*
      LN = LEN(NAME)
      STR(1:LN) = NAME
      PRINT*,'$$$ ERROR IN '//STR(1:LN)//', VARIABLE EXHAUSTED'
      CALL FOXDEB
*
      RETURN
      END
*
      SUBROUTINE SCRLEN(INA)
*     **********************
*
*     THIS SUBROUTINE SETS THE SPACE A SRATCH VARIABLE IS ALLOCATED WITH
*     IF THE PARAMETER IS NEGATIVE, IT RETURNS THE MOMENTARY VALUE
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      COMMON /CODEX/ ILIS,NSCR
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      NS = NINT(CC(IA))
*
      IF(NS.GE.0) THEN
         NSCR = NS
      ELSE
         CC(IA) = DBFLOAT(NSCR)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE FOXKEY(CFOX)
*     ***********************
*
*     THIS SUBROUTINE CHECKS THAT THE PRECOMPILER VERSION AND THE DA PACKAGE
*     VERSION MATCH.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER CFOX*9,CDA*9
*
      DATA CDA / 'FOX V10.0' /
*
      IF(CFOX.NE.CDA) THEN
         PRINT*,'@@@ ERROR, FOXY VERSION AND LIBRARY INCOMPATIBLE'
         CALL FOXSTP(1)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE FOXDEB
*     *****************
*
*     THIS SUBROUTINE SERVES AS A DEBUGGING TOOL.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      COMMON /CODEX/ ILIS,NSCR
*
      ISTOP = 1
      IF(ISTOP.EQ.0) RETURN
*
*     PRODUCING SYSTEM ERROR TRACEBACK BY PERFORMING ILLEGAL OPERATION
*
      PRINT*,' '
      IF(ILIS.LT.0) PRINT*,'    ERROR OCCURED IN .LIS LINE ',-ILIS
      PRINT*,'    PRODUCING TRACEBACK BY DELIBERATE ILLEGAL OPERATION'
     *     //' SQRT(-1.D0)'
      PRINT*,' '
      CALL CRASH
*
      END
*
      SUBROUTINE FOXREA(INA,IUNIT)
*     ****************************
*
*     THIS SUBROUTINE READS INA FROM IUNIT.
*     GUI UNITS ARE SUPPORTED AND SPECIAL HANDLING IS DONE.
*
      IF(IUNIT.LE.-201.AND.IUNIT.GE.-210) THEN
         CALL GUIREA(INA,IUNIT)
      ELSEIF(IUNIT.EQ.-999) THEN
         CALL GUICMD('-999 read')
         CALL FOXREAO(INA,5)
      ELSE
         IF(IUNIT.EQ.5) CALL GUICMD('0 read')
         CALL FOXREAO(INA,IUNIT)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE FOXREAO(INA,IUNIT)
*     *****************************
*
*     THIS SUBROUTINE READS INA (RE OR ST) FROM IUNIT.
*     BLANKS AT THE END OF A STRING ARE OMITTED.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER A*4096
*
      IA = NBEG(INA)
      READ(IUNIT,'(A)',END=20,ERR=20) A
      CALL VALCH(A,1,4096,VAL,IER)
      IF(IER.EQ.0) THEN
         NTYP(INA) = NRE
         CC(IA) = VAL
      ELSE
         LA = ILAST(A,1,MIN(4096,NMAX(INA)-IA+1))
         DO 10 I=1,LA
         NC(IA+I-1) = ICHAR(A(I:I))
  10     CONTINUE
         NTYP(INA) = NST
         NEND(INA) = IA+LA-1
      ENDIF
*
      RETURN
*
*     ERROR AND END OF FILE HANDLING, RETURN AN EMPTY STRING
  20  NTYP(INA) = NST
      NEND(INA) = IA-1
*
      RETURN
      END
*
      SUBROUTINE GUIREA(INA,IUNIT)
*     ****************************
*
*     THIS SUBROUTINE READS VALUES FOR INA FROM THE GUI UNIT NUMBER IUNIT
*     USING DELAYED READS.  -201 >= IUNIT >= -210.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      CHARACTER STR*4096,STRE*(2*4096+2),SNR*3
*
      IWND = -200-IUNIT
      WRITE(SNR,'(I3)') IWND
      IV = NGVAR(IWND,1)
*
*   * CLAIM THE LAST GUI FIELD FOR THIS READ IF AVAILABLE AND UNCLAIMED
      IF(IV.GT.0.AND.IGVAR(IWND,IV,1).LT.0) THEN
         IGVAR(IWND,IV,1) = INA
*
*   * NO PREVIOUS GUI FIELD TO READ, SO WE MAKE OUR OWN
      ELSE
         IV = IV+1
         IF(IV.GT.LGV) THEN
            PRINT*,'$$$ ERROR IN GUIREA, TOO MANY FIELDS'
            CALL FOXDEB
         ENDIF
         IGVAR(IWND,IV,1) = INA
         NGVAR(IWND,1) = IV
         CALL GUICMD(SNR//' just')
*
*      * CONVERT TO STRING IF POSSIBLE
         IF(NTYP(INA).EQ.NST .OR. NTYP(INA).EQ.NRE .OR.
     *      NTYP(INA).EQ.NLO .OR. NTYP(INA).EQ.NCM) THEN
            CALL FOXSTR(INA,STR,LS)
            CALL GUIESC(STR,LS,STRE,LSE)
            CALL GUICMD(SNR//' readfield '//STRE(1:LSE))
*
*      * NO CONVERSION POSSIBLE
         ELSE
            CALL GUICMD(SNR//' readfield')
         ENDIF
         CALL GUICMD(SNR//' newline')
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE FOXPRI(INA,IUNIT,IWRT,NWRT)
*     **************************************
*
*     THIS SUBROUTINE PRINTS A COSY OBJECT INA OF A COSY WRITE COMMAND
*     TO IUNIT. INA IS THE IWRT-TH OBJECT AMONG NWRT OBJECTS
*     IN THE COSY WRITE COMMAND.
*     GUI UNITS ARE SUPPORTED AND SPECIAL HANDLING IS DONE.
*
      IF(IUNIT.LE.-201.AND.IUNIT.GE.-210) THEN
         CALL GUIPRI(INA,IUNIT,IWRT,NWRT)
      ELSE
         CALL FOXPRIO(INA,IUNIT)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE FOXPRIO(INA,IUNIT)
*     *****************************
*
*     THIS SUBROUTINE PRINTS A COSY OBJECT INA TO IUNIT
*     BY INVOKING THE GENERIC PRINTS.
*     THIS CAN BE CALLED FROM THE OTHER FORTRAN SUBROUTINES.
*     GUI UNITS ARE SUPPORTED AND SPECIAL HANDLING IS DONE.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CALL PRDIRB(IUNIT,IUNITA)
*
      NTY = NTYP(INA)
      IF(NTY.EQ.NRE) THEN
         CALL REPRI(INA,IUNITA)
      ELSEIF(NTY.EQ.NST) THEN
         CALL STPRI(INA,IUNITA)
      ELSEIF(NTY.EQ.NLO) THEN
         CALL LOPRI(INA,IUNITA)
      ELSEIF(NTY.EQ.NDA) THEN
         CALL DAPRI(INA,IUNITA)
      ELSEIF(NTY.EQ.NVE) THEN
         CALL VEPRI(INA,IUNITA)
      ELSEIF(NTY.EQ.NCM) THEN
         CALL CMPRI(INA,IUNITA)
      ELSEIF(NTY.EQ.NGR) THEN
         CALL GRPRI(INA,IUNITA)
      ELSEIF(NTY.EQ.NCD) THEN
         CALL CDPRI(INA,IUNITA)
C     ELSEIF(NTY.EQ.NTM) THEN
C        CALL TMPRI(INA,IUNITA)
      ELSE
         PRINT*,'$$$ ERROR IN FOXPRIO, ILLEGAL TYPE ',NTY
         CALL FOXDEB
      ENDIF
*
      CALL PRDIRE(IUNIT,IUNITA)
*
      RETURN
      END
*
      SUBROUTINE PRDIRB(IUNIT,IUNITA)
*     *******************************
*
*     THIS SUBROUTINE CONTROLS THE OUTPUT CHANNEL FOR IUNIT TO BEGIN PRINTING.
*     THE ACTUAL UNIT NUMBER IUNITA IS RETURNED TO BE USED FOR PRINTING.
*     FOR GUI OUTPUT, SPECIAL HANDLING FOR PROPER REDIRECTS IS DONE.
*     AFTER PRINTING, PRDIRE(IUNIT,IUNITA) HAS TO BE CALLED.
*
*     THIS ALLOWS PRINTING COMPLEX COSY OBJECTS SUCH AS DA TO THE GUI.
*
      CHARACTER SNR*3
*
      IUNITA = IUNIT
*
      IF(IUNIT.LE.-201.AND.IUNIT.GE.-210) THEN
         IWND = -200-IUNIT
         WRITE(SNR,'(I3)') IWND
         CALL GUICMD(SNR//' redirect')
         IUNITA = 6
      ELSEIF(IUNIT.EQ.-999) THEN
         CALL GUICMD('-999 redirect')
         IUNITA = 6
      ELSEIF(IUNIT.EQ.6) THEN
         CALL GUICMD('0 redirect')
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE PRDIRE(IUNIT,IUNITA)
*     *******************************
*
*     THIS SUBROUTINE CONTROLS THE OUTPUT CHANNEL FOR IUNIT TO END PRINTING.
*     FOR THE UNIT NUMBER IUNIT OUTPUT, DUE TO THE SPECIAL HANDLING FOR GUI,
*     IUNITA, THE ACTUAL UNIT NUMBER USED, MAY DIFFER FROM IUNIT.
*     THIS SHOULD BE CALLED AS A PAIR WITH PRDIRB(IUNIT,IUNITA).
*
      IF(IUNITA.NE.IUNIT) CALL GUICMD('0 redirect')
*
      RETURN
      END
*
      SUBROUTINE GCMD(INA,IDC,CST,NCST)
*     *********************************
*
*     THIS SUBROUTINE IDENTIFIES IF INA IS A STRING FOR A COSY-GUI COMMAND.
*
*     A COSY-GUI COMMAND IS A STRING
*       - STARTING WITH A BACKSLASH '\'
*       - FOLLOWED BY ALPHABETICAL LETTERS FORMING ONE WORD FOR THE COMMAND
*     ANY TRAILING SPACE CHARACTERS WILL BE IGNORED.
*
*     IF INA REPRESENTS A COSY-GUI COMMAND,
*     THE COMMAND NAME IS RETURNED BY CST(1:NCST),
*     WHERE THE LENGTH OF CST IS LIMITED BY 256. IF THE COMMAND IS
*     LONGER THAN 256, THE FIRST 256 LETTERS ARE RETURNED.
*     IDC IS THE IDENTIFICATION CODE,
*         1: A GUI COMMAND THAT BEGINS WITH R FOR READING
*         0: A GUI COMMAND THAT DOES NOT BEGIN WITH R
*        -1: NOT A GUI COMMAND
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----CHARACTER HANDLING ----------------------------------------------------
      INTEGER ICSA,ICSZ,IDIF,ICCA,ICCZ
      COMMON /TRCHAR/ ICSA,ICSZ,IDIF,ICCA,ICCZ
      CHARACTER SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
      COMMON /SCHAR/ SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
*     SCNUL={NUL}, SCHT={HORIZONTAL TAB}, SCLF={LINE FEED},
*     SCCR={CARRIAGE RETURN}, SCSP={SPACE}=' ',
*     SCQM={QUOTATION MARK}='"', SCBSL={BACK SLASH}='\'
*----------------------------------------------------------------------------
*
      CHARACTER CST*256,A*4096
*
      IDC = -1
*
      IF(NTYP(INA).NE.NST.OR.NBEG(INA).GE.NEND(INA)) RETURN
*
      CALL STSTR(INA,A,LA)
      IF(A(1:1).NE.SCBSL) RETURN
*
      LA = ILAST(A,1,LA)
      DO 10 I=2,LA
      ICI = ICHAR(A(I:I))
      IF(ICI.GE.ICCA.AND.ICI.LE.ICCZ) THEN
         ICI = ICI-IDIF
         A(I:I) = CHAR(ICI)
      ENDIF
      IF(ICI.LT.ICSA.OR.ICI.GT.ICSZ) RETURN
  10  CONTINUE
*
      IDC = 0
      NCST = LA-1
*
      IF(NCST.GT.256) THEN
         CST = A(2:256+1)
      ELSE
         CST = '                    '
         CST(1:NCST) = A(2:LA)
      ENDIF
*
      IF(CST(1:4).EQ.'read') IDC = 1
*
      RETURN
      END
*
      SUBROUTINE GUIPRI(INA,IUNIT,IWRT,NWRT)
*     **************************************
*
*     THIS SUBROUTINE OUTPUTS ARGUMENTS OF THE COSY WRITE COMMAND
*     TO THE GUI UNIT NUMBER IUNIT.  -201 >= IUNIT >= -210.
*     THE ARGUMENTS ARE SUBJECT TO INTERPRETATION
*     ACCORDING TO GUI WRITE RULES DESCRIBED IN THE MANUAL.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      CHARACTER CST*256,A*(2*4096),STR*4096,STRE*(2*4096+2),SNR*3
*
*     A:     THE COMMAND AND THE ARGUMENTS    A(1:LA) IS OUTPUT TO GUI
*     LA:    LOCATION OF THE END OF THE COMMAND WITH ARGUMENTS
*     IF LA>0 WHEN ENTERING INTO THIS ROUTINE, INA IS AN ARGUMENT
*     LCONS: CONSOLE FLAG    1: ON   0: OFF
*
      SAVE A,LA,LCONS
      DATA LA,LCONS / 0, 0 /
*
      CALL GCMD(INA,IDC,CST,NCST)
*
      IWND = -200-IUNIT
      WRITE(SNR,'(I3)') IWND
*
*   * GUI COMMAND
      IF(IDC.GE.0) THEN
*
*      * OUTPUT THE PREVIOUS GUI COMMAND WITH THE ARGUMENTS TO THE GUI.
         IF(LA.GT.0) CALL GUICMD(SNR//' '//A(1:LA))
         LCONS = 0
*
*      * NEW GUI COMMAND NAME
         A(1:NCST) = CST(1:NCST)
         LA = NCST
*
*      * CHECK IF IT RETURNS A VALUE AND STORE IN DELAYED READ ARRAY
         IF(IDC.GE.1) THEN
            NGVAR(IWND,1) = NGVAR(IWND,1)+1
            IGVAR(IWND,NGVAR(IWND,1),1) = -1
         ENDIF
*
*      * SPECIAL COMMANDS
*      * NOT A "REAL" GUI COMMAND, SET CONSOLE FLAG
         IF(CST.EQ.'console             ') THEN
            LCONS = 1
            LA = 0
*
*      * UPDATE THE DELAYED READ STRUCTURES TO MATCH NEW WINDOW
         ELSEIF(CST.EQ.'show                ') THEN
            N = NGVAR(IWND,1)
            NGVAR(IWND,2) = N
            NGVAR(IWND,1) = 0
            DO 10 I=1,N
            IGVAR(IWND,I,2) = IGVAR(IWND,I,1)
  10        CONTINUE
*
*      * UPDATE THE DELAYED READ STRUCTURES TO MATCH CLOSED WINDOW
         ELSEIF(CST.EQ.'close               ') THEN
            NGVAR(IWND,2) = 0
         ENDIF
*
*   * NEW ARGUMENT TO A PREVIOUS GUI COMMAND
*   * CONVERT TO TEXT, ESCAPE, AND APPEND TO GUI COMMAND
      ELSEIF(LA.GT.0) THEN
         CALL FOXSTR(INA,STR,LS)
         CALL GUIESC(STR,LS,STRE,LSE)
         A = A(1:LA)//' '//STRE(1:LSE)
         LA = LA+LSE+1
*
*   * OUTPUT TO THE GUI DIRECTLY
      ELSE
*
*      * STRING - CONVERT TO TEXT COMMAND IN GUI OUTPUT
*      *          UNLESS IN CONSOLE MODE DUE TO A \CONSOLE COMMAND
         IF(LCONS.EQ.0.AND.NTYP(INA).EQ.NST) THEN
            CALL STSTR(INA,STR,LS)
            CALL GUIESC(STR,LS,STRE,LSE)
            CALL GUICMD(SNR//' text '//STRE(1:LSE))
            CALL GUICMD(SNR//' newline')
*
*      * GRAPHICS - SET UP GUI FOR GRPRI AND STORE IN DELAYED READ ARRAY
         ELSEIF(LCONS.EQ.0.AND.NTYP(INA).EQ.NGR) THEN
            CALL GUICMD(SNR//' canvas')
            CALL GUICMD(SNR//' redirect')
            CALL GRPRI(INA,IUNIT)
            CALL GUICMD('0 redirect')
            CALL GUICMD(SNR//' newline')
            NGVAR(IWND,1) = NGVAR(IWND,1)+1
            IGVAR(IWND,NGVAR(IWND,1),1) = -1
*
*      * OTHERS - USE THE ASCII REPRESENTATION
         ELSE
            CALL FOXPRIO(INA,IUNIT)
         ENDIF
      ENDIF
*
*   * THE LAST ARGUMENT
*   * OUTPUT THE GUI COMMAND WITH THE ARGUMENTS TO THE GUI
*
      IF(IWRT.EQ.NWRT) THEN
         IF(LA.GT.0) CALL GUICMD(SNR//' '//A(1:LA))
         LCONS = 0
         LA = 0
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GREAD(IUNIT,INB)
*     ***************************
*
*     THIS SUBROUTINE WAITS FOR A BUTTON TO BE PUSHED, AND INITIATES
*     A DELAYED READ FOR ALL GUI VARIABLES
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      CHARACTER SNR*3,A*4096
*
      IF(IUNIT.LT.-210.OR.IUNIT.GT.-201) THEN
         PRINT*,'$$$ ERROR IN GREAD, ILLEGAL UNIT NUMBER ',IUNIT
         CALL FOXDEB
      ENDIF
*
* *** BOTH MODES - WAIT FOR THE BUTTON AND READ IT FROM STANDARD INPUT
      IWND = -200-IUNIT
      WRITE(SNR,'(I3)') IWND
      CALL GUICMD(SNR//' wait')
      CALL FOXREAO(INB,5)
      CALL GUICMD(SNR//' read')
*
* *** GUI MODE - READ FIELDS FROM INPUT UNIT NOW.
* *** IN ASCII GUI MODE THE READING OF FIELDS IS DONE IN ASCGUI DIRECTLY
* *** DURING read COMMAND PROCESSING
      IF(LGUI.EQ.1) THEN
         READ(5,*) NFLD
         DO 10 J=1,MIN(NFLD,NGVAR(IWND,2))
            IF(IGVAR(IWND,J,2).GT.IVAR) THEN
               PRINT*,'$$$ ERROR IN GREAD, '//
     *         'DELAYED READ VARIABLE NUMBER INVALID: ',IGVAR(IWND,J,2)
               CALL FOXDEB
            ELSEIF(IGVAR(IWND,J,2).GT.0) THEN
               CALL FOXREAO(IGVAR(IWND,J,2),5)
            ELSE
               READ(5,'(A)') A
            ENDIF
  10     CONTINUE
*
*        HANDLE ERRORS WHERE # OF LINES EXPECTED != # OF LINES SENT
*
         IF(NFLD.NE.NGVAR(IWND,2)) THEN
            PRINT*,'$$$ WARNING IN GREAD, GUI OUT OF SYNC '//
     *       '(RECEIVED,EXPECTED): ', NFLD, NGVAR(IWND,2)
            DO 20 J=NGVAR(IWND,2)+1, NFLD
               READ(5,'(A)') A
  20        CONTINUE
            DO 30 J=NFLD+1,NGVAR(IWND,2)
               IF(IGVAR(IWND,J,2).GT.IVAR) THEN
                  PRINT*,'$$$ ERROR IN GREAD, '//
     *          'DELAYED READ VARIABLE NUMBER INVALID: ',IGVAR(IWND,J,2)
                  CALL FOXDEB
               ELSEIF(IGVAR(IWND,J,2).GT.0) THEN
                  NTYP(IGVAR(IWND,J,2)) = NRE
                  CC(NBEG(IGVAR(IWND,J,2))) = 0.D0
               ENDIF
  30        CONTINUE
         ENDIF
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GSHOW(IUNIT)
*     ***********************
*
*     THIS SUBROUTINE ADDS AN OK BUTTON IF THERE IS NO OTHER BUTTON,
*     SHOWS THE NEW WINDOW, CALLS GREAD, AND DEACTIVATES THE WINDOW.
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      CHARACTER SNR*3
*
      IF(IUNIT.LT.-210.OR.IUNIT.GT.-201) THEN
         PRINT*,'$$$ ERROR IN GSHOW, ILLEGAL UNIT NUMBER ',IUNIT
         CALL FOXDEB
      ENDIF
*
      IWND = -200-IUNIT
      WRITE(SNR,'(I3)') IWND
      CALL GUICMD(SNR//' finish')
      CALL GUICMD(SNR//' show')
*
*     UPDATE THE DELAYED READ STRUCTURES
*
      N = NGVAR(IWND,1)
      NGVAR(IWND,2) = N
      NGVAR(IWND,1) = 0
      DO 10 I=1,N
      IGVAR(IWND,I,2) = IGVAR(IWND,I,1)
10    CONTINUE
*
      CALL FOXALL(ISC,1,1)
      CALL GREAD(IUNIT,ISC)
      CALL FOXDAL(ISC,1)
      CALL GUICMD(SNR//' deactivate')
*
      RETURN
      END
*
      SUBROUTINE GUICMD(CMD)
*     **********************
*
*     THIS SUBROUTINE WRITES A GUI COMMAND
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      CHARACTER CMD*(*),CSCR*(2*4096+35)
*
      IF(LGUI.EQ.1) THEN
C        WRITE(6,'(A)') '<*%GUI%*>'//CMD//'</*%GUI%*>'
C        DUE TO A G77 PROBLEM, THE ABOVE LINE IS CHANGED TO THE NEXT TWO LINES.
         CSCR = '<*%GUI%*>'//CMD//'</*%GUI%*>'
         WRITE(6,'(A)') CSCR(1:(LEN(CMD)+19))
      ELSE
         CALL ASCGUI(CMD)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE ASCGUI(A)
*     ********************
*
*     THIS SUBROUTINE EMULATES A SIMPLE ASCII BASED GUI PROGRAM.
*     THE ARGUMENT ARE GUI COMMAND AS WRITTEN TO OUTSIDE GUIS.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      PARAMETER(LWORD=20)
      CHARACTER A*(*),AA*(2*4096+16),CST*256,CHA
      INTEGER IWORD(LWORD,2),IBUTTON(LGW)
*
*     HAS THERE BEEN A BUTTON ADDED TO THIS WINDOW YET
      DATA IBUTTON / LGW*0 /
      SAVE IBUTTON
*
*     INTERPRET GUI COMMAND STRING
      CALL GUISTR(A,LEN(A),AA,IWND,CST,NCST,NWORD,IWORD)
*
*     IGNORE COMMANDS SENT TO ANYTHING BUT A GUI WINDOW (IWND = 1...10)
      IF(IWND.LT.1.OR.IWND.GT.10) RETURN
*
*     HANDLE THE MOST ESSENTIAL COMMANDS FOR THE ASCII GUI
      IF    (CST.EQ.'text                ') THEN
         IF(NWORD.GT.0) THEN
            WRITE(6,'(A)') AA(IWORD(1,1):IWORD(1,2))
         ELSE
            WRITE (6,'(A)') ''
         ENDIF
      ELSEIF(CST.EQ.'line                ') THEN
            WRITE (6,'(A)') '******************************************'
      ELSEIF(CST.EQ.'button              ') THEN
         IF(NWORD.EQ.0) THEN
*           THERE IS 1 REQUIRED ARGUMENT, THIS IS AN ERROR
            WRITE(6,'(A)') 'NOT ENOUGH ARGUMENTS FOR BUTTON'
         ELSE
            WRITE(6,'(A)') '[['//AA(IWORD(1,1):IWORD(1,2))//']]'
            IBUTTON(IWND) = 1
         ENDIF
      ELSEIF(CST.EQ.'finish              ') THEN
         IF(IBUTTON(IWND).EQ.0) THEN
            IF(NWORD.EQ.0) THEN
               WRITE(6,'(A)') '[[OK]]'
            ELSE
               WRITE(6,'(A)') '[['//AA(IWORD(1,1):IWORD(1,2))//']]'
            ENDIF
            IBUTTON(IWND) = 1
         ENDIF
      ELSEIF(CST.EQ.'readcheckbox        ') THEN
         IF(NWORD.EQ.0) THEN
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),' [ ]'
         ELSEIF(NWORD.EQ.1) THEN
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),
     *         ' [ ] '//AA(IWORD(1,1):IWORD(1,2))
         ELSE
            I = IBEGIN(AA,IWORD(2,1),IWORD(2,2))
            IF((I.LE.IWORD(2,2)).AND.
     *         ((AA(I:I).EQ.'1').OR.(AA(I:I).EQ.'T'))) THEN
               CHA = 'X'
            ELSE
               CHA = ' '
            ENDIF
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),
     *         ' ['//CHA//'] '//AA(IWORD(1,1):IWORD(1,2))
         ENDIF
      ELSEIF(CST.EQ.'readoption          ') THEN
         IF(NWORD.EQ.0) THEN
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),' ( )'
         ELSEIF(NWORD.EQ.1) THEN
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),
     *         ' ( ) '//AA(IWORD(1,1):IWORD(1,2))
         ELSE
            I = IBEGIN(AA,IWORD(2,1),IWORD(2,2))
            IF((I.LE.IWORD(2,2)).AND.
     *         ((AA(I:I).EQ.'1').OR.(AA(I:I).EQ.'T'))) THEN
               CHA = '*'
            ELSE
               CHA = ' '
            ENDIF
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),
     *         ' ('//CHA//') '//AA(IWORD(1,1):IWORD(1,2))
         ENDIF
      ELSEIF(CST.EQ.'readfield           '.OR.
     *       CST.EQ.'readfilename        ') THEN
         IF(NWORD.EQ.0) THEN
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1), ' | |'
         ELSE
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),
     *         ' |'//AA(IWORD(1,1):IWORD(1,2))//'|'
         ENDIF
      ELSEIF(CST.EQ.'readnumber          ') THEN
         IF(NWORD.LT.3) THEN
*           THERE ARE 3 REQUIRED ARGUMENTS, THIS IS AN ERROR
            WRITE(6,'(A)') 'NOT ENOUGH ARGUMENTS FOR READNUMBER'
         ELSE
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),
     *         ' #'//AA(IWORD(1,1):IWORD(1,2))//'#'//
     *         ' ['//AA(IWORD(2,1):IWORD(2,2))//', '//
     *         ' '//AA(IWORD(3,1):IWORD(3,2))//']'
         ENDIF
      ELSEIF(CST.EQ.'readlist            ') THEN
         IF(NWORD.LT.2) THEN
*           THERE ARE 2 REQUIRED ARGUMENTS, THIS IS AN ERROR
            WRITE(6,'(A)') 'NOT ENOUGH ARGUMENTS FOR READLIST'
         ELSEIF(NWORD.EQ.2) THEN
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),
     *         ' !'//AA(IWORD(2,1):IWORD(2,2))//'! {'//
     *         AA(IWORD(1,1):IWORD(1,2))//'}'
         ELSE
            I = IBEGIN(AA,IWORD(3,1),IWORD(3,2))
            IF((I.LE.IWORD(3,2)).AND.
     *         ((AA(I:I).EQ.'1').OR.(AA(I:I).EQ.'T'))) THEN
               CHA = '|'
            ELSE
               CHA = '/'
            ENDIF
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),
     *         ' '//CHA//AA(IWORD(2,1):IWORD(2,2))//CHA//' {'//
     *         AA(IWORD(1,1):IWORD(1,2))//'}'
         ENDIF
      ELSEIF(CST.EQ.'readprogress        ') THEN
         IF(NWORD.EQ.0) THEN
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1), ' >>>>'
         ELSE
            WRITE(6,'(I2,A,I0,A)') IWND, '.', NGVAR(IWND,1),
     *         ' >>>> '//AA(IWORD(1,1):IWORD(1,2))//'%'
         ENDIF
      ELSEIF(CST.EQ.'show                ') THEN
*        THE NEW HIDDEN WINDOW DOES NOT HAVE A BUTTON YET
         IBUTTON(IWND) = 0
      ELSEIF(CST.EQ.'wait                ') THEN
*        PROMPT USER TO ENTER A BUTTON NAME, GREAD WILL READ IT FOR US
         WRITE(6,'(/,A)') 'Please enter the text on the button to push:'
      ELSEIF(CST.EQ.'read                ') THEN
*        PROMPT USER TO ENTER VALUES FOR EACH FIELD
*        VALUES ARE READ DIRECTLY INTO COSY VARIABLES HERE
         DO 10 J=1,NGVAR(IWND,2)
            IF(IGVAR(IWND,J,2).GT.IVAR) THEN
               PRINT*,'$$$ ERROR IN ASCGUI, '//
     *         'DELAYED READ VARIABLE NUMBER INVALID: ',IGVAR(IWND,J,2)
               CALL FOXDEB
            ELSEIF(IGVAR(IWND,J,2).GT.0) THEN
               WRITE(6,'(A,I2,A,I0)') 'Please enter value for field ',
     *            IWND, '.', J
               CALL FOXREAO(IGVAR(IWND,J,2),5)
            ENDIF
  10     CONTINUE
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GUISTR(A,LA,AA,IWND,CST,NCST,NWORD,IWORD)
*     ****************************************************
*
*     THIS SUBROUTINE ANALYZES A GUI OUTPUT STRING A(1:LA).
*     SEE GUI OUTPUT RELATED ROUTINES, ESPECIALLY GUIPRI, GUIESC,
*     FOR THE CONSTRUCTION RULES ON GUI OUTPUT STRINGS.
*
*     A GUI OUTPUT STRING CONSISTS OF THE FOLLOWING COMPONENTS.
*     1. GUI WINDOW NUMBER              IWND
*     2. COSY-GUI COMMAND               CST(1:NCST)
*     3. ARGUMENTS TO THE COMMNAD
*        NUMBER OF ARGUMENTS            NWORD
*        STRING AS THE I-TH ARGUMENT    AA(IWORD(I,1):IWORD(I,2))
*
*     NOTE: A(1:LA) IS ASSUMED TO HAVE AT LEAST TWO WORDS,
*     NAMELY THE GUI WINDOW NUMBER AND THE COSY-GUI COMMAND.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----CHARACTER HANDLING ----------------------------------------------------
      INTEGER ICSA,ICSZ,IDIF,ICCA,ICCZ
      COMMON /TRCHAR/ ICSA,ICSZ,IDIF,ICCA,ICCZ
      CHARACTER SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
      COMMON /SCHAR/ SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
*     SCNUL={NUL}, SCHT={HORIZONTAL TAB}, SCLF={LINE FEED},
*     SCCR={CARRIAGE RETURN}, SCSP={SPACE}=' ',
*     SCQM={QUOTATION MARK}='"', SCBSL={BACK SLASH}='\'
*----------------------------------------------------------------------------
*
      PARAMETER(LWORD=20)
      CHARACTER A*(*),AA*(*),CST*256,ESC*2,AAA,CSCR*(2*4096+17)
      DIMENSION IWORD(LWORD,2)
*
      SAVE ICALL,ESC
      DATA ICALL / 0 /
*
      IF(ICALL.EQ.0) THEN
         ICALL = 1
         ESC = SCQM//SCBSL
      ENDIF
*
*     IA:  THE STARTING POSITION IN A TO BE ANALYZED.
*          AS THE TASK PROGRESSES, IA IS UPDATED.
*     IS:  THE POSITION COUNTER TO LOCATE A DELIMITTER IN A(IA:LA)
*          THE IS COUNTING STARTS AT THE POSITION IA, NAMELY AT IA, IS=1.
*          THE DELIMITTER FOR THE BEGINNING TWO WORDS IS ' ' {SPACE, BLANK}.
*
*   * TAKE THE FIRST WORD, WHICH IS A NUMBER FOR IWND
      IA = IBEGIN(A,1,LA)
      IS = INDEX(A(IA:LA),SCSP)
      CALL VALCH(A,IA,IA+IS-2,VAL,IER)
      IF(IER.NE.0) THEN
         PRINT*,'$$$ ERROR IN GUISTR, 1ST WORD IS NOT A NUMBER'
         CALL FOXDEB
      ENDIF
      IWND = NINT(VAL)
*
      NWORD = 0
*
*   * TAKE THE SECOND WORD, WHICH IS A STRING FOR THE COMMAND
      IA = IBEGIN(A,IA+IS,LA)
C     IS = INDEX(A(IA:LA)//SCSP,SCSP)
C     DUE TO A G77 PROBLEM, THE ABOVE LINE IS CHANGED TO THE NEXT TWO LINES.
      CSCR = A(IA:LA)//SCSP
      IS = INDEX(CSCR,SCSP)
      NCST = IS-1
*
      IF(NCST.GT.256) THEN
         PRINT*,'$$$ ERROR IN GUISTR, 2ND WORD IS TOO LONG'
         CALL FOXDEB
      ENDIF
      CST = '                    '
      CST(1:NCST) = A(IA:(IA+NCST-1))
*
      IA = IA+IS
      IF(IA.GT.LA) RETURN
*
*   * PROCESS THE REST OF THE STRING, CONSISTING OF ARGUMENTS
      IA = IBEGIN(A,IA,LA)
      A(LA+1:LA+1) = SCSP
*
*     IP:  THE POINTER IN AA
*     IW:  THE COUNTER OF CHARACTERS IN ONE ARGUMENT
*     IQM: FLAG TO INDICATE THAT IT IS INSIDE AN ARGUMENT OPENED AND
*          CLOSED BY '"' {QUOTATION MARK}.   1: INSIDE   0 : OUTSIDE
*     IES: FLAG TO INDICATE A GUI ESCAPING BY '\' {BACK SLACH}
*          PRECEDING '\' OR '"'.             1: GUI ESCAPING   0: NON
*
      IQM = 0
      IES = 0
      IP = 0
*
      DO 10 I=IA,LA
      AAA = A(I:I)
      IF(IQM.EQ.0.AND.AAA.EQ.SCQM) THEN
         IQM = 1
         NWORD = NWORD+1
         IW = 0
         IWORD(NWORD,1) = IP+1
      ELSEIF(IQM.EQ.1) THEN
         IF(AAA.EQ.SCBSL) THEN
            IF(IES.EQ.0.AND.INDEX(ESC,A(I+1:I+1)).NE.0) THEN
               IES = 1
               GOTO 10
            ELSE
               IES = 0
            ENDIF
         ELSEIF(AAA.EQ.SCQM) THEN
            IF(IES.EQ.1) THEN
               IES = 0
            ELSE
               IQM = 0
               IWORD(NWORD,2) = IP
               GOTO 10
            ENDIF
         ENDIF
         IP = IP+1
         AA(IP:IP) = AAA
         IW = IW+1
         IF(IW.EQ.1) IWORD(NWORD,1) = IP
      ENDIF
  10  CONTINUE
*
      IF(NWORD.GE.0.AND.IQM.EQ.1) IWORD(NWORD,2) = IP
*
      RETURN
      END
*
      SUBROUTINE GUISET(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE SETS GUI FIELD INB IN GUI WINDOW UNIT INA TO THE VALUE
*     IN INC. ALSO WORKS FOR GRAPHS.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      CHARACTER STR*4096,STRE*(2*4096+2),SNR*4,GNR*3
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NINT(CC(NBEG(INA)))
      IF((IA.LT.-210).OR.(IA.GT.-201)) THEN
         PRINT*,'$$$ ERROR IN GUISET, ILLEGAL GUI UNIT: ',IA
         CALL FOXDEB
      ENDIF
*
      IB = NINT(CC(NBEG(INB)))
      IF(IB.LT.1) THEN
         PRINT*,'$$$ ERROR IN GUISET, ILLEGAL GUI ELEMENT NUMBER: ',IB
         CALL FOXDEB
      ENDIF
*
      WRITE(SNR,'(I3)') -200-IA
      WRITE(GNR,'(I3)') IB
      IF(NTYP(INC).EQ.NGR) THEN
         CALL GUICMD(SNR//' set '//GNR//' 0')
         WRITE(SNR,'(I4)') 200+IA
         CALL GUICMD(SNR//' redirect')
         CALL GRPRI( INC, IA )
         CALL GUICMD('0 redirect')
      ELSE
         CALL FOXSTR(INC,STR,LS)
         CALL GUIESC(STR,LS,STRE,LSE)
         CALL GUICMD(SNR//' set '//GNR//' '//STRE(1:LSE))
      ENDIF
*
      RETURN
      END
*
*%%
      INTEGER FUNCTION ILAST(A,IA1,IA2)
*     *********************************
*
*     THIS FUNCTION DETERMINES THE LAST NON SPACE CHARACTER POSITION
*     IN CHARACTER A BETWEEN POSITIONS IA1 AND IA2
*
*-----CHARACTER HANDLING ----------------------------------------------------
      INTEGER ICSA,ICSZ,IDIF,ICCA,ICCZ
      COMMON /TRCHAR/ ICSA,ICSZ,IDIF,ICCA,ICCZ
      CHARACTER SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
      COMMON /SCHAR/ SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
*     SCNUL={NUL}, SCHT={HORIZONTAL TAB}, SCLF={LINE FEED},
*     SCCR={CARRIAGE RETURN}, SCSP={SPACE}=' ',
*     SCQM={QUOTATION MARK}='"', SCBSL={BACK SLASH}='\'
*----------------------------------------------------------------------------
*
      CHARACTER A*(*)
*
      DO 10 I=IA2,IA1,-1
      IF(A(I:I).NE.SCNUL.AND.A(I:I).NE.SCSP.AND.A(I:I).NE.SCHT) THEN
         ILAST = I
         RETURN
      ENDIF
  10  CONTINUE
*
      ILAST = IA1-1
      RETURN
*
      END
*
      INTEGER FUNCTION IBEGIN(A,IA1,IA2)
*     **********************************
*
*     THIS FUNCTION DETERMINES THE FIRST NON SPACE CHARACTER POSITION
*     IN CHARACTER A BETWEEN POSITIONS IA1 AND IA2
*
*-----CHARACTER HANDLING ----------------------------------------------------
      INTEGER ICSA,ICSZ,IDIF,ICCA,ICCZ
      COMMON /TRCHAR/ ICSA,ICSZ,IDIF,ICCA,ICCZ
      CHARACTER SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
      COMMON /SCHAR/ SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
*     SCNUL={NUL}, SCHT={HORIZONTAL TAB}, SCLF={LINE FEED},
*     SCCR={CARRIAGE RETURN}, SCSP={SPACE}=' ',
*     SCQM={QUOTATION MARK}='"', SCBSL={BACK SLASH}='\'
*----------------------------------------------------------------------------
*
      CHARACTER A*(*)
*
      DO 10 I=IA1,IA2
      IF(A(I:I).NE.SCNUL.AND.A(I:I).NE.SCSP.AND.A(I:I).NE.SCHT) THEN
         IBEGIN = I
         RETURN
      ENDIF
  10  CONTINUE
*
      IBEGIN = IA2+1
      RETURN
*
      END
*
      SUBROUTINE CAPSTR(A,IA1,IA2)
*     ****************************
*
*     THIS SUBROUTINE CAPITALIZES THE LETTERS IN CHARACTER A(IA1:IA2)
*
*-----CHARACTER HANDLING ----------------------------------------------------
      INTEGER ICSA,ICSZ,IDIF,ICCA,ICCZ
      COMMON /TRCHAR/ ICSA,ICSZ,IDIF,ICCA,ICCZ
      CHARACTER SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
      COMMON /SCHAR/ SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
*     SCNUL={NUL}, SCHT={HORIZONTAL TAB}, SCLF={LINE FEED},
*     SCCR={CARRIAGE RETURN}, SCSP={SPACE}=' ',
*     SCQM={QUOTATION MARK}='"', SCBSL={BACK SLASH}='\'
*----------------------------------------------------------------------------
*
      CHARACTER A*(*)
*
      DO 10 I=IA1,IA2
      ICI = ICHAR(A(I:I))
      IF(ICI.GE.ICSA.AND.ICI.LE.ICSZ) A(I:I) = CHAR(ICI+IDIF)
  10  CONTINUE
*
      RETURN
      END
*
      SUBROUTINE GUIESC(AI,LAI,AG,LAG)
*     ********************************
*
*     THIS SUBROUTINE CONVERTS A STRING AI WITH LENGTH LAI
*     TO A GUI ESCAPED STRING AG WITH LENTH LAG.   1 =< LAG =< 2*LAI+2
*     WHEN CALLING GUIESC, MAKE SURE AI HAS THE LENGTH AT LEAST LAI.
*
*-----CHARACTER HANDLING ----------------------------------------------------
      INTEGER ICSA,ICSZ,IDIF,ICCA,ICCZ
      COMMON /TRCHAR/ ICSA,ICSZ,IDIF,ICCA,ICCZ
      CHARACTER SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
      COMMON /SCHAR/ SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
*     SCNUL={NUL}, SCHT={HORIZONTAL TAB}, SCLF={LINE FEED},
*     SCCR={CARRIAGE RETURN}, SCSP={SPACE}=' ',
*     SCQM={QUOTATION MARK}='"', SCBSL={BACK SLASH}='\'
*----------------------------------------------------------------------------
*
      CHARACTER AI*(*),AG*(*),ESC*2
      SAVE ICALL,ESC
      DATA ICALL / 0 /
*
      IF(ICALL.EQ.0) THEN
         ICALL = 1
         ESC = SCQM//SCBSL
      ENDIF
*
      LAG = 1
      AG(LAG:LAG) = SCQM
*
      DO 10 I=1,LAI
      IF(INDEX(ESC,AI(I:I)).NE.0) THEN
         LAG = LAG+1
         AG(LAG:LAG) = SCBSL
      ENDIF
      LAG = LAG+1
      AG(LAG:LAG) = AI(I:I)
  10  CONTINUE
*
      LAG = LAG+1
      AG(LAG:LAG) = SCQM
*
      RETURN
      END
*
      SUBROUTINE STSTR(INA,A,LA)
*     **************************
*
*     THIS SUBROUTINE STORES THE STRING INA WITH LENGTH LA IN A(1:LA).
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER A*(*)
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      LA = NEND(INA)-IA+1
      DO 10 I=1,LA
      A(I:I) = CHAR(NC(IA+I-1))
  10  CONTINUE
*
      RETURN
      END
*
      SUBROUTINE FOXSTR(INA,A,LA)
*     ***************************
*
*     THIS SUBROUTINE CONVERTS INA TO A STRING A WITH THE COSY DEAULT FORMAT.
*     THE RESULTING LENGTH LA OF A IS RETURNED.
*     FOR ST TYPE:      LA IS THE LENGTH OF INA.
*                       DECLARE THE INCOMING CHARCTER A WITH SIZE 4096.
*     FOR LO TYPE:      LA = 5.
*     FOR NON ST TYPES: LA EXCLUDES BLANKS AT THE END.
*                       DECLARE THE INCOMING CHARCTER A WITH SIZE AT LEAST 50.
*     REFER TO RECSTC FOR THE COSY DEFAULT FORMAT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER A*(*)
*
      NTY = NTYP(INA)
      IA = NBEG(INA)
      IF(NTY.EQ.NRE) THEN
         A(1:25) = '                         '
         WRITE(A,'(G23.16E3)') CC(IA)
         LA = ILAST(A,1,23)
      ELSEIF(NTY.EQ.NST) THEN
         LA = NEND(INA)-IA+1
         DO 10 I=1,LA
         A(I:I) = CHAR(NC(IA+I-1))
  10     CONTINUE
      ELSEIF(NTY.EQ.NLO) THEN
         IF(NC(IA).EQ.0) THEN
            A(1:5) = 'FALSE'
         ELSE
            A(1:5) = 'TRUE '
         ENDIF
         LA = 5
      ELSEIF(NTY.EQ.NCM) THEN
         A(1:40) = '                                        '
         CALL CMCSTO(CC(IA),CC(IA+1),'(G17.9E3)           ',A)
         LA = 37
      ELSE
         PRINT*,'$$$ ERROR IN FOXSTR, ILLEGAL TYPE ',NTY
         CALL FOXDEB
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE STRIM(INA,INC)
*     *************************
*
*     THIS SUBROUTINE REMOVES SPACE CHARACTERS FROM THE END OF
*     A STRING INA AND STORES IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----CHARACTER HANDLING ----------------------------------------------------
      INTEGER ICSA,ICSZ,IDIF,ICCA,ICCZ
      COMMON /TRCHAR/ ICSA,ICSZ,IDIF,ICCA,ICCZ
      CHARACTER SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
      COMMON /SCHAR/ SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
*     SCNUL={NUL}, SCHT={HORIZONTAL TAB}, SCLF={LINE FEED},
*     SCCR={CARRIAGE RETURN}, SCSP={SPACE}=' ',
*     SCQM={QUOTATION MARK}='"', SCBSL={BACK SLASH}='\'
*----------------------------------------------------------------------------
*
      CHARACTER AA
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      LA = NEND(INA)-IA+1
      IC = NBEG(INC)-1
*
      DO 10 I=LA,1,-1
      AA = CHAR(NC(IA+I-1))
      IF(AA.NE.SCNUL.AND.AA.NE.SCSP.AND.AA.NE.SCHT) THEN
         IAL = I
         GOTO 20
      ENDIF
  10  CONTINUE
      GOTO 40
*
  20  CONTINUE
      DO 30 I=1,IAL
      IC = IC+1
      NC(IC) = NC(IA+I-1)
  30  CONTINUE
*
  40  CONTINUE
      NEND(INC) = IC
      NTYP(INC) = NST
      IF(IC.GT.NMAX(INC)) CALL FOXERV('STRIM')
*
      RETURN
      END
*
      SUBROUTINE SLTRIM(INA,INC)
*     **************************
*
*     THIS SUBROUTINE REMOVES SPACE CHARACTERS FROM THE BEGINNING OF
*     A STRING INA AND STORES IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----CHARACTER HANDLING ----------------------------------------------------
      INTEGER ICSA,ICSZ,IDIF,ICCA,ICCZ
      COMMON /TRCHAR/ ICSA,ICSZ,IDIF,ICCA,ICCZ
      CHARACTER SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
      COMMON /SCHAR/ SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
*     SCNUL={NUL}, SCHT={HORIZONTAL TAB}, SCLF={LINE FEED},
*     SCCR={CARRIAGE RETURN}, SCSP={SPACE}=' ',
*     SCQM={QUOTATION MARK}='"', SCBSL={BACK SLASH}='\'
*----------------------------------------------------------------------------
*
      CHARACTER AA
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      LA = NEND(INA)-IA+1
      IC = NBEG(INC)-1
*
      DO 10 I=1,LA
      AA = CHAR(NC(IA+I-1))
      IF(AA.NE.SCNUL.AND.AA.NE.SCSP.AND.AA.NE.SCHT) THEN
         IAS = I
         GOTO 20
      ENDIF
  10  CONTINUE
      GOTO 40
*
  20  CONTINUE
      DO 30 I=IAS,LA
      IC = IC+1
      NC(IC) = NC(IA+I-1)
  30  CONTINUE
*
  40  CONTINUE
      NEND(INC) = IC
      NTYP(INC) = NST
      IF(IC.GT.NMAX(INC)) CALL FOXERV('SLTRIM')
*
      RETURN
      END
*
      SUBROUTINE DFCHAR
*     *****************
*
*     THIS SUBROUTINE DEFINES THE CONTENTS OF COMMON TRCHAR AND SCHAR
*
*-----CHARACTER HANDLING ----------------------------------------------------
      INTEGER ICSA,ICSZ,IDIF,ICCA,ICCZ
      COMMON /TRCHAR/ ICSA,ICSZ,IDIF,ICCA,ICCZ
      CHARACTER SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
      COMMON /SCHAR/ SCNUL,SCHT,SCLF,SCCR,SCSP,SCQM,SCBSL
*     SCNUL={NUL}, SCHT={HORIZONTAL TAB}, SCLF={LINE FEED},
*     SCCR={CARRIAGE RETURN}, SCSP={SPACE}=' ',
*     SCQM={QUOTATION MARK}='"', SCBSL={BACK SLASH}='\'
*----------------------------------------------------------------------------
*
      ICSA = ICHAR('a')
      ICSZ = ICHAR('z')
      ICCA = ICHAR('A')
      ICCZ = ICHAR('Z')
      IDIF = ICCA-ICSA
C
C     IF ACHAR IS NOT AVAILABLE (SUCH AS PGI F77 PGF77), USE CHAR INSTEAD.
C     MAKE SURE THE RESULTING CHARACTERS SCBSL ETC. ARE AS INTENDED ABOVE.
C     NOTE: PGF95 SUPPORTS ACHAR.
C
      SCNUL = ACHAR(0)
      SCHT  = ACHAR(9)
      SCLF  = ACHAR(10)
      SCCR  = ACHAR(13)
      SCSP  = ACHAR(32)
      SCQM  = ACHAR(34)
      SCBSL = ACHAR(92)
*
      RETURN
      END
*
*%%
      SUBROUTINE MEMDPV(INUNIT,INA)
*     *****************************
*
*     THIS SUBROUTINE INTERFACES TO PRINT DUMP OF INA TO UNIT IUNIT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INUNIT).NE.NRE) CALL FOXNTY(INUNIT)
*
      IUNIT = NINT(CC(NBEG(INUNIT)))
      CALL FOXDMP(IUNIT,INA)
*
      RETURN
      END
*
      SUBROUTINE FOXDMP(IUNIT,INA)
*     ****************************
*
*     THIS SUBROUTINE PRINTS DUMP OF INA TO UNIT IUNIT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      WRITE(IUNIT,'(A,I8,A,I3,A,I9)') ' ** FOXDMP ** DUMP OF INA (=',
     * INA,'),    NTYP() =',NTYP(INA),',  NBEG() =',NBEG(INA)
      WRITE(IUNIT,'(A,I3,A,I5,A,I6,A,15X,A,I9)') ' NOMAX =',NOMAX,
     * ', NVMAX =',NVMAX,', NMMAX =',NMMAX,',','NEND() =',NEND(INA)
      WRITE(IUNIT,'(A,I3,A,D13.6,A,22X,A,I9)') ' NOCUT =',NOCUT,
     * ', EPS   =',EPS,',','NMAX() =',NMAX(INA)
*
      WRITE(IUNIT,'(A)') ' ADDRESS I     NC(I)          CC(I)'
      DO 100 I=NBEG(INA),NEND(INA)
      WRITE(IUNIT,'(2I10,3X,G23.16E3)') I,NC(I),CC(I)
 100  CONTINUE
*
      RETURN
      END
*
      SUBROUTINE MEMWRT(IUNIT)
*     ************************
*
*     THIS SUBROUTINE WRITES MEMORY TO UNIT
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      LOGICAL FL
      FL = .FALSE.
*
      IF(NTYP(IUNIT).NE.NRE) CALL FOXNTY(IUNIT)
      IU = NINT(CC(NBEG(IUNIT)))
*
      DO 10 I=1,IVAR
        WRITE(IU,'(1X,5I7,G27.16E3,I5)')
     *  I,NBEG(I),NEND(I),NMAX(I),NTYP(I),CC(NBEG(I)),NC(NBEG(I))
        DO 20 J=NBEG(I)+1,NMAX(I)
          IF(CC(J).EQ.0.D0.AND.NC(J).EQ.0) THEN
            IF(FL) THEN
              NUM = NUM+1
            ELSE
              NUM = 1
              FL = .TRUE.
            ENDIF
          ELSE
            IF(FL) THEN
              WRITE(IU,'(8X,I7,A53)') NUM,
     *'  TIMES                   **********************   **'
              FL = .FALSE.
            ENDIF
            WRITE(IU,'(8X,I7,G48.16E3,I5)') J,CC(J),NC(J)
          ENDIF
   20   CONTINUE
        IF(FL) THEN
          WRITE(IU,'(8X,I7,A53)') NUM,
     *'  TIMES                   **********************   **'
          FL = .FALSE.
        ENDIF
   10   CONTINUE
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     DIFFERENTIAL ALGEBRA ROUTINES                                           *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE DAINI(INO,INV,IIU,INM)
*     *********************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----ROUNDING --------------------------------------------------------------
      INTEGER INTINI
      DOUBLE PRECISION RINUP(2),RINDN(2),FINSR
      COMMON /INCOM/ RINUP,RINDN,FINSR,INTINI
*----------------------------------------------------------------------------
*
      COMMON /MONDA/ ID1(LNV),ID2(LNV)
      COMMON /TREEF/ NFATH(LEA)
      COMMON /INBCOM/ IJB(LEA)
      DIMENSION N(LNV+1),J(LNV),JJ(LNV),IEW0(LNV+1)
*
      IF(INTINI.EQ.-1) CALL INRNDO
*
      IF(NTYP(INO).NE.NRE) CALL FOXNTY(INO)
      IF(NTYP(INV).NE.NRE) CALL FOXNTY(INV)
      IF(NTYP(IIU).NE.NRE) CALL FOXNTY(IIU)
*
      CNO = CC(NBEG(INO))
      CNV = CC(NBEG(INV))
      NO = NINT(CNO)
      NV = NINT(CNV)
      IU = NINT(CC(NBEG(IIU)))
*
      IF(CNO.LT.1.D0.OR.CNV.LT.1.D0) THEN
         PRINT*, '$$$ ERROR IN DAINI, NO OR NV < 1'
         CALL FOXSTL
      ELSEIF(CNV.GT.DBLE(LNV).OR.CNO.GT.DBLE(LNO)) THEN
         PRINT*, '!!! MEMORY EXHAUSTION, INCREASE LNO, LNV BEYOND ',
     *        NO,NV
         CALL FOXSTL
      ENDIF
*
      NOMAX = NO
      NOCUT = NO
      NVMAX = NV
      CALL DANUM(NO,NV,NMMAX)
      CALL DANUMD(NO,NV,CNMMAX)
*
      IF(IU.NE.0) THEN
         DO 10 I=0,LIA
         IA1(I) = 0
         IA2(I) = 0
  10     CONTINUE
      ENDIF
*
      DO 20 I=1,LNV+1
      IEW0(I) = 1
  20  CONTINUE
*
      IF(LEWI.EQ.1) THEN
         LEWI = 0
         LEW  = 1
         IEWM = 1
         DO 30 I=1,NV
         IF(IEW(I).GT.NO) THEN
            PRINT*,'$$$ WARNING IN DAINI, THE WEIGHT OF VARIABLE',
     *             I,' EXCEEDS ORDER NO'
         ENDIF
         IEWM = IEWM*IEW(I)
         IEW0(NV+1-I) = IEW(I)
  30     CONTINUE
         IF(IEWM.EQ.1) LEW = 0
      ELSE
         DO 40 I=1,NV
         IEW(I) = 1
  40     CONTINUE
      ENDIF
*
      IBASE = NO+1
      IF(LEW.EQ.0) THEN
         RIA = FLOAT(IBASE)**((NV+1)/2)
         IF(RIA.GT.DBLE(LIA)) THEN
            IF(RIA.GT.2147483647.D0) THEN
               PRINT*,'!!! MEMORY EXHAUSTION, IT WILL CAUSE LIA LIMIT.'
            ELSE
               PRINT*,'!!! MEMORY EXHAUSTION, INCREASE LIA BEYOND ',
     *                 NINT(RIA)
            ENDIF
            CALL FOXSTL
         ELSEIF(CNMMAX.GT.DBLE(LEA)) THEN
            IF(CNMMAX.GT.2147483647.D0) THEN
               PRINT*,'!!! MEMORY EXHAUSTION, IT WILL CAUSE LEA LIMIT.'
            ELSE
               PRINT*,'!!! MEMORY EXHAUSTION, INCREASE LEA BEYOND ',
     *                 NMMAX
            ENDIF
            CALL FOXSTL
         ENDIF
      ENDIF
*
      IESP = (NV+1)/2
      JS = NV/2
      IF(JS.EQ.0) THEN
         NO2 = 0
      ELSE
         NO2 = NO
      ENDIF
*
      ICMAX = 0
      NN = 0
*
      DO 100 IO2=0,NO2
*
      N(1) = IO2
      JL = 0
      JD = 1
*
  50  JL = JL+JD
*
      IF(JL.LE.0) THEN
         GOTO 100
      ELSEIF(JD.EQ.1) THEN
         J(JL) = 0
      ELSE
         J(JL) = J(JL)+IEW0(JL)
      ENDIF
*
      N(JL+1) = N(JL)-J(JL)
*
      IF(J(JL).GT.N(JL)) THEN
         JD = -1
         GOTO 50
      ELSEIF(JL.LT.JS) THEN
         JD = 1
         GOTO 50
      ELSE
         J(JL) = (N(JL)/IEW0(JL))*IEW0(JL)
*
         IF(J(JL).EQ.N(JL)) THEN
            DO 60 IJ2=1,NV
            JJ(IJ2) = J(NV+1-IJ2)
  60        CONTINUE
            CALL DADEC(JJ,ICC1,ICC2)
            IC2 = ICC2
            ICMAX = MAX(ICMAX,IC2)
            IF(LEW.EQ.1) THEN
               IF(ICMAX.GT.LIA) THEN
                  PRINT*,'!!! MEMORY EXHAUSTION, INCREASE LIA'
                  CALL FOXSTL
               ENDIF
            ENDIF
*
            IA2(IC2) = NN
*
            DO 90 IO1=0,NO-IO2
*
            N(JS+1) = IO1
            JD = 1
*
  70        JL = JL+JD
*
            IF(JL.EQ.JS) THEN
               GOTO 90
            ELSEIF(JD.EQ.1) THEN
               J(JL) = 0
            ELSE
               J(JL) = J(JL)+IEW0(JL)
            ENDIF
*
            N(JL+1) = N(JL)-J(JL)
*
            IF(J(JL).GT.N(JL)) THEN
               JD = -1
               GOTO 70
            ELSEIF(JL.LT.NV) THEN
               JD = 1
               GOTO 70
            ELSE
               JD = -1
               J(JL) = (N(JL)/IEW0(JL))*IEW0(JL)
*
               IF(J(JL).EQ.N(JL)) THEN
                  DO 80 IJ2=1,NV
                  JJ(IJ2) = J(NV+1-IJ2)
  80              CONTINUE
                  CALL DADEC(JJ,ICC1,ICC2)
                  IC1 = ICC1
                  ICMAX = MAX(ICMAX,IC1)
                  NN = NN+1
                  IF(LEW.EQ.1) THEN
                     IF(NN.GT.LEA) THEN
                        PRINT*,'!!! MEMORY EXHAUSTION, INCREASE LEA'
                        CALL FOXSTL
                     ELSEIF(ICMAX.GT.LIA) THEN
                        PRINT*,'!!! MEMORY EXHAUSTION, INCREASE LIA'
                        CALL FOXSTL
                     ENDIF
                  ENDIF
*
                  IE1(NN) = IC1
                  IE2(NN) = IC2
                  IEO(NN) = IO1+IO2
*
                  IF(IC2.EQ.0) IA1(IC1) = NN
               ENDIF
*
               GOTO 70
            ENDIF
*
  90        CONTINUE
*
         ENDIF
         JD = -1
         GOTO 50
      ENDIF
*
 100  CONTINUE
*
      IF(NN.NE.NMMAX) THEN
         IF(LEW.EQ.0) THEN
            PRINT*,'@@@ CATASTROPHIC ERROR IN DAINI, NN # NMMAX '
            IU = -1
         ELSE
            NMMAX = NN
         ENDIF
      ENDIF
*
      CC(NBEG(INM)) = NMMAX
      NTYP(INM) = NRE
*
      ICN = 0
      IED1 = 1
      DO 110 I=1,IESP
      IF(ICN.GT.0) IED1 = IED1*IBASE
      IED(I) = IED1
      ICN = ICN+1
 110  CONTINUE
*
      ICN = 0
      IED1 = 1
      DO 120 I=IESP+1,NV
      IF(ICN.GT.0) IED1 = IED1*IBASE
      IED(I) = IED1
      ICN = ICN+1
 120  CONTINUE
*
      DO 140 I=1,NMMAX
      JJJ = IA1(IE1(I))+IA2(IE2(I))
      IF(JJJ.NE.I) THEN
         PRINT*,'@@@ CATASTROPHIC ERROR IN DAINI, '//
     *          'ARRAYS IE1,IE2,IA1,IA2 '
         PRINT*,'    INCOMPATIBLE (CF DAINI.DAT) AT I = ',I
         IU = -1
      ENDIF
*
      CALL DAENC(IE1(I),IE2(I),JJ)
*
      IJB(I) = 0
      DO 130 L=1,NV
      IF(JJ(L).NE.IJJ(L,I)) THEN
         PRINT*,'@@@ CATASTROPHIC ERROR IN DAINI, '//
     *          'FUNCTION IJJ(J,I) INCOMPATIBLE AT J = ',L,', I = ',I
         IU = -1
      ENDIF
      IF(MOD(JJ(L),2).EQ.1) GOTO 140
      IJB(I) = IJB(I)+1
 130  CONTINUE
 140  CONTINUE
*
      DO 150 I=1,LNV
 150  JJ(I) = 0
      DO 160 I=1,NV
      JJ(I) = IEW(I)
      CALL DADEC(JJ,ID1(I),ID2(I))
      JJ(I) = 0
 160  CONTINUE
*
      NFATH(1) = 0
*     PRINT*,'I, NFATH',1,NFATH(0)
      DO 180 I=2,NMMAX
      IC1 = IE1(I)
      IC2 = IE2(I)
      CALL DAENC(IC1,IC2,JJ)
      DO 170 L=NV,1,-1
      IF(JJ(L).EQ.0) GOTO 170
      IT1 = IC1-ID1(L)
      IT2 = IC2-ID2(L)
      NFATH(I) = IA1(IT1)+IA2(IT2)
*     PRINT*,'I, NFATH',I,NFATH(I)
      IF(NFATH(I).EQ.0) PRINT*,'$$$ WARNING IN DAINI, NFATH = 0'
      GOTO 180
 170  CONTINUE
 180  CONTINUE
*
      LENTM = NMMAX+2*(NVMAX+NOMAX+2)+1
*
      LFLT = 0
*
      IF(IU.EQ.0) RETURN
*
      PRINT*,'   ARRAY SETUP DONE, BEGIN PRINTING'
*
      IOUT = 32
      OPEN(IOUT,FILE='DAINI.DAT',STATUS='UNKNOWN')
*
      WRITE(IOUT,'(/A/A/)') ' ARRAYS IE1, IE2, IEO, I(1) THROUGH I(LNV)'
     *                      ,' ****************************************'
      DO 190 I=1,NMMAX
      CALL DAENC(IE1(I),IE2(I),JJ)
      WRITE(IOUT,'(1X,I5,2X,3I6,2X,100(5I2,1X))')
     *                        I,IE1(I),IE2(I),IEO(I),(JJ(L),L=1,NV)
 190  CONTINUE
*
      WRITE(IOUT,'(/A/A/)') ' ARRAYS IA1,IA2',' **************'
      DO 200 I=0,ICMAX
      WRITE(IOUT,'(3I10)') I,IA1(I),IA2(I)
 200  CONTINUE
*
      IF(IU.EQ.-1) CALL FOXSTP(1)
*
      RETURN
      END
*
      SUBROUTINE INRNDO
*     *****************
*
*     THIS SUBROUTINE COMPUTES ROUNDING PARAMETERS
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----ROUNDING --------------------------------------------------------------
      INTEGER INTINI
      DOUBLE PRECISION RINUP(2),RINDN(2),FINSR
      COMMON /INCOM/ RINUP,RINDN,FINSR,INTINI
*----------------------------------------------------------------------------
*
      SAVE RINMIN,MANTIS,MINEXP
*
      DOUBLE PRECISION EPSI(16500)
*     REPLACE 16500 BY 1100 IF QUADRUPLE PRECISION IS NOT NEEDED
*
      IF(INTINI.EQ.0) RETURN
*
*     INITIALIZING THE VARIABLES FOR ROUNDING
*
      IF(INTINI.EQ.-1) THEN
         INTINI = 1
         MANTIS = 1
         EPSI(MANTIS) = 1.D0
 40      IF(1.D0.LT.1.D0+EPSI(MANTIS)) THEN
            EPSI(MANTIS+1) = EPSI(MANTIS)/2.D0
            MANTIS = MANTIS+1
            GOTO 40
         ENDIF
         MANTIS = MANTIS-1
*
*        TEST NUMBERS TO CHECK THE ROUNDING
*
         A = SQRT(2.D0)
         B = 1/EXP(1.D0)
         C = EXP(20.D0)
*
         RINUP(1) = 1.D0+EPSI(MANTIS)
         RINDN(1) = 1.D0-EPSI(MANTIS)
*
 50      IF((A*RINUP(1).LE.A).OR.(B*RINUP(1).LE.B).OR.
     *      (C*RINUP(1).LE.C).OR.(C*RINDN(1).GE.C).OR.
     *      (B*RINDN(1).GE.B).OR.(A*RINDN(1).GE.A)) THEN
            MANTIS = MANTIS-1
            RINUP(1) = 1.D0+EPSI(MANTIS)
            RINDN(1) = 1.D0-EPSI(MANTIS)
            GOTO 50
         ENDIF
*
         MINEXP = 1
         EPSI(MINEXP) = 1.D0
 60      IF(EPSI(MINEXP)/2.D0.GT.0.D0) THEN
            EPSI(MINEXP+1) = EPSI(MINEXP)/2.D0
            MINEXP = MINEXP+1
            GOTO 60
         ENDIF
*
         MINEXP = MINEXP-MANTIS
         RINMIN = EPSI(MINEXP)/RINDN(1)
         EPSMAC = RINMIN
*
*        OPEN(17,FILE='INRNDO.DAT',STATUS='UNKNOWN')
*        WRITE(17,*) '*** MANTISSA,EXPMIN,MINFLOAT,ROUNDUP,ROUNDDOWN'
*        WRITE(17,*) MANTIS
*        WRITE(17,*) MINEXP
*        WRITE(17,'(G27.20E3)') RINMIN
*        WRITE(17,'(G27.20E3)') RINUP(1)
*        WRITE(17,'(G27.20E3)') RINDN(1)
*        CLOSE(17)
*
         EPSM = EPSI(MANTIS)
         EPS  = EPSM*1.D-4
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DANOTW(IMJ,INJ)
*     **************************
*
*     THIS SUBROUTINE SETS THE ORDER WEIGHTING FACTORS, WHICH ARE STORED
*     IN AN ARRAY IMJ WITH THE SIZE INJ
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INJ).NE.NRE) CALL FOXNTY(INJ)
      NJ = NINT(CC(NBEG(INJ)))
*
      IF(NJ.LE.0) THEN
         PRINT*,'$$$ ERROR IN DANOTW, 2ND ARGUMENT IS NOT POSITIVE'
         CALL FOXDEB
      ENDIF
*
      DO 10 J=1,LNV
      IEW(J) = 1
  10  CONTINUE
*
      DO 20 J=1,MIN(NJ,LNV)
      IF(NTYP(IMJ+J).NE.NRE) CALL FOXNTY(IMJ+J)
      IEW(J) = NINT(CC(NBEG(IMJ+J)))
      IF(IEW(J).LE.0) THEN
         PRINT*,
     *   '$$$ ERROR IN DANOTW, ARRAY ELEMENT (',J,') IS NOT POSITIVE'
         CALL FOXDEB
      ENDIF
  20  CONTINUE
*
      LEWI = 1
*
      RETURN
      END
*
      SUBROUTINE DANOT(INON)
*     **********************
*
*     THIS SUBROUTINE SETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INON).NE.NRE) CALL FOXNTY(INON)
      NON = NINT(CC(NBEG(INON)))
      NOCUT = MIN(NON,NOMAX)
*
      RETURN
      END
*
      SUBROUTINE DAEPS(IEPS)
*     **********************
*
*     THIS SUBROUTINE SETS THE DA ZERO TOLERANCE EPSILON TO A NEW VALUE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(IEPS).NE.NRE) CALL FOXNTY(IEPS)
      EPS = CC(NBEG(IEPS))
*
      RETURN
      END
*
      SUBROUTINE DAEPSM(INC)
*     **********************
*
*     THIS SUBROUTINE REPORTS THE CURRENT DA ZERO TOLERANCE EPSILON IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = EPS
*
      RETURN
      END
*
      SUBROUTINE EPSMIN(INC)
*     **********************
*
*     THIS SUBROUTINE REPORTS THE SMALLEST NUMBER REPRESENTABLE ON MACINE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = EPSMAC
*
      RETURN
      END
*
      SUBROUTINE TMTOL(ITOL)
*     **********************
*
*     THIS SUBROUTINE SETS THE TM REMAINDER MAX TOLERANCE TO A NEW VALUE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(ITOL).NE.NRE) CALL FOXNTY(ITOL)
      TOLTMR = CC(NBEG(ITOL))
*
      RETURN
      END
*
      SUBROUTINE TMTOLR(INC)
*     **********************
*
*     THIS SUBROUTINE REPORTS THE CURRENT TM REMAINDER MAX TOLERANCE IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = TOLTMR
*
      RETURN
      END
*%%
      SUBROUTINE DA(IIV,INC)
*     **********************
*
*     THIS SUBROUTINE CREATES AN IDENTITY DA VARIABLE INC OF IV-TH
*     INDEPENDENT VARIABLE.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(IIV).NE.NRE) CALL FOXNTY(IIV)
*
      IV = NINT(CC(NBEG(IIV)))
      IF(IV.LT.1.OR.IV.GT.NVMAX) THEN
         PRINT*,'$$$ ERROR IN DA, VARIABLE IS OUT OF RANGE'
         PRINT*,'    VARIABLE INDEX, BOUND = ',IV,',',NVMAX
         CALL FOXDEB
      ENDIF
*
      CALL DAVAR(INC,0.D0,IV)
*
      RETURN
      END
*
      SUBROUTINE DACNST(INA,INC)
*     **************************
*
*     THIS SUBROUTINE STORES THE CONSTANT PART OF DA VARIABLE INA IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
*
      NTYP(INC) = NRE
      IC = NBEG(INC)
      IPOA = NBEG(INA)
*
      IF(NEND(INA).EQ.IPOA-1) THEN
         CC(IC) = 0.D0
         RETURN
      ELSE
         IF(NC(IPOA).EQ.1) THEN
            CC(IC) = CC(IPOA)
            RETURN
         ELSE
            CC(IC) = 0.D0
            RETURN
         ENDIF
      ENDIF
      END
*
      DOUBLE PRECISION FUNCTION RE(INA)
*     *********************************
*
*     THIS FUNCTION RETURNS THE ZEROTH ORDER COMPONENT OF THE DA VARIABLE
*     INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
*
      IPOA = NBEG(INA)
      IF(NEND(INA).EQ.IPOA-1) THEN
         RE = 0.D0
         RETURN
      ELSE
         IF(NC(IPOA).EQ.1) THEN
            RE = CC(IPOA)
            RETURN
         ELSE
            RE = 0.D0
            RETURN
         ENDIF
      ENDIF
      END
*
      SUBROUTINE DANFI(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PRINT*,'$$$ ERROR, OPERATION NOT POSSIBLE: ',INA,INB,INC
      CALL FOXDEB
      RETURN
      END
*
      SUBROUTINE DACOP(INA,INB)
*     *************************
*
*     THIS SUBROUTINE COPIES THE DA VECTOR A TO THE DA VECTOR B
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IB = NBEG(INB)-1
*
      TMS = 0.D0
*
      DO 10 IA = NBEG(INA),NEND(INA)
      IF(IEO(NC(IA)).GT.NOCUT) GOTO 10
      CCC = CC(IA)
      IF(ABS(CCC).LT.EPS) THEN
         TMS = TMS+ABS(CCC)
         GOTO 10
      ENDIF
      IB = IB+1
      CC(IB) = CCC
      NC(IB) = NC(IA)
 10   CONTINUE
*
      NTYP(INB) = NDA
      NEND(INB) = IB
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('DACOP')
*
      RETURN
      END
*
      SUBROUTINE LDA(INA,INC)
*     ***********************
*
*     THIS SUBROUTINE RETURNS THE ALLOCATION LENGTH OF DA VECTOR IN INC.
*     THE VECTOR INA SUPPLIES THE ORDER AND THE DIMENSION OF DA.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      CNO = CC(IA)
      CNV = CC(IA+1)
      CMIN = MIN(CNO,CNV)
      CMAX = MAX(CNO,CNV)
*
      A0 = CC(NBEG(INA))
      IF(CMIN.LE.0.9D0) THEN
         PRINT*,'$$$ WARNING IN LDA, '//
     *          'TWO COMPONENTS OF ARGUMENT MUST BE POSITIVE INTEGERS.'
         CC(NBEG(INC)) = 1.D0
      ELSEIF(CMAX.LE.14.D0) THEN
         NO = NINT(CNO)
         NV = NINT(CNV)
         CALL DANUM(NO,NV,NUMDA)
         CALL DANUMD(NO,NV,CNUMDA)
         IF(ABS(NUMDA-CNUMDA).GT.1.D-12) THEN
            CC(NBEG(INC)) = DBLE(NUMDA)
         ELSE
            CC(NBEG(INC)) = CNUMDA
         ENDIF
      ELSEIF(CMIN.LE.500.D0) THEN
C        CNO=CNV=514 RESULTS IN 0.7156051054877889+308
C        CNO=CNV=515 RESULTS IN INFINITY IN DOUBLE PRECISION
C        THERE MAY BE CASES TO RESULT IN INFINITY
         CNUMDA = 1.D0
         DO 10 I=1,NINT(CMIN)
         XI = I
 10      CNUMDA = (CNUMDA*(CMAX+XI))/XI
         CC(NBEG(INC)) = CNUMDA
      ELSE
         PRINT*,'$$$ WARNING IN LDA, '//
     *          'GIVEN NUMBERS ARE TOO LARGE TO COMPUTE. RETURNING 1.'
         CC(NBEG(INC)) = 1.D0
      ENDIF
*
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE DAZERO(INA,INB)
*     **************************
*
*     THIS SUBROUTINE SETS B TO DA ZERO
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      NEND(INB) = NBEG(INB)-1
      NTYP(INB) = NDA
*
      RETURN
      END
*
      SUBROUTINE DAADA(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS TWO DA VECTORS INA AND INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NDA) CALL FOXNTY(INB)
      IF((INA.EQ.INC).OR.(INB.EQ.INC)) CALL DANFI(INA,INB,INC)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)-1
      IAMAX = NEND(INA)
      IBMAX = NEND(INB)
      NA = NC(IA)
      NB = NC(IB)
*
      TMT = 0.D0
      TMS = 0.D0
*
      IF(IA.GT.IAMAX) THEN
         DO 17 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = CC(IS)
  17     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAADA')
         RETURN
      ELSEIF(IB.GT.IBMAX) THEN
         DO 18 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  18     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAADA')
         RETURN
      ENDIF
*
      IF(NA-NB) 30,20,40
*
*     ADDING TWO TERMS
*     ****************
*
  20  CONTINUE
      CCC = CC(IA)+CC(IB)
      TMT = TMT+ABS(CC(IA))+ABS(CC(IB))
      IF(ABS(CCC).LT.EPS) THEN
         TMS = TMS+ABS(CCC)
         GOTO 25
      ENDIF
      IC = IC+1
      CC(IC) = CCC
      NC(IC) = NA
  25  CONTINUE
      IA = IA+1
      IB = IB+1
      IF(IA.GT.IAMAX) THEN
         DO 27 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = CC(IS)
  27     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAADA')
         RETURN
      ELSEIF(IB.GT.IBMAX) THEN
         DO 28 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  28     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAADA')
         RETURN
      ENDIF
      NA = NC(IA)
      NB = NC(IB)
      IF(NA-NB) 30,20,40
*
*     STORING TERM A
*     **************
*
  30  CONTINUE
      IC = IC+1
      CC(IC) = CC(IA)
      NC(IC) = NA
      IA = IA+1
      IF(IA.GT.IAMAX) THEN
         DO 35 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = CC(IS)
  35     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAADA')
         RETURN
      ENDIF
      NA = NC(IA)
      IF(NA-NB) 30,20,40
*
*     STORING TERM B
*     **************
*
  40  CONTINUE
      IC = IC+1
      CC(IC) = CC(IB)
      NC(IC) = NB
      IB = IB+1
      IF(IB.GT.IBMAX) THEN
         DO 45 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  45     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAADA')
         RETURN
      ENDIF
      NB = NC(IB)
      IF(NA-NB) 30,20,40
*
      END
*%%
      SUBROUTINE DASDA(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE SUBTRACTS TWO DA VECTORS INA AND INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NDA) CALL FOXNTY(INB)
      IF((INA.EQ.INC).OR.(INB.EQ.INC)) CALL DANFI(INA,INB,INC)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)-1
      IAMAX = NEND(INA)
      IBMAX = NEND(INB)
      NA = NC(IA)
      NB = NC(IB)
*
      TMT = 0.D0
      TMS = 0.D0
*
      IF(IA.GT.IAMAX) THEN
         DO 17 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = -CC(IS)
  17     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DASDA')
         RETURN
      ELSEIF(IB.GT.IBMAX) THEN
         DO 18 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  18     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DASDA')
         RETURN
      ENDIF
*
      IF(NA-NB) 30,20,40
*
*     SUBTRACTING TWO TERMS
*     *********************
*
  20  CONTINUE
      CCC = CC(IA)-CC(IB)
      TMT = TMT+ABS(CC(IA))+ABS(CC(IB))
      IF(ABS(CCC).LT.EPS) THEN
         TMS = TMS+ABS(CCC)
         GOTO 25
      ENDIF
      IC = IC+1
      CC(IC) = CCC
      NC(IC) = NA
  25  CONTINUE
      IA = IA+1
      IB = IB+1
      IF(IA.GT.IAMAX) THEN
         DO 27 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = -CC(IS)
  27     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DASDA')
         RETURN
      ELSEIF(IB.GT.IBMAX) THEN
         DO 28 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  28     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DASDA')
         RETURN
      ENDIF
      NA = NC(IA)
      NB = NC(IB)
      IF(NA-NB) 30,20,40
*
*     STORING TERM A
*     **************
*
  30  CONTINUE
      IC = IC+1
      CC(IC) = CC(IA)
      NC(IC) = NA
      IA = IA+1
      IF(IA.GT.IAMAX) THEN
         DO 35 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = -CC(IS)
  35     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DASDA')
         RETURN
      ENDIF
      NA = NC(IA)
      IF(NA-NB) 30,20,40
*
*     STORING TERM B
*     **************
*
  40  CONTINUE
      IC = IC+1
      CC(IC) = -CC(IB)
      NC(IC) = NB
      IB = IB+1
      IF(IB.GT.IBMAX) THEN
         DO 45 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  45     NC(IC) = NC(IS)
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DASDA')
         RETURN
      ENDIF
      NB = NC(IB)
      IF(NA-NB) 30,20,40
*
      END
*%%
      SUBROUTINE DALIN(INA,AFAC,INB,BFAC,INC)
*     ***************************************
*
*     THIS SUBROUTINE COMPUTES THE LINEAR COMBINATION
*     C = AFAC*INA + BFAC*INB OF THE DA VECTORS INA AND INB.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION AFAC,BFAC,CCC
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NDA) CALL FOXNTY(INB)
      IF((INA.EQ.INC).OR.(INB.EQ.INC)) CALL DANFI(INA,INB,INC)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)-1
      IAMAX = NEND(INA)
      IBMAX = NEND(INB)
      NA = NC(IA)
      NB = NC(IB)
*
      IF(IA.GT.IAMAX) THEN
         ISMIN = IB
         ISMAX = IBMAX
         COPF  = BFAC
         GOTO 50
      ELSEIF(IB.GT.IBMAX) THEN
         ISMIN = IA
         ISMAX = IAMAX
         COPF  = AFAC
         GOTO 50
      ENDIF
*
      IF(NA-NB) 30,20,40
*
*     ADDING TWO TERMS
*     ****************
*
  20  CONTINUE
      CCC = CC(IA)*AFAC+CC(IB)*BFAC
      IF(ABS(CCC).LT.EPS) GOTO 25
      IC = IC+1
      CC(IC) = CCC
      NC(IC) = NA
  25  CONTINUE
      IA = IA+1
      IB = IB+1
      IF(IA.GT.IAMAX) THEN
         ISMIN = IB
         ISMAX = IBMAX
         COPF  = BFAC
         GOTO 50
      ELSEIF(IB.GT.IBMAX) THEN
         ISMIN = IA
         ISMAX = IAMAX
         COPF  = AFAC
         GOTO 50
      ENDIF
      NA = NC(IA)
      NB = NC(IB)
      IF(NA-NB) 30,20,40
*
*     STORING TERM A
*     **************
*
  30  CONTINUE
      CCC = CC(IA)*AFAC
      IF(ABS(CCC).LT.EPS) GOTO 35
      IC = IC+1
      CC(IC) = CCC
      NC(IC) = NA
  35  CONTINUE
      IA = IA+1
      IF(IA.GT.IAMAX) THEN
         ISMIN = IB
         ISMAX = IBMAX
         COPF  = BFAC
         GOTO 50
      ENDIF
      NA = NC(IA)
      IF(NA-NB) 30,20,40
*
*     STORING TERM B
*     **************
*
  40  CONTINUE
      CCC = CC(IB)*BFAC
      IF(ABS(CCC).LT.EPS) GOTO 45
      IC = IC+1
      CC(IC) = CCC
      NC(IC) = NB
  45  CONTINUE
      IB = IB+1
      IF(IB.GT.IBMAX) THEN
         ISMIN = IA
         ISMAX = IAMAX
         COPF  = AFAC
         GOTO 50
      ENDIF
      NB = NC(IB)
      IF(NA-NB) 30,20,40
*
*     COPYING THE REST
*     ****************
*
  50  CONTINUE
      DO 60 IS=ISMIN,ISMAX
      CCC = CC(IS)*COPF
      IF(ABS(CCC).LT.EPS) GOTO 60
      IC = IC+1
      CC(IC) = CCC
      NC(IC) = NC(IS)
  60  CONTINUE
*
      NTYP(INC) = NDA
      NEND(INC) = IC
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DALIN')
*
      RETURN
      END
*%%
      SUBROUTINE DAMDA(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION NBNO(0:LNO),NENO(0:LNO)
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NDA) CALL FOXNTY(INB)
*
      CALL DACLR
*
*     PRESORTING B BY ORDER
*     *********************
*
      JMEM = IMEM
      DO 20 I=0,NOMAX
      NBNO(I) = JMEM+1
      NENO(I) = JMEM
      JMEM = JMEM+NMMAX
  20  CONTINUE
*
      CALL MEMCHK(JMEM)
*
      DO 40 IB=NBEG(INB),NEND(INB)
*
      NCIB = NC(IB)
      NOIB = IEO(NCIB)
      IPOS = NENO(NOIB)+1
      NENO(NOIB) = IPOS
*
      CC(IPOS) = CC(IB)
      NC(IPOS) = NC(IB)
*
  40  CONTINUE
*
      TMT = 0.D0
      TMS = 0.D0
*
*     PERFORMING ACTUAL MULTIPLICATION
*     ********************************
*
      IF(ITM.NE.1) THEN
*
         DO 100 IA=NBEG(INA),NEND(INA)
*
*        CHANGING ORDER OF NESTED LOOPS DOES NOT IMPROVE DA COEFFICIENT ACCURACY
*        SINCE ALL CONTRIBUTIONS TO TARGET COEFFICIENTS ARE SIMILAR SIZE
*
         NCIA = NC(IA)
         I1IA = IE1(NCIA)
         I2IA = IE2(NCIA)
         CCIA = CC(IA)
*
         DO 100 NOIB = 0,NOCUT-IEO(NCIA)
         DO 100 IB = NBNO(NOIB),NENO(NOIB)
*
         NCIB = NC(IB)
         ICC = IA2(I2IA+IE2(NCIB))+IA1(I1IA+IE1(NCIB))
         CDA(ICC) = CDA(ICC)+CCIA*CC(IB)
*
 100     CONTINUE
*
      ELSEIF(ITM.EQ.1) THEN
*
         DO 200 IA=NBEG(INA),NEND(INA)
*
*        CHANGING ORDER OF NESTED LOOPS DOES NOT IMPROVE DA COEFFICIENT ACCURACY
*        SINCE ALL CONTRIBUTIONS TO TARGET COEFFICIENTS ARE SIMILAR SIZE
*
         NCIA = NC(IA)
         I1IA = IE1(NCIA)
         I2IA = IE2(NCIA)
         CCIA = CC(IA)
*
         DO 200 NOIB = 0,NOCUT-IEO(NCIA)
         DO 200 IB = NBNO(NOIB),NENO(NOIB)
*
         NCIB = NC(IB)
         ICC = IA2(I2IA+IE2(NCIB))+IA1(I1IA+IE1(NCIB))
         CK = CCIA*CC(IB)
         CDAA = CDA(ICC)
         TMT = TMT+ABS(CK)+MAX(ABS(CDAA),ABS(CK))
*
*        CHANGING ORDER OF NESTED LOOP DOES NOT IMPROVE TM REMAINDER ACCURACY
*        IT ONLY IMRPOVES ACCURACY OF TMT, BUT NOT ITS SIZE
*
C        THE FOLLOWING GIVES SMALLER TMT, BUT TAKES MORE CPU TIME.
C
C        IF(CDAA.EQ.0.D0) THEN
C           TMT = TMT+ABS(CK)
C        ELSE
C           TMT = TMT+ABS(CK)+MAX(ABS(CDAA),ABS(CK))
C        ENDIF
         CDA(ICC) = CDAA+CK
*
 200     CONTINUE
*
      ENDIF
*
      CALL DAPAC(INC)
*
      RETURN
      END
*
      SUBROUTINE DADDA(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE DIVIDES THE DA VECTORS INA AND INB.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(IDIV,1,NMMAX)
      CALL DAMUI(INB,IDIV)
      CALL DAMDA(INA,IDIV,INC)
      CALL FOXDAL(IDIV,1)
*
      RETURN
      END
*
      SUBROUTINE READA(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS THE REAL INA AND THE DA VECTOR INB.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL DAARE(INB,INA,INC)
      RETURN
      END
*
      SUBROUTINE DAARE(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE ADDS THE DA VECTOR INA AND THE REAL INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(INA.EQ.INC)       CALL DANFI(INA,INB,INC)
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IBEA = NBEG(INA)
      IENA = NEND(INA)
      IBEC = NBEG(INC)
*
      TMT = 0.D0
      TMS = 0.D0
*
      CKON = CC(NBEG(INB))
*
      IF(IENA.LT.IBEA) THEN
         IF(ABS(CKON).GT.EPS) THEN
            CC(IBEC) = CKON
            TMT = ABS(CKON)
            NC(IBEC) = 1
            NEND(INC) = IBEC
            NTYP(INC) = NDA
            RETURN
         ELSE
            TMS = ABS(CKON)
            NEND(INC) = IBEC-1
            NTYP(INC) = NDA
            RETURN
         ENDIF
      ELSEIF(NC(IBEA).EQ.1) THEN
         CCC = CC(IBEA)+CKON
         TMT = ABS(CC(IBEA))+ABS(CKON)
         IBBA = IBEA+1
      ELSE
         CCC = CKON
         IBBA = IBEA
      ENDIF
*
      IF(ABS(CCC).GT.EPS) THEN
         TMT = TMT+ABS(CCC)
         CC(IBEC) = CCC
         NC(IBEC) = 1
         IC = IBEC
      ELSE
         TMS = ABS(CCC)
         IC = IBEC-1
      ENDIF
*
      DO 10 I=IBBA,IENA
      IC = IC+1
      CC(IC) = CC(I)
      NC(IC) = NC(I)
  10  CONTINUE
*
      IF(IC.GT.NMAX(INC)) CALL FOXERV('DAARE')
*
      NEND(INC) = IC
      NTYP(INC) = NDA
*
      RETURN
      END
*
      SUBROUTINE DASRE(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE SUBTRACTS THE REAL INB FROM THE DA VECTOR INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      ICC = NBEG(INB)
      CC(ICC) = -CC(ICC)
      CALL DAARE(INA,INB,INC)
      CC(ICC) = -CC(ICC)
*
      RETURN
      END
*
      SUBROUTINE RESDA(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE SUBTRACTS THE VECTOR INB FROM THE CONSTANT INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      ICC = NBEG(INA)
      CC(ICC) = -CC(ICC)
      CALL DAARE(INB,INA,INC)
      CALL DAADI(INC,INC)
      CC(ICC) = -CC(ICC)
*
      RETURN
      END
*%%
      SUBROUTINE DAADI(INA,INB)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE ADDITIVE INVERSE OF INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
*
      IB = NBEG(INB)-1
*
      DO 100 IA = NBEG(INA),NEND(INA)
      IB = IB+1
      CC(IB) = -CC(IA)
      NC(IB) = NC(IA)
 100  CONTINUE
*
      NTYP(INB) = NDA
      NEND(INB) = IB
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('DAADI')
*
      RETURN
      END
*
      SUBROUTINE REMDA(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL DAMRE(INB,INA,INC)
*
      RETURN
      END
*
      SUBROUTINE DAMRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MULTIPLIES DA VECTOR INA WITH REAL INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      CKON = CC(NBEG(INB))
*
      IC = NBEG(INC)-1
*
      TMT = 0.D0
      TMS = 0.D0
*
      IF(ITM.NE.1) THEN
*
         DO 100 IA = NBEG(INA),NEND(INA)
         CCC = CC(IA)*CKON
         IF(ABS(CCC).LT.EPS) GOTO 100
         IC = IC+1
         CC(IC) = CCC
         NC(IC) = NC(IA)
 100     CONTINUE
*
      ELSEIF(ITM.EQ.1) THEN
*
         DO 200 IA = NBEG(INA),NEND(INA)
         CCC = CC(IA)*CKON
         IF(ABS(CCC).LT.EPS) THEN
            TMS = TMS+ABS(CCC)
            GOTO 200
         ELSE
            TMT = TMT+ABS(CCC)
         ENDIF
         IC = IC+1
         CC(IC) = CCC
         NC(IC) = NC(IA)
 200     CONTINUE
*
      ENDIF
*
      NTYP(INC) = NDA
      NEND(INC) = IC
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAMRE')
*
      RETURN
      END
*%%
      SUBROUTINE DADRE(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE DIVIDES THE DA VECTOR INA BY THE REAL INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      CKON = CC(NBEG(INB))
      IF(CKON.EQ.0.D0) THEN
         PRINT*, '$$$ ERROR IN DADRE, DIVISOR IS ZERO'
         CALL FOXDEB
      ENDIF
*
      IC = NBEG(INC)-1
*
      TMT = 0.D0
      TMS = 0.D0
*
      IF(ITM.NE.1) THEN
*
         DO 100 IA = NBEG(INA),NEND(INA)
         CCC = CC(IA)/CKON
         IF(ABS(CCC).LT.EPS) GOTO 100
         IC = IC+1
         CC(IC) = CCC
         NC(IC) = NC(IA)
 100     CONTINUE
*
      ELSEIF(ITM.EQ.1) THEN
*
         DO 200 IA = NBEG(INA),NEND(INA)
         CCC = CC(IA)/CKON
         IF(ABS(CCC).LT.EPS) THEN
            TMS = TMS+ABS(CCC)
            GOTO 200
         ELSE
            TMT = TMT+ABS(CCC)
         ENDIF
         IC = IC+1
         CC(IC) = CCC
         NC(IC) = NC(IA)
 200     CONTINUE
*
      ENDIF
*
      NTYP(INC) = NDA
      NEND(INC) = IC
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DADRE')
*
      RETURN
      END
*
      SUBROUTINE REDDA(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE DIVIDES THE REAL INA BY THE DA VECTOR INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL DAMUI(INB,INC)
      CALL DAMRE(INC,INA,INC)
*
      RETURN
      END
*%%
      SUBROUTINE DAMUI(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE MULTIPLICATIVE INVERSE TO INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
      IF(A0.EQ.0.D0) THEN
         PRINT*,'$$$ ERROR IN DAMUI, '//
     *          'MULTIPLICATIVE INVERSE DOES NOT EXIST FOR ',INA
         CALL FOXDEB
      ENDIF
*
      XF(0) = 1.D0/A0
      DO 10 I=1,NOCUT
  10  XF(I) = -XF(I-1)
*
      CALL DAFUN(INA,A0,XF,NOCUT,INC)
*
      RETURN
      END
*
      SUBROUTINE DASQR(INA,INC)
*     *************************
*
*     THIS SUBROUTINE SQUARES INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL DAMDA(INA,INA,INC)
*
      RETURN
      END
*
      SUBROUTINE DASQRT(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE SQUARE ROOT OF A DA QUANTITY
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
*
      IF(A0.LE.0.D0) THEN
         PRINT*,'$$$ ERROR IN DASQRT, SQUARE ROOT DOES NOT EXIST FOR '
     *         ,INA
         CALL FOXDEB
      ENDIF
*
      XF(0) = SQRT(A0)
      DO 10 I=1,NOCUT
  10  XF(I) = -XF(I-1)/DBFLOAT(2*I)*DBFLOAT(2*I-3)
*
      CALL DAFUN(INA,A0,XF,NOCUT,INC)
      RETURN
      END
*
      SUBROUTINE DAISRT(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE INVERSE SQUARE ROOT OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
*
      IF(A0.LE.0.D0) THEN
         PRINT*,'$$$ ERROR IN DAISRT, '//
     *          'INVERSE OF SQUARE ROOT DOES NOT EXIST FOR ',INA
         CALL FOXDEB
      ENDIF
*
      XF(0) = 1.D0/SQRT(A0)
      DO 10 I=1,NOCUT
  10  XF(I) = -XF(I-1)/DBFLOAT(2*I)*DBFLOAT(2*I-1)
*
      CALL DAFUN(INA,A0,XF,NOCUT,INC)
      RETURN
      END
*
      SUBROUTINE DAISR3(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE INVERSE TO THE POWER 3/2 OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
*
      IF(A0.LE.0.D0) THEN
         PRINT*,'$$$ ERROR IN DAISR3, '//
     *          'INVERSE TO THE POWER 3/2 DOES NOT EXIST FOR ',INA
         CALL FOXDEB
      ENDIF
*
      XF(0) = 1.D0/(A0*SQRT(A0))
      DO 10 I=1,NOCUT
  10  XF(I) = -XF(I-1)/DBFLOAT(2*I)*DBFLOAT(2*I+1)
*
      CALL DAFUN(INA,A0,XF,NOCUT,INC)
      RETURN
      END
*
      SUBROUTINE DAEXP(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE EXPONENTIAL OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
*
      XF(0) = DEXP(A0)
      DO 10 I=1,NOCUT
  10  XF(I) = XF(I-1)/DBFLOAT(I)
*
      CALL DAFUN(INA,1.D0,XF,NOCUT,INC)
      RETURN
      END
*
      SUBROUTINE DALOG(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE LOG OF INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
*
      IF(A0.LE.0.D0) THEN
         PRINT*,'$$$ ERROR IN DALOG, LOGARITHM DOES NOT EXIST FOR ',INA
         CALL FOXDEB
      ENDIF
*
      XF(0) = DLOG(A0)
      XF(1) = 1.D0
      DO 10 I=2,NOCUT
  10  XF(I) = -XF(I-1)/DBFLOAT(I)*DBFLOAT(I-1)
*
      CALL DAFUN(INA,A0,XF,NOCUT,INC)
*
      RETURN
      END
*
      SUBROUTINE DASINE(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE SINE OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
*
      XF(0) = DSIN(A0)
      XF(1) = DCOS(A0)
      DO 10 I=2,NOCUT
  10  XF(I) = -XF(I-2)/DBFLOAT(I)/DBFLOAT(I-1)
*
      CALL DAFUN(INA,1.D0,XF,NOCUT,INC)
*
      RETURN
      END
*
      SUBROUTINE DACOSE(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE COSINE OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
*
      XF(0) = DCOS(A0)
      XF(1) = -DSIN(A0)
      DO 10 I=2,NOCUT
  10  XF(I) = -XF(I-2)/DBFLOAT(I)/DBFLOAT(I-1)
*
      CALL DAFUN(INA,1.D0,XF,NOCUT,INC)
*
      RETURN
      END
*
      SUBROUTINE DADTAN(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE TAN OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ISIN,1,NMMAX)
      CALL FOXALL(ICOS,1,NMMAX)
      CALL DASINE(INA,ISIN)
      CALL DACOSE(INA,ICOS)
      CALL DADDA(ISIN,ICOS,INC)
      CALL FOXDAL(ICOS,1)
      CALL FOXDAL(ISIN,1)
*
      RETURN
      END
*
      SUBROUTINE DADCOT(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE COTAN OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ISIN,1,NMMAX)
      CALL FOXALL(ICOS,1,NMMAX)
      CALL DASINE(INA,ISIN)
      CALL DACOSE(INA,ICOS)
      CALL DADDA(ICOS,ISIN,INC)
      CALL FOXDAL(ICOS,1)
      CALL FOXDAL(ISIN,1)
*
      RETURN
      END
*
      SUBROUTINE DAASIN(INA,INC)
*     **************************
*     USE 4.4.32, AND 4.4.40 PAGE 80 HANDBOOK MATH. FUNCTIONS. LET R AND N BE
*     THE REAL AND NILPOTENT PARTS OF INA. CHOOSE
*     Z1 = R, Z2 = SQRT(1-R*R)*(R+N) - R*SQRT(1-(R+N)*(R+N))
*     THEN Z2 IS NILPOTENT, INA CAN BE WRITTEN AS
*     INA = Z1*SQRT(1-Z2*Z2) + Z2*SQRT(1-Z1*Z1),
*     AND THUS INC = ASIN(INA) = ASIN(Z1) + ASIN(Z2)
*
*     ANOTHER WAY IS TO USE THE FORMULAS OF HIGHER ORDER DERIVATIVES.
*     (FORMULAS FROM P.39,41 OF MATHEMATICAL FORMULAS I, IWANAMI-ZENSHO, 1956.)
*     F = ASIN(X)
*     THE N-TH DERIVATIVES (N>=1):
*     F(N) = 1/2^(N-1)/(1-X^2)^((2*N-1)/2) * SUM_{K=0}^{N-1} (-1)^K
*            (N-1)!/K!/(N-1-K)! * (2*K-1)!! (2*N-2*K-3)!! (1-X)^K (1+X)^(N-1-K)
*     THE RECURSION FORMULA FOR N>=0:
*     F(N+2) = 1/(1-X^2) * ( (2*N+1)*X*F(N+1) + N^2*F(N) )
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      Z1 = RE(INA)
      IF(ABS(Z1).GE.1.D0) THEN
         PRINT*,'$$$ ERROR IN DAASIN, ARC SIN DOES NOT EXIST FOR ',INA
         CALL FOXDEB
      ENDIF
*
      CALL FOXALL( IZ1,   1, 1     )
      CALL FOXALL( IZ2,   1, NMMAX )
      CALL FOXALL( IASZ1, 1, 1     )
      CALL FOXALL( IREZ12,1, 1     )
      CALL FOXALL( IS1,   1, NMMAX )
      CALL FOXALL( IS2,   1, NMMAX )
      CALL FOXALL( IONE,  1, 1     )
*
      CC(NBEG(IZ1))    = Z1
      CC(NBEG(IASZ1))  = ASIN(Z1)
      CC(NBEG(IREZ12)) = SQRT(1.D0-Z1*Z1)
      CC(NBEG(IONE))    = 1.D0
      NTYP(IZ1)        = NRE
      NTYP(IASZ1)      = NRE
      NTYP(IREZ12)     = NRE
      NTYP(IONE)       = NRE
*
      CALL DAMDA(INA,INA,IS1)
      CALL RESDA(IONE,IS1,IS2)
      CALL DASQRT(IS2,IS1)
      CALL DAMRE(IS1,IZ1,IS2)
      CALL REMDA(IREZ12,INA,IS1)
      CALL DASDA(IS1,IS2,IZ2)
      IF(ABS(RE(IZ2)).GT.1.D-14)
     *   PRINT*, ' @@@ ERROR IN DAASIN AT 2 ', RE(IZ2)
*
      XF(0) = 0.D0
      XF(1) = 1.D0
*
      DO 10 I=3,NOCUT,2
  10  XF(I) = XF(I-2)*DBFLOAT(I-2)*DBFLOAT(I-2)/DBFLOAT(I-1)/DBFLOAT(I)
      DO 20 I=2,NOCUT,2
  20  XF(I) = 0.D0
*
      CALL DAFUN(IZ2,1.D0,XF,NOCUT,IS1)
      CALL READA(IASZ1,IS1,INC)
*
      CALL FOXDAL( IONE,  1)
      CALL FOXDAL( IS2,   1)
      CALL FOXDAL( IS1,   1)
      CALL FOXDAL( IREZ12,1)
      CALL FOXDAL( IASZ1, 1)
      CALL FOXDAL( IZ2,   1)
      CALL FOXDAL( IZ1,   1)
*
      RETURN
      END
*
      SUBROUTINE DAACOS(INA,INC)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      PI2 = 1.57079632679489661923132169163975144209858469968755291048D0
*
      CALL FOXALL(ISIN,1,NMMAX)
      CALL FOXALL(IPI2,1,1    )
*
      CC(NBEG(IPI2)) = PI2
      NTYP(IPI2) = NRE
      CALL DAASIN(INA,ISIN)
      CALL RESDA(IPI2,ISIN,INC)
*
      CALL FOXDAL(IPI2,1)
      CALL FOXDAL(ISIN,1)
*
      RETURN
      END
*
      SUBROUTINE DAATAN(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE ARC TAN OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      CALL FOXALL( INUM, 1, NMMAX )
      CALL FOXALL( IDEN, 1, NMMAX )
      CALL FOXALL( IARG, 1, NMMAX )
      CALL FOXALL( IA,   1, 1 )
      CALL FOXALL( IONE, 1, 1 )
*
      NTYP(IA  ) = NRE
      NTYP(IONE) = NRE
*
      A0 = RE(INA)
      CC(NBEG(IA))   = A0
      CC(NBEG(IONE)) = 1.D0
*
      CALL DASRE(INA, IA,  INUM)
      CALL DAMRE(INA, IA,  IARG)
      CALL DAARE(IARG,IONE,IDEN)
      CALL DADDA(INUM,IDEN,IARG)
*
      XF(0) = ATAN(A0)
      DO 5  K=1,NOCUT
  5   XF(K) = 0.D0
*
      T = -1.D0
      DO 10 K=1,NOCUT,2
      T = -T
      XF(K) = T/DBFLOAT(K)
 10   CONTINUE
*
      CALL DAFUN(IARG, 1.D0, XF, NOCUT, INC )
*
      CALL FOXDAL(IONE,1)
      CALL FOXDAL(IA,  1)
      CALL FOXDAL(IARG,1)
      CALL FOXDAL(IDEN,1)
      CALL FOXDAL(INUM,1)
*
      RETURN
      END
*
      SUBROUTINE DASINH(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE HYPERBOLIC SINE OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
*
      XF(0) = DSINH(A0)
      XF(1) = DCOSH(A0)
      DO 10 I=2,NOCUT
  10  XF(I) = XF(I-2)/DBFLOAT(I)/DBFLOAT(I-1)
*
      CALL DAFUN(INA,1.D0,XF,NOCUT,INC)
*
      RETURN
      END
*
      SUBROUTINE DACOSH(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE HYPERBOLIC COSINE OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
*
      XF(0) = DCOSH(A0)
      XF(1) = DSINH(A0)
      DO 10 I=2,NOCUT
  10  XF(I) = XF(I-2)/DBFLOAT(I)/DBFLOAT(I-1)
*
      CALL DAFUN(INA,1.D0,XF,NOCUT,INC)
*
      RETURN
      END
*
      SUBROUTINE DATANH(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE HYPERBOLIC TANGENT OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ISINH,1,NMMAX)
      CALL FOXALL(ICOSH,1,NMMAX)
      CALL DASINH(INA,ISINH)
      CALL DACOSH(INA,ICOSH)
      CALL DADDA(ISINH,ICOSH,INC)
      CALL FOXDAL(ICOSH,1)
      CALL FOXDAL(ISINH,1)
*
      RETURN
      END
*
      SUBROUTINE DACOTH(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE HYPERBOLIC COTANGENT OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ISINH,1,NMMAX)
      CALL FOXALL(ICOSH,1,NMMAX)
      CALL DASINH(INA,ISINH)
      CALL DACOSH(INA,ICOSH)
      CALL DADDA(ICOSH,ISINH,INC)
      CALL FOXDAL(ICOSH,1)
      CALL FOXDAL(ISINH,1)
*
      RETURN
      END
*
      SUBROUTINE DAERF(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE ERROR FUNCTION ERF OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      A0 = RE(INA)
*
      CALL FOXALL(ISR,1,1)
      CC(NBEG(ISR)) = A0
      NTYP(ISR) = NRE
      CALL REERF(ISR,ISR)
      XF(0) = CC(NBEG(ISR))
      CALL FOXDAL(ISR,1)
*
      E1 = DEXP(-A0*A0)/SQRT(DATAN(1.D0))
      XF(1) = E1
      XF(2) = -A0*E1
      XF(3) = (-2.D0+4.D0*A0*A0)/6.D0*E1
      XF(4) = (12.D0*A0-8.D0*A0*A0*A0)/24.D0*E1
      XF(5) = (16.D0*A0*A0*A0*A0-48.D0*A0*A0+12.D0)/120.D0*E1
      XF(6) =
     *  (-32.D0*A0*A0*A0*A0*A0+160.D0*A0*A0*A0-120.D0*A0)/720.D0*E1
      IF(NOCUT.GT.6) THEN
         PRINT*,'$$$ ERROR IN DAERF, SUPPORTED ONLY UP TO NOCUT = 6'
         CALL FOXDEB
      ENDIF
*
      CALL DAFUN(INA,1.D0,XF,NOCUT,INC)
*
      PRINT*,'$$$ WARNING IN DAERF, THIS ROUTINE IS UNDER DEVELOPMENT.'
      PRINT*,'    CONTACT US AT BERZ@MSU.EDU'
*
      RETURN
      END
*
      SUBROUTINE DAFUN(INA,CNORM,XF,NP,INC)
*     *************************************
*
*     THIS SUBROUTINE EVALUATES
*     C = XF(0) + XF(1)*AA + ...+ XF(NP)*AA^NP
*     WHERE AA IS THE NON-CONSTANT PART OF (A/CNORM)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
*
      CALL FOXALL(ISC,1,NMMAX+NOMAX)
      CALL FOXALL(ICA,1,NMMAX)
      CALL FOXALL(ISR,1,1)
*
      EPSO = EPS
      EPS  = EPSMAC
*
      NTYP(ISR) = NRE
      NBSR = NBEG(ISR)
      CC(NBSR) = 1.D0/CNORM
      CALL DAMRE(INA,ISR,ICA)
*
      NBICA = NBEG(ICA)
      IF(NC(NBICA).EQ.1) NBEG(ICA) = NBICA+1
*
      NBISC = NBEG(ISC)+NP-1
      NBEG(ISC) = NBISC
*
      CC(NBSR) = XF(NP)
      CALL DAMRE(ICA,ISR,ISC)
*
      DO 10 I=NP-1,1,-1
      NBISC = NBISC-1
      NBEG(ISC) = NBISC
      CC(NBISC) = XF(I)
      NC(NBISC) = 1
      IF(I.EQ.1) EPS = EPSO
      CALL DAMDA(ICA,ISC,ISC)
  10  CONTINUE
*
      EPS = EPSO
      CC(NBSR) = XF(0)
      CALL DAARE(ISC,ISR,INC)
*
      CALL FOXDAL(ISR,1)
      CALL FOXDAL(ICA,1)
      CALL FOXDAL(ISC,1)
*
      RETURN
      END
*
      SUBROUTINE DANOR(INA,INORM)
*     ***************************
*
*     THIS SUBROUTINE COMPUTES THE NORM OF THE DA VECTOR A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      ANORM = 0.D0
      DO 100 I=NBEG(INA),NEND(INA)
      IF(IEO(NC(I)).GT.NOCUT) GOTO 100
      ANORM = MAX(ANORM,ABS(CC(I)))
 100  CONTINUE
*
      NTYP(INORM) = NRE
      CC(NBEG(INORM)) = ANORM
      RETURN
      END
*%%
      SUBROUTINE DANOW(INA,INW,INORM)
*     *******************************
*
*     THIS SUBROUTINE COMPUTES THE WEIGHTED NORM OF THE DA VECTOR A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XW(0:LNO)
*
      IF(NTYP(INW).NE.NRE) CALL FOXNTY(INW)
      IF(NTYP(INA).EQ.NDA) THEN
        LRC = 1
      ELSEIF(NTYP(INA).EQ.NCD) THEN
        LRC = 2
      ELSE
        CALL FOXNTY(INA)
      ENDIF
*
      XW(0) = 1.D0
      XW(1) = CC(NBEG(INW))
      DO 10 I=2,NOMAX
  10  XW(I) = XW(I-1)*XW(1)
*
      ANORM = 0.D0
      DO 100 I=NBEG(INA),NEND(INA),LRC
      IO = IEO(NC(I))
      IF(IO.GT.NOCUT) GOTO 100
      IF(LRC.EQ.1) THEN
        ANORM = MAX(ANORM,ABS(CC(I)*XW(IO)))
      ELSE
        ANORM = MAX(ANORM,SQRT(CC(I)*CC(I)+CC(I+1)*CC(I+1))*XW(IO))
      ENDIF
 100  CONTINUE
*
      NTYP(INORM) = NRE
      CC(NBEG(INORM)) = ANORM
      RETURN
      END
*
      SUBROUTINE DANORO(INA,IIV,INMN,IMN,INN)
*     ***************************************
*
*     THIS SUBROUTINE COMPUTES THE MAXIMUM NORMS OF POWER SORTED PARTS OF
*     THE DA INA. THE POWER SORTING IS PERFORMED WITH RESPECT TO X_IV.
*     IF IV=0, A MERE POWER SORTING IS PERFORMED.
*     THE RESULT IS RETURNED VIA AN ARRAY IMN WITH LENGTH NMN, AND NN+1
*     ELEMENTS WILL BE RETURNED, WHERE THE K+1 ST ELEMENT IS THE NORM OF THE
*     K-TH POWER PART IN X_IVT, AND THE HIGHEST POWER IN X_IV IS RETURNED IN NN.
*     FOR WEIGHTED DA COMPUTATION, THE "POWER" IS DIVIDED BY THE WEIGHT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CN(0:LNO)
      INTEGER JJ(LNV)
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(IIV).NE.NRE) CALL FOXNTY(IIV)
      IF(NTYP(INMN).NE.NRE) CALL FOXNTY(INMN)
      IV = NINT(CC(NBEG(IIV)))
      NMN = NINT(CC(NBEG(INMN)))
*
      IF(NTYP(IMN).GT.0) THEN
         PRINT*,'$$$ ERROR IN DANORO, 4TH ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(IV.LT.0.OR.IV.GT.NVMAX) THEN
         PRINT*,'$$$ ERROR IN DANORO, VARIABLE IS OUT OF RANGE'
         PRINT*,'    VARIABLE INDEX, BOUND = ',IV,',',NVMAX
         CALL FOXDEB
      ENDIF
*
      NN = 0
      DO 10 I=0,LNO
      CN(I+1) = 0.D0
 10   CONTINUE
*
      IF(IV.EQ.0) THEN
         IF(LEW.EQ.0) THEN
            DO 20 I=NBEG(INA),NEND(INA)
            IE = IEO(NC(I))
            NN = MAX(NN,IE)
            CN(IE+1) = MAX(CN(IE+1),ABS(CC(I)))
 20         CONTINUE
         ELSE
            DO 40 I=NBEG(INA),NEND(INA)
            IA = NC(I)
            CALL DAENCW(IE1(IA),IE2(IA),JJ)
            IE = 0
            DO 30 J=1,NVMAX
            IE = IE+JJ(J)
 30         CONTINUE
            NN = MAX(NN,IE)
            CN(IE+1) = MAX(CN(IE+1),ABS(CC(I)))
 40         CONTINUE
         ENDIF
      ELSEIF(IEW(IV).GT.NOMAX) THEN
         DO 50 I=NBEG(INA),NEND(INA)
         CN(1) = MAX(CN(1),ABS(CC(I)))
 50      CONTINUE
      ELSE
         DO 60 I=NBEG(INA),NEND(INA)
         IE = IJJ(IV,NC(I))/IEW(IV)
         NN = MAX(NN,IE)
         CN(IE+1) = MAX(CN(IE+1),ABS(CC(I)))
 60      CONTINUE
      ENDIF
*
      IF(NMN.LT.NN+1) THEN
         PRINT*,'$$$ ERROR IN DANORO, '//
     *          '3RD ARGUMENT < (HIGHEST POWER)+1 = ',NN+1
         CALL FOXDEB
      ENDIF
*
      NTYP(INN) = NRE
      CC(NBEG(INN)) = DBFLOAT(NN)
*
      DO 70 I=1,NN+1
      IF(IMN+I.EQ.INN) THEN
         PRINT*,'$$$ ERROR IN DANORO, MEMORY FOR THE 4TH ARGUMENT '//
     *              'IS SHARED WITH THE 5TH ARGUMENT.'
         CALL FOXDEB
      ENDIF
      NTYP(IMN+I) = NRE
      CC(NBEG(IMN+I)) = CN(I)
 70   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE DANORS(INA,IIV,INMN,IMN,INN)
*     ***************************************
*
*     THIS SUBROUTINE COMPUTES THE SUMMATION NORMS OF POWER SORTED PARTS OF
*     THE DA INA. THE POWER SORTING IS PERFORMED WITH RESPECT TO X_IV.
*     IF IV=0, A MERE POWER SORTING IS PERFORMED.
*     THE RESULT IS RETURNED VIA AN ARRAY IMN WITH LENGTH NMN, AND NN+1
*     ELEMENTS WILL BE RETURNED, WHERE THE K+1 ST ELEMENT IS THE NORM OF THE
*     K-TH POWER PART IN X_IVT, AND THE HIGHEST POWER IN X_IV IS RETURNED IN NN.
*     FOR WEIGHTED DA COMPUTATION, THE "POWER" IS DIVIDED BY THE WEIGHT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CN(0:LNO)
      INTEGER JJ(LNV)
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(IIV).NE.NRE) CALL FOXNTY(IIV)
      IF(NTYP(INMN).NE.NRE) CALL FOXNTY(INMN)
      IV = NINT(CC(NBEG(IIV)))
      NMN = NINT(CC(NBEG(INMN)))
*
      IF(NTYP(IMN).GT.0) THEN
         PRINT*,'$$$ ERROR IN DANORS, 4TH ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(IV.LT.0.OR.IV.GT.NVMAX) THEN
         PRINT*,'$$$ ERROR IN DANORS, VARIABLE IS OUT OF RANGE'
         PRINT*,'    VARIABLE INDEX, BOUND = ',IV,',',NVMAX
         CALL FOXDEB
      ENDIF
*
      NN = 0
      DO 10 I=0,LNO
      CN(I+1) = 0.D0
 10   CONTINUE
*
      IF(IV.EQ.0) THEN
         IF(LEW.EQ.0) THEN
            DO 20 I=NBEG(INA),NEND(INA)
            IE = IEO(NC(I))
            NN = MAX(NN,IE)
            CN(IE+1) = CN(IE+1)+ABS(CC(I))
 20         CONTINUE
         ELSE
            DO 40 I=NBEG(INA),NEND(INA)
            IA = NC(I)
            CALL DAENCW(IE1(IA),IE2(IA),JJ)
            IE = 0
            DO 30 J=1,NVMAX
            IE = IE+JJ(J)
 30         CONTINUE
            NN = MAX(NN,IE)
            CN(IE+1) = CN(IE+1)+ABS(CC(I))
 40         CONTINUE
         ENDIF
      ELSEIF(IEW(IV).GT.NOMAX) THEN
         DO 50 I=NBEG(INA),NEND(INA)
         CN(1) = CN(1)+ABS(CC(I))
 50      CONTINUE
      ELSE
         DO 60 I=NBEG(INA),NEND(INA)
         IE = IJJ(IV,NC(I))/IEW(IV)
         NN = MAX(NN,IE)
         CN(IE+1) = CN(IE+1)+ABS(CC(I))
 60      CONTINUE
      ENDIF
*
      IF(NMN.LT.NN+1) THEN
         PRINT*,'$$$ ERROR IN DANORS, '//
     *          '3RD ARGUMENT < (HIGHEST POWER)+1 = ',NN+1
         CALL FOXDEB
      ENDIF
*
      NTYP(INN) = NRE
      CC(NBEG(INN)) = DBFLOAT(NN)
*
      DO 70 I=1,NN+1
      IF(IMN+I.EQ.INN) THEN
         PRINT*,'$$$ ERROR IN DANORS, MEMORY FOR THE 4TH ARGUMENT '//
     *              'IS SHARED WITH THE 5TH ARGUMENT.'
         CALL FOXDEB
      ENDIF
      NTYP(IMN+I) = NRE
      CC(NBEG(IMN+I)) = CN(I)
 70   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE DAEST(INA,IIV,INN,INC)
*     *********************************
*
*     THIS SUBROUTINE ESTIMATES THE SIZE OF THE NN-TH ORDER TERMS OF THE DA INA
*     BY A LEAST SQUARE FIT OF THE SUMMATION NORMS OF EACH ORDER TERMS
*     WITH RESPECT TO X_IV IF IV>0.
*     THE DATA POINTS FOR A LEAST SQUARE FIT ARE LIMITED TO NON-ZERO NORMS.
*     FOR WEIGHTED DA COMPUTATION, THE "POWER" IS DIVIDED BY THE WEIGHT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION YL(LNO)
*
      IF(NTYP(IIV).NE.NRE) CALL FOXNTY(IIV)
      IF(NTYP(INN).NE.NRE) CALL FOXNTY(INN)
*
      IV = NINT(CC(NBEG(IIV)))
      IF(IV.LT.0.OR.IV.GT.NVMAX) THEN
         PRINT*,'$$$ ERROR IN DAEST, VARIABLE IS OUT OF RANGE'
         PRINT*,'    VARIABLE INDEX, BOUND = ',IV,',',NVMAX
         CALL FOXDEB
      ENDIF
*
      CALL FOXAAL(ISJ1,LNO+1,1)
      CALL FOXALL(ISR1,1,1)
      CALL FOXALL(ISR2,1,1)
*
      NTYP(ISR1) = NRE
      CC(NBEG(ISR1)) = DBFLOAT(LNO+1)
*
      CALL DANORS(INA,IIV,ISR1,ISJ1,ISR2)
      NM = NINT(CC(NBEG(ISR2)))
*
      YL(1) = 0.D0
      DO 10 I=1,NM
      YL(I) = CC(NBEG(ISJ1+I+1))
 10   CONTINUE
*
      IF(MIN(NOCUT,NM).LE.0) THEN
         CC(NBEG(INC)) = 0.D0
      ELSEIF(NOCUT.EQ.1) THEN
         CC(NBEG(INC)) = YL(1)
      ELSE
         CALL FOXAAL(ISJ2,LNO,1)
         CALL FOXALL(ISR3,1,1)
*
         IC = 0
         DO 20 I=1,NM
         IF(YL(I).GT.0.D0) THEN
            IC = IC+1
            NTYP(ISJ1+IC) = NRE
            CC(NBEG(ISJ1+IC)) = DBFLOAT(I)
            NTYP(ISJ2+IC) = NRE
            CC(NBEG(ISJ2+IC)) = LOG(YL(I))
            IF(IC.EQ.1) THEN
               I1 = I
               Y1 = YL(I)
            ENDIF
         ENDIF
 20      CONTINUE
*
         IF(IC.GE.2) THEN
            CC(NBEG(ISR1)) = DBFLOAT(IC)
            CALL LSLINE(ISJ1,ISJ2,ISR1,ISR2,ISR3)
            CC(NBEG(INC)) =
     *         EXP( CC(NBEG(ISR2))*CC(NBEG(INN))+CC(NBEG(ISR3)) )
         ELSE
            IF(I1.EQ.NINT(CC(NBEG(INN)))) THEN
               CC(NBEG(INC)) = Y1
            ELSE
               PRINT*,'$$$ WARNING IN DAEST, THERE IS ONLY ONE '//
     *                'NON-ZERO TERM.'
               PRINT*,'    NORM ',Y1,' AT ORDER',I1
               PRINT*,'    TOO INACCURATE TO MAKE AN ESTIMATE, '//
     *                'THUS RETURNING 0.D0.'
               CC(NBEG(INC)) = 0.D0
            ENDIF
         ENDIF
*
         CALL FOXDAL(ISR3,1)
         CALL FOXADA(ISJ2,LNO)
      ENDIF
*
      CALL FOXDAL(ISR2,1)
      CALL FOXDAL(ISR1,1)
      CALL FOXADA(ISJ1,LNO+1)
*
      NTYP(INC) = NRE
*
      RETURN
      END
*%%
      SUBROUTINE MTREE(IMB,INB,INC,INL,INV,INA,LTREE)
*     ***********************************************
*
*     THIS SUBROUTINE INTERFACES FROM COSY PROGRAM TO SUBROUTINE MTREEF
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----POLYNOMIAL EVALUATION -------------------------------------------------
      PARAMETER(LPS=5,LPL=16,LPM=10000)
      INTEGER IL(LPM,LPS),IV(LPM,LPS),NNP(LPS),LSW,LSWS(LPS),IPOL
      DOUBLE PRECISION CF(LPM,LPL,LPS)
      COMMON /POLDAT/ CF,IL,IV,NNP,LSW,LSWS,IPOL
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      LARR = NINT(CC(NBEG(INA)))
      IB   = NINT(CC(NBEG(INB)))
*
      IF(NTYP(INC).GT.0.OR.NTYP(INL).GT.0.OR.NTYP(INV).GT.0) THEN
         PRINT*,
     *   '$$$ ERROR IN MTREE, 3RD, 4TH OR 5TH ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(LARR.GT.LPM) THEN
         PRINT*,'$$$ ERROR IN MTREE, ARRAYS EXHAUSTED. INCREASE LPM > '
     *         ,LARR-1
         CALL FOXDEB
      ELSEIF(IB.GT.LPL) THEN
         PRINT*,'$$$ ERROR IN MTREE, ARRAY EXHAUSTED. INCREASE LPL > '
     *         ,IB-1
         CALL FOXDEB
      ENDIF
*
      CALL MTREEF(1,IMB,IB)
      NTERM = NNP(1)
*
      DO 20 I=1,NTERM
      NTYP(INL+I) = NRE
      NTYP(INV+I) = NRE
      CC(NBEG(INL+I)) = IL(I,1)
      CC(NBEG(INV+I)) = IV(I,1)
*
      IA = INC+I
      DO 10 J=1,IB
      NTYP(IA) = NRE
      CC(NBEG(IA)) = CF(I,J,1)
      IA = IA+LARR
  10  CONTINUE
  20  CONTINUE
*
      I=NTERM+1
      NTYP(INL+I) = NRE
      CC(NBEG(INL+I)) = IL(I,1)
*
      NTYP(LTREE) = NRE
      CC(NBEG(LTREE)) = NTERM
*
      RETURN
      END
*
      SUBROUTINE MTREEF(ISW,IMB,IB)
*     *****************************
*
*     THIS SUBROUTINE EFFICIENTLY DETERMINES A NEARLY OPTIMALLY SHORT TREE THAT
*     REACHES EVERY NONZERO ELEMENT IN ONE OF THE IB DA VECTORS MB.
*
*     REFER BLOCKDATA POLDT FOR THE INFORMATION ON COMMON VARIABLES IN POLDAT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----POLYNOMIAL EVALUATION -------------------------------------------------
      PARAMETER(LPS=5,LPL=16,LPM=10000)
      INTEGER IL(LPM,LPS),IV(LPM,LPS),NNP(LPS),LSW,LSWS(LPS),IPOL
      DOUBLE PRECISION CF(LPM,LPL,LPS)
      COMMON /POLDAT/ CF,IL,IV,NNP,LSW,LSWS,IPOL
*----------------------------------------------------------------------------
*-----NON-ABORTING ARITHMETIC -----------------------------------------------
      INTEGER LARI
      COMMON /ARIST/ LARI
*----------------------------------------------------------------------------
*
      INTEGER JV(0:LNO)
      COMMON /MONDA/ ID1(LNV),ID2(LNV)
      COMMON /TREEF/ NFATH(LEA)
      INTEGER JA(LNV),IEWS(LNV)
*
      IF(NTYP(IMB).GT.0) THEN
         PRINT*,'$$$ ERROR IN MTREEF, FIRST ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(ISW.LE.0.OR.ISW.GT.LPS) THEN
         PRINT*,'$$$ ERROR IN MTREEF, SWITCH ISW IS OUT OF RANGE.'
         PRINT*,'    ISW, BOUND = ',ISW,',',LPS
         CALL FOXDEB
      ENDIF
*
*     CONSISTENCY CHECKS AND PREPS
*     ****************************
*
      JMEM = IMEM+IB*NMMAX
      CALL MEMCHK(JMEM)
      DO 10 J=IMEM+1,JMEM
  10  CC(J) = 0.D0
*
      JC = IMEM
      DO 20 J=JC+2,JC+NMMAX
  20  NC(   J) = 0
      NC(JC+1) = 1
*
*     FIND ALL THE NODES IN TREE
*     **************************
*
      JMEM = IMEM
      DO 50 I=1,IB
*
      IF(NTYP(IMB+I).EQ.NDA) THEN
         NEDA = NEND(IMB+I)
C     ELSEIF(NTYP(IMB+I).EQ.NTM) THEN
C        NVTM = NC(NEND(IMB+I)-1)
C        NETM = NEND(IMB+I)
C        NEDA = NBEG(IMB+I)+NC(NETM)-1
C        IF(NVTM.LT.0) THEN
C           IF(LARI.EQ.1) THEN
C              PRINT*,'$$$ WARNING IN MTREEF, TM ARITHMETIC FAILURE'
C           ELSEIF(LARI.EQ.0) THEN
C              PRINT*,'$$$ ERROR IN MTREEF, TM ARITHMETIC FAILURE'
C              CALL FOXDEB
C           ENDIF
C        ENDIF
      ELSE
         CALL FOXNTY(IMB+I)
      ENDIF
*
      DO 40 J=NBEG(IMB+I),NEDA
      NCJ = NC(J)
      IF(IEO(NCJ).GT.NOCUT) GOTO 40
      CC(JMEM+NCJ) = CC(J)
  30  IF(NC(JC+NCJ).EQ.1) GOTO 40
      NC(JC+NCJ) = 1
      NCJ = NFATH(NCJ)
      GOTO 30
  40  CONTINUE
      JMEM = JMEM+NMMAX
  50  CONTINUE
*
*     SETTING UP TREE STRUCTURE
*     *************************
*
      DO 60 I=1,IB
      CF(1,I,ISW) = CC(IMEM+(I-1)*NMMAX+1)
  60  CONTINUE
*
      IL(1,ISW) = 0
      IV(1,ISW) = 0
*
      IF(NC(JC+1).EQ.1) THEN
         NC(JC+1) = -1
      ELSE
         PRINT*,'@@@ ERROR IN MTREEF, ZEROTH ORDER TERM OF ICHK IS ZERO'
         CALL FOXDEB
      ENDIF
*
      DO 70 JL=1,NOMAX
  70  JV(JL) = 0
*
      DO 80 I=1,NVMAX
      JA(I) = I
      IEWS(I) = IEW(I)
  80  CONTINUE
*
      DO 90 I=1,NVMAX-1
      DO 90 J=I+1,NVMAX
      IF(IEWS(I).GT.IEWS(J)) THEN
         IEWS0 = IEWS(I)
         IEWS(I) = IEWS(J)
         IEWS(J) = IEWS0
         JA0 = JA(I)
         JA(I) = JA(J)
         JA(J) = JA0
      ENDIF
  90  CONTINUE
*
      IT1 = 0
      IT2 = 0
*
      JL     = 0
      JLW    = 0
      NCHKJJ = 1
      NTERM  = 1
*
 100  CONTINUE
      IF(JL.EQ.0.AND.NCHKJJ.LE.0) GOTO 200
      IF(JLW+IEW(JA(1)).LE.NOCUT.AND.NCHKJJ.EQ.1) THEN
         JL = JL+1
         IT1 = IT1+ID1(JA(1))
         IT2 = IT2+ID2(JA(1))
         JLW = JLW+IEW(JA(1))
         JV(JL) = 1
      ELSEIF(JV(JL).EQ.NVMAX) THEN
         IT1 = IT1-ID1(JA(NVMAX))
         IT2 = IT2-ID2(JA(NVMAX))
         JLW = JLW-IEW(JA(NVMAX))
         JV(JL) = 0
         JL = JL-1
         NCHKJJ = 0
         GOTO 100
      ELSE
         IT1 = IT1-ID1(JA(JV(JL)))
         IT2 = IT2-ID2(JA(JV(JL)))
         JLW = JLW-IEW(JA(JV(JL)))
         JV(JL) = JV(JL)+1
         IF(JLW+IEW(JA(JV(JL))).LE.NOCUT) THEN
            IT1 = IT1+ID1(JA(JV(JL)))
            IT2 = IT2+ID2(JA(JV(JL)))
            JLW = JLW+IEW(JA(JV(JL)))
         ELSE
            JV(JL) = 0
            JL = JL-1
            NCHKJJ = 0
            GOTO 100
         ENDIF
      ENDIF
*
      JCN = IA1(IT1)+IA2(IT2)
      NCHKJJ = NC(JC+JCN)
*
      IF(NCHKJJ.LE.0) GOTO 100
*
      NTERM = NTERM+1
      IF(NTERM.GT.LPM) THEN
         PRINT*,'$$$ ERROR IN MTREEF, OUTPUT ARRAYS EXHAUSTED.'
     *          //' INCREASE LPM > ',NTERM-1
         CALL FOXDEB
      ENDIF
*
      NC(JC+JCN) = -1
*
      IP = IMEM+JCN
      DO 110 I=1,IB
      CF(NTERM,I,ISW) = CC(IP)
      IP = IP+NMMAX
 110  CONTINUE
*
      IL(NTERM,ISW) = JL
      IV(NTERM,ISW) = JA(JV(JL))
*
      GOTO 100
*
 200  CONTINUE
      IF(NTERM+1.GT.LPM) THEN
         PRINT*,'$$$ ERROR IN MTREEF, OUTPUT ARRAYS EXHAUSTED.'
     *          //' INCREASE LPM > ',NTERM
         CALL FOXDEB
      ENDIF
*
      IL(NTERM+1,ISW) = 1
      NNP(ISW) = NTERM
*
*     PERFORMING CROSS CHECKS
*     ***********************
*
      NTERMF = 0
      DO 300 I=1,NMMAX
      IF(NC(JC+I).EQ.-1) THEN
         NTERMF = NTERMF+1
      ELSEIF(NC(JC+I).NE.0) THEN
         PRINT*,'@@@ ERROR IN MTREEF, NOT ALL TERMS IN ICHK ARE -1'
         CALL FOXSTP(1)
      ENDIF
 300  CONTINUE
*
      IF(NTERMF.NE.NTERM) THEN
         PRINT*,'@@@ ERROR IN MTREEF, NTERMF # NTERM'
         CALL FOXSTP(1)
      ENDIF
      RETURN
      END
*%%
      SUBROUTINE INTPOL(IP,IL)
*     ************************
*
*     DETERMINES THE COEFFICIENTS OF THE ANTISYMMETRIC 2*L+1 ST ORDER
*     POLYNOMIAL P THAT SATISFIES P(1) = 1, P^(I)(1) = 0 (I=1,...L).
*     THIS IS A HAMILTON-TYPE POLYNOMIAL FOR DERIVATIVE PRESERVING
*     INTERPOLATION.
*
      PARAMETER(LSP=10)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(2*LSP+2,2*LSP+2),AI(2*LSP+2,2*LSP+2),CFAC(0:2*LSP+2)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      L = CC(NBEG(IL))
      IF(NTYP(IL).NE.NRE) THEN
         PRINT*,'$$$ ERROR IN INTPOL, 2ND ARGUMENT IS NOT REAL NUMBER'
         CALL FOXDEB
      ELSEIF(L.GT.LSP.OR.L.LE.0) THEN
         PRINT*,'$$$ ERROR IN INTPOL, ORDER IS OUT OF BOUND'
         PRINT*,'    ORDER, BOUND = ',L,',',LSP
         CALL FOXDEB
      ENDIF
*
      INP = NBEG(IP+1)-1
      IF(NTYP(IP).GT.0) THEN
         PRINT*,'$$$ ERROR IN INTPOL, FIRST ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(NMAX(IP+1).NE.NBEG(IP+1)) THEN
         PRINT*,'$$$ ERROR IN INTPOL, ARRAY IS NOT '//
     *          'ALLOCATED WITH LENGTH ONE'
         CALL FOXDEB
      ELSEIF(NDIM(NMAX(IP)).LT.2*L+2) THEN
         PRINT*,'$$$ ERROR IN INTPOL, ARRAY HAS INSUFFICIENT LENGTH'
         CALL FOXDEB
      ENDIF
*
      CFAC(0) = 1.D0
      DO 10 I=1,2*L+2
  10  CFAC(I) = CFAC(I-1)*I
*
*     SETTING UP MATRIX
*
      DO 20 I=1,2*L+2
      DO 20 K=1,2*L+2
  20  A(I,K) = 0.D0
*
      DO 30 N=0,L
      DO 30 K=1,2*L+1,2
      IF(N.GT.K) THEN
         A(N+1,K+1) = 0
      ELSE
         A(N+1,K+1) = 1.D0/CFAC(K-N)
      ENDIF
  30  CONTINUE
*
      DO 40 I=L+2,2*L+2
  40  A(I,2*(I-L-1)-1) = 1.D0
*
*     WRITE(*,'(6G12.4)') ((A(I,K),K=1,6),I=1,6)
      CALL MATINV(A,AI,2*L+2,2*LSP+2,IER)
*     WRITE(*,'(6G12.4)') ((AI(I,K),K=1,6),I=1,6)
*
      IF(IER.NE.0) THEN
         PRINT*,'@@@ ERROR IN INTPOL, MATRIX INVERSION FAILED. 1'
         CALL FOXDEB
      ENDIF
*
      DO 50 I=1,2*L+2
      DO 50 K=1,2*L+2
      S = 0
      DO 45 J=1,2*L+2
  45  S = S+A(I,J)*AI(J,K)
      IF(I.EQ.K) THEN
         IF(ABS(S-1D0).GT.2.D-5) THEN
            WRITE(6,*) 'ERROR:',I,K,S
            IER = 1
         ENDIF
      ELSE
         IF(ABS(S).GT.2.D-5) THEN
            WRITE(6,*) 'ERROR:',I,K,S
            IER = 1
         ENDIF
      ENDIF
  50  CONTINUE
*
      IF(IER.NE.0) THEN
         PRINT*,'@@@ ERROR IN INTPOL, MATRIX INVERSION FAILED. 2'
      ENDIF
*
      DO 60 I=1,2*L+2
      CC(INP+I) = AI(I,1)/CFAC(I-1)
  60  CONTINUE
*
      RETURN
      END
*%%
      SUBROUTINE LINV(JA,JAI,JN,JNMX,JIER)
*     ************************************
*
*     THIS SUBROUTINE SERVES AS AN INTERFACE BETWEEN FOXY AND THE LIBRARY
*     ROUTINE MATINV. IT CAN BE USED AS AN EXAMPLE FOR OTHER INTERFACE
*     PROBLEMS FROM FOXY TO FORTRAN.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      PARAMETER(LS=100)
      DIMENSION A(LS,LS),AI(LS,LS)
*
      IF(NTYP(JN  ).NE.NRE) CALL FOXNTY(JN  )
      IF(NTYP(JNMX).NE.NRE) CALL FOXNTY(JNMX)
*
      IF(NBEG(JA+1).NE.NMAX(JA+1).OR.NBEG(JAI+1).NE.NMAX(JAI+1)) THEN
         PRINT*,'$$$ ERROR IN LINV, MEMORY SIZE OF MATRICES MUST BE 1.'
         CALL FOXDEB
      ENDIF
*
      N   = NINT(CC(NBEG(JN  )))
      NMX = NINT(CC(NBEG(JNMX)))
*
      IF(NMX.GT.LS) THEN
         PRINT*,'$$$ ERROR IN LINV, LS EXHAUSTED. INCREASE LS >= ',NMX
         CALL FOXDEB
      ELSEIF(N.GT.NMX) THEN
         PRINT*,'$$$ ERROR IN LINV, THE 3RD ARGUMENT EXCEEDS THE 4TH'
         PRINT*,'    THE 3RD: ',N,' , THE 4TH: ',NMX
         CALL FOXDEB
      ENDIF
*
      IA  = NBEG(JA +1)-1-NMX
      IAI = NBEG(JAI+1)-1-NMX
*
      DO 10 I=1,N
      DO 10 K=1,N
      A(I,K) = CC(IA+I+K*NMX)
  10  CONTINUE
*
      CALL MATINV(A,AI,N,LS,IER)
*
      DO 20 I=1,N
      DO 20 K=1,N
      CC(IAI+I+K*NMX) = AI(I,K)
  20  CONTINUE
*
      NTYP(JIER) = NRE
      CC(NBEG(JIER)) = IER
*
      RETURN
      END
*
      SUBROUTINE LDET(JA,JN,JNMX,JDT)
*     *******************************
*
*     THIS SUBROUTINE COMPUTES THE DETERMINANT OF A MATRIX
*     INPUT A: NXN MATRIX
*           NMX: PHYSICAL DIMENSION OF A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      PARAMETER(LS=100)
      DIMENSION A(LS,LS),INDX(LS)
*
      IF(NTYP(JN  ).NE.NRE) CALL FOXNTY(JN  )
      IF(NTYP(JNMX).NE.NRE) CALL FOXNTY(JNMX)
*
      IF(NBEG(JA+1).NE.NMAX(JA+1)) THEN
         PRINT*,'$$$ ERROR IN LDET, MEMORY SIZE OF MATRIX MUST BE 1.'
         CALL FOXDEB
      ENDIF
*
      N   = NINT(CC(NBEG(JN  )))
      NMX = NINT(CC(NBEG(JNMX)))
*
      IF(NMX.GT.LS) THEN
         PRINT*,'$$$ ERROR IN LDET, LS EXHAUSTED. INCREASE LS >= ',NMX
         CALL FOXDEB
      ELSEIF(N.GT.NMX) THEN
         PRINT*,'$$$ ERROR IN LDET, THE 2ND ARGUMENT EXCEEDS THE 3RD'
         PRINT*,'    THE 2ND: ',N,' , THE 3RD: ',NMX
         CALL FOXDEB
      ENDIF
*
      IA  = NBEG(JA+1)-1-NMX
*
      DO 10 I=1,N
      DO 10 K=1,N
      A(I,K) = CC(IA+I+K*NMX)
  10  CONTINUE
*
      CALL LUDCMP(A,N,LS,INDX,DT,IER)
*
      IF(IER.EQ.132) THEN
        DT = 0.D0
      ELSE
        DO 2 I=1,N
    2   DT = DT*A(I,I)
      ENDIF
*
      CC(NBEG(JDT)) = DT
      NTYP(JDT) = NRE
*
      END
*
      SUBROUTINE MATINV(A,AI,N,NMX,IER)
*     *********************************
*
*     THIS SUBROUTINE INVERTS THE MATRIX A AND STORES THE RESULT IN AI
*     INPUT  A   - SAVED
*            N   - ORDER OF MATRIX < 50
*     OUTPUT AI  - A INVERSE
*            IER - 0 NO ERROR
*                  132 ZERO DETERMINANT
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NMAX=50,TINY=1.D-24)
      DIMENSION A(NMX,NMX),AI(NMX,NMX),AW(NMAX,NMAX),INDX(NMAX)
*
      IF(N.GE.NMAX) THEN
         PRINT*,'$$$ ERROR IN MATINV, MATRIX TOO LARGE '
         PRINT*,'    MATRIX SIZE, BOUND = ',N,',',NMAX
         CALL FOXDEB
      ENDIF
*
      DO 12 I=1,N
         DO 11 J=1,N
            AW(I,J) = A(I,J)
11       AI(I,J) = 0.D0
12    AI(I,I) = 1.D0
*
      CALL LUDCMP(AW,N,NMAX,INDX,D,IER)
*
      IF (IER .EQ. 132) RETURN
*
      DO 14 J=1,N
14    CALL LUBKSB(AW,N,NMAX,INDX,AI(1,J),NMX)
*
      RETURN
      END
*
      SUBROUTINE LUDCMP(A,N,NP,INDX,D,IER)
*     ************************************
*
*     THIS SUBROUTINE DECOMPOSES A MATRIX INTO LU FORMAT
*     INPUT A: NXN MATRIX - WILL BE OVERWRITTEN BY THE LU DECOMP.
*           NP: PHYSICAL DIMENSION OF A
*           INDX: ROW PERMUTATION VECTOR
*           D: EVEN OR ODD ROW INTERCHANGES
*           IER: 132 IF MATRIX A SINGULAR, 0 IF NOT
*
*     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 35
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NMAX=50,TINY=1.D-24)
      DIMENSION A(NP,NP),INDX(NP),VV(NMAX)
      IER=0
      D=1.D0
      DO 12 I=1,N
         AAMAX=0.D0
         DO 11 J=1,N
            IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11       CONTINUE
         IF(ABS(AAMAX).LT.TINY) THEN
            IER=132
            RETURN
         ENDIF
         VV(I)=1.D0/AAMAX
12    CONTINUE
      DO 19 J=1,N
         IF(J.GT.1) THEN
            DO 14 I=1,J-1
               SUM=A(I,J)
               IF(I.GT.1) THEN
                  DO 13 K=1,I-1
                     SUM=SUM-A(I,K)*A(K,J)
13                CONTINUE
                  A(I,J)=SUM
               ENDIF
14          CONTINUE
         ENDIF
         AAMAX=0.D0
         DO 16 I=J,N
            SUM=A(I,J)
            IF (J.GT.1) THEN
               DO 15 K=1,J-1
                  SUM=SUM-A(I,K)*A(K,J)
15             CONTINUE
               A(I,J)=SUM
            ENDIF
            DUM=VV(I)*ABS(SUM)
            IF(DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            ENDIF
16       CONTINUE
         IF (J.NE.IMAX) THEN
            DO 17 K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
17          CONTINUE
            D=-D
            VV(IMAX)=VV(J)
         ENDIF
         INDX(J)=IMAX
         IF(J.NE.N) THEN
            IF(ABS(A(J,J)).LT.TINY) THEN
               IER = 132
               RETURN
            ENDIF
            DUM=1.D0/A(J,J)
            DO 18 I=J+1,N
               A(I,J)=A(I,J)*DUM
18          CONTINUE
         ENDIF
19    CONTINUE
      IF(ABS(A(N,N)).LT.TINY) IER = 132
      RETURN
      END
*%%
      SUBROUTINE LUBKSB(A,N,NP,INDX,B,NMX)
*     ************************************
*
*     THIS SUBROUTINE SOLVES SET OF LINEAR EQUATIONS AX=B,
*     INPUT A: NXN MATRIX IN lu FORM GIVEN BY LUDCMP
*           NP: PHYSICAL DIMENSION OF A
*           INDX: ROW PERMUTATION VECTOR
*           D: EVEN OR ODD ROW INTERCHANGES
*           B: RHS OF LINEAR EQUATION - WILL BE OVERWRITTEN BY X
*
*     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 36
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NP,NP),INDX(NP),B(NMX)
      II = 0
      DO 12 I=1,N
         LL = INDX(I)
         SUM = B(LL)
         B(LL) = B(I)
         IF(II.NE.0) THEN
            DO 11 J=II,I-1
               SUM = SUM-A(I,J)*B(J)
11          CONTINUE
         ELSEIF(SUM.NE.0.D0) THEN
            II = I
         ENDIF
         B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
         SUM=B(I)
         IF(I.LT.N) THEN
            DO 13 J=I+1,N
               SUM = SUM-A(I,J)*B(J)
13          CONTINUE
         ENDIF
         B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
*%%
      SUBROUTINE LEV(IMA,IMER,IMEI,IMT,INN,INMX)
*     ******************************************
*
*     THIS SUBROUTINE COMPUTES EIGENVALUES AND EIGENVECTORS OF A MATRIX A.
*     SEE ETY2 FOR THE EIGENVECTOR ASSIGNMENT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      PARAMETER(LS=100)
      DIMENSION A(LS,LS),T(LS,LS),AA(LS,LS),EVR(LS),EVI(LS),ORT(LS)
*
      IF(NTYP(INN) .NE.NRE) CALL FOXNTY(INN)
      IF(NTYP(INMX).NE.NRE) CALL FOXNTY(INMX)
*
      IF(NTYP(IMA).GT.0.OR.NTYP(IMT).GT.0.OR.
     *   NTYP(IMER).GT.0.OR.NTYP(IMEI).GT.0) THEN
         PRINT*,
     *   '$$$ ERROR IN LEV, 1ST TO 4TH ARGUMENTS MUST BE ARRAYS.'
         CALL FOXDEB
      ELSEIF(IMT+1.EQ.IMER+1.OR.IMER+1.EQ.IMEI+1.OR.IMT+1.EQ.IMEI+1)
     *   THEN
         PRINT*,'$$$ WARNING IN LEV, '//
     *          'VARIABLES FOR THE 2ND TO 4TH ARGUMENTS MAY BE SAME.'
      ENDIF
*
      N   = NINT(CC(NBEG(INN)))
      NMX = NINT(CC(NBEG(INMX)))
      IF(NMX.GT.LS) THEN
         PRINT*,'$$$ ERROR IN LEV, LS EXHAUSTED. INCREASE LS >= ',NMX
         CALL FOXDEB
      ELSEIF(N.GT.NMX) THEN
         PRINT*,'$$$ ERROR IN LEV, THE 5TH ARGUMENT EXCEEDS THE 6TH'
         PRINT*,'    THE 5TH: ',N,' , THE 6TH: ',NMX
         CALL FOXDEB
      ENDIF
*
      DO 10 I=1,N
      JA = IMA+I
      DO 10 J=1,N
      A(I,J) = CC(NBEG(JA))
      AA(I,J) = A(I,J)
      JA = JA+NMX
  10  CONTINUE
*
      CALL ETY(LS,N,1,N,AA,ORT)
      CALL ETYT(LS,N,1,N,AA,ORT,T)
      CALL ETY2(LS,N,1,N,AA,EVR,EVI,T,INFO)
      IF (INFO.NE.0) THEN
         PRINT*,'$$$ ERROR IN LEV, EIGENVECTOR COMPUTATION UNSUCCESSFUL'
      ENDIF
*
      DO 30 I=1,N
      JT = IMT+I
      DO 20 J=1,N
      NTYP(JT) = NRE
      CC(NBEG(JT)) = T(I,J)
      JT = JT+NMX
  20  CONTINUE
      NTYP(IMER+I) = NRE
      NTYP(IMEI+I) = NRE
      CC(NBEG(IMER+I)) = EVR(I)
      CC(NBEG(IMEI+I)) = EVI(I)
  30  CONTINUE
*
      RETURN
      END
*%%
      SUBROUTINE MBLOCK(JA,JT,JTI,JNMX,JN)
*     ************************************
*
*     THIS SUBROUTINE INTERFACES BETWEEN FOXY AND THE ROUTINE MBLOC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(LS=100)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      DIMENSION A(LS,LS),T(LS,LS),TI(LS,LS)
*
      IF(NTYP(JNMX).NE.NRE) CALL FOXNTY(JNMX)
      IF(NTYP(JN).NE.NRE) CALL FOXNTY(JN)
*
      IF(NTYP(JA).GT.0.OR.NTYP(JT).GT.0.OR.NTYP(JTI).GT.0) THEN
         PRINT*,
     *   '$$$ ERROR IN MBLOCK, 1ST, 2ND OR 3RD ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(JT+1.EQ.JTI+1) THEN
         PRINT*,'$$$ WARNING IN MBLOCK, '//
     *          'VARIABLES FOR THE 2ND AND 3RD ARGUMENTS ARE SAME.'
      ELSEIF(NBEG(JA+1).NE.NMAX(JA+1).OR.
     *   NBEG(JT+1).NE.NMAX(JT+1).OR.NBEG(JTI+1).NE.NMAX(JTI+1)) THEN
         PRINT*,
     *   '$$$ ERROR IN MBLOCK, MEMORY SIZE OF MATRICES MUST BE 1.'
         CALL FOXDEB
      ENDIF
*
      N   = NINT(CC(NBEG(JN)))
      NMX = NINT(CC(NBEG(JNMX)))
*
      IF(NMX.GT.LS) THEN
         PRINT*,'$$$ ERROR IN MBLOCK, LS EXHAUSTED. INCREASE LS >= ',NMX
         CALL FOXDEB
      ELSEIF(N.GT.NMX) THEN
         PRINT*,'$$$ ERROR IN MBLOCK, THE 5TH ARGUMENT EXCEEDS THE 4TH'
         PRINT*,'    THE 5TH: ',N,' , THE 4TH: ',NMX
         CALL FOXDEB
      ENDIF
*
      IA  = NBEG(JA +1)-1-NMX
      IT  = NBEG(JT +1)-1-NMX
      ITI = NBEG(JTI+1)-1-NMX
*
      DO 10 I=1,N
      DO 10 K=1,N
      A(I,K) = CC(IA+I+K*NMX)
  10  CONTINUE
*
      CALL MBLOC(A,T,TI,LS,N)
*
      DO 20 I=1,N
      DO 20 K=1,N
      CC(IT +I+K*NMX) = T (I,K)
      CC(ITI+I+K*NMX) = TI(I,K)
  20  CONTINUE
*
      RETURN
      END
*
      SUBROUTINE MBLOC(A,T,TI,ND,N)
*     *****************************
*
*     THIS SUBROUTINE DECOUPLES A REAL N BY N MATRIX WITH DISTINCT EIGENVALUES
*     INTO 2 BY 2 BLOCKS OR 1 BY 1 BLOCKS ALONG THE DIAGONAL. THE 2 BY 2
*     BLOCKS WILL APPEAR AT THE TOP, BUT THE OTHER ARRANGEMENT IS ARBITRARY.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER(LS=100,EPS=1.D-12)
      DIMENSION A(LS,LS),T(LS,LS),TI(LS,LS)
      DIMENSION AA(LS,LS),EVR(LS),EVI(LS),ORT(LS)
*
*      CALL MATOUT(' A IN MBLOC',A,LS)
*
*     CHECK IF MATRIX IS ALREADY IN PROPER FORM
*
*      DO 20 I=1,N
*      DO 10 J=1,2*((I-1)/2)
*  10  IF(ABS(A(I,J)).GT.EPS) GOTO 40
*      DO 20 J=2*((I-1)/2)+3,N
*  20  IF(ABS(A(I,J)).GT.EPS) GOTO 40
**
*      DO 30 I=1,N
*      DO 25 J=1,N
*      T (I,J) = 0.D0
*      TI(I,J) = 0.D0
*  25  CONTINUE
*      T (I,I) = 1.D0
*      TI(I,I) = 1.D0
*  30  CONTINUE
*      RETURN
*  40  CONTINUE
*
*     PERFORM BLOCKING
*
      DO 50 I=1,N
      DO 50 J=1,N
  50  AA(I,J) = A(I,J)
*
      CALL ETY(LS,N,1,N,AA,ORT)
      CALL ETYT(LS,N,1,N,AA,ORT,T)
      CALL ETY2(LS,N,1,N,AA,EVR,EVI,T,INFO)
      IF (INFO.NE.0) THEN
         PRINT*,'$$$ ERROR IN MBLOC, '//
     *          'EIGENVECTOR COMPUTATION UNSUCCESSFUL'
      ENDIF
*
*     SORT AND INVERT EIGENVECTORS
*
      CALL EVPREP(EVR,EVI,T,LS,N)
      CALL MATINV(T,TI,N,LS,IER)
*
*      CALL MATOUT(' T IN MBLOC',T,LS)
*      CALL MATOUT(' TI IN MBLOC',TI,LS)
*
*      CHECK SYMPLECTICITY
*
*      DO 60 I=1,N
*      DO 60 J=1,N
*  60  RJ(I,J) = 0
*      DO 70 I=1,N,2
*      RJ(I,I+1) =  1.D0
*  70  RJ(I+1,I) = -1.D0
*
*      CALL MATCOP(T,SCR1,N,LS)
*      CALL MATMUL(SCR1,T,SCR2,N,LS)
*      CALL MATMUL(SCR1,RJ,SCR2,N,LS)
*      CALL MATMUL(SCR2,T,SCR1,N,LS)
*      CALL MATOUT(' MT J M ',SCR1,LS)
*
*      ISYM = 1
*      DO 80 I=1,N
*      DO 80 J=1,N
*      IF(ABS(SCR1(I,J)-RJ(I,J)).GT.1E-6) ISYM = 0
*  80  CONTINUE
*
      RETURN
      END
*
      SUBROUTINE EVPREP(EVR,EVI,EVV,ND,N)
*     ***********************************
*
*     SORTS EIGENVALUES AND EIGENVECTORS SUCH THAT ALL COMPLEX PAIRS OCCUR
*     AT THE BEGINNING AND SCALES THEM TO PROVIDE A SYMPLECTIC TRANSFORMATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(LS=100,TINY=1.D-24)
      DIMENSION EVR(ND),EVI(ND),EVV(ND,ND)
      DIMENSION AVR(LS),AVI(LS),AVV(LS,LS)
*
      IC = 0
*
*     PICKING COMPLEX PAIRS FIRST
*
      DO 20 I=1,N
      IF(ABS(EVI(I)).GE.1E-10) THEN
         IC = IC+1
         AVR(IC) = EVR(I)
         AVI(IC) = EVI(I)
         DO 10 J=1,N
  10     AVV(J,IC) = EVV(J,I)
      ENDIF
  20  CONTINUE
*
*     PICKING REAL PAIRS SECOND
*
      DO 40 I=1,N
      IF(ABS(EVI(I)).LT.1.D-10) THEN
         IC = IC+1
         AVR(IC) = EVR(I)
         AVI(IC) = EVI(I)
         DO 30 J=1,N
  30     AVV(J,IC) = EVV(J,I)
      ENDIF
  40  CONTINUE
*
      DO 50 I=1,N
      EVR(I) = AVR(I)
      EVI(I) = AVI(I)
      DO 50 J=1,N
  50  EVV(I,J) = AVV(I,J)
*
*     SCALING EIGENVECTOR PAIRS
*
      DO 100 J=1,N,2
      FAC = 0.D0
      DO 60 I=1,N,2
  60  FAC = FAC+EVV(I,J)*EVV(I+1,J+1)-EVV(I,J+1)*EVV(I+1,J)
      IF(ABS(FAC).LT.1.D-10) FAC = 1.D-10
      DO 70 I=1,N
  70  EVV(I,J+1) = EVV(I,J+1)/FAC
 100  CONTINUE
*
      RETURN
      END
*
C     SUBROUTINE MATMUL(A,B,C,N,LS)
*     *****************************
*
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     DIMENSION A(LS,*), B(LS,*), C(LS,*)
*
C     DO 10 I=1,N
C     DO 10 J=1,N
C     C(I,J) = 0.D0
C     DO 10 K=1,N
C 10  C(I,J) = C(I,J)+A(I,K)*B(K,J)
*
C     RETURN
C     END
*
C     SUBROUTINE MATCOP(A,B,N,LS)
*     ***************************
*
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     DIMENSION A(LS,*),B(LS,*)
*
C     DO 10 I=1,N
C     DO 10 J=1,N
C 10  B(I,J) = A(J,I)
*
C     RETURN
C     END
*
C     SUBROUTINE MATOUT(AC,A,N)
*     *************************
*
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     DIMENSION A(N,N)
C     CHARACTER AC*(*)
*
*     WRITE(7,'(A/6(6E12.4/))') AC,((A(I,J),J=1,6),I=1,6)
*
C     RETURN
C     END
*
      SUBROUTINE ETY(NM,N,LOW,IGH,A,ORT)
*     **********************************
*
*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTHES,
*     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
*
*     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
*     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
*     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
*     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
*
*     ON INPUT-
*
*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
*          DIMENSION STATEMENT,
*
*        N IS THE ORDER OF THE MATRIX,
*
*        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
*          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
*          SET LOW=1, IGH=N,
*
*        A CONTAINS THE INPUT MATRIX.
*
*     ON OUTPUT-
*
*        A CONTAINS THE HESSENBERG MATRIX.  INFORMATION ABOUT
*          THE ORTHOGONAL TRANSFORMATIONS USED IN THE REDUCTION
*          IS STORED IN THE REMAINING TRIANGLE UNDER THE
*          HESSENBERG MATRIX,
*
*        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
*          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
*
*     FORTRAN ROUTINE BY B. S. GARBOW
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(NM,N),ORT(IGH)
*
      LA = IGH-1
      KP1 = LOW+1
      IF(LA.LT.KP1) GOTO 200
*
      DO 180 M=KP1,LA
         H = 0.D0
         ORT(M) = 0.D0
         SCALE = 0.D0
*     SCALE COLUMN (ALGOL TOL THEN NOT NEEDED)
         DO 90 I=M,IGH
   90    SCALE = SCALE+ABS(A(I,M-1))
*
         IF(SCALE.EQ.0.D0) GOTO 180
         MP = M+IGH
*     FOR I=IGH STEP -1 UNTIL M DO --
         DO 100 II=M,IGH
            I = MP-II
            ORT(I) = A(I,M-1)/SCALE
            H = H+ORT(I)*ORT(I)
  100    CONTINUE
*
         G = -SIGN(SQRT(H),ORT(M))
         H = H-ORT(M)*G
         ORT(M) = ORT(M)-G
*     FORM (I-(U*UT)/H)*A
         DO 130 J=M,N
            F = 0.D0
            DO 110 II=M,IGH
               I = MP-II
               F = F+ORT(I)*A(I,J)
  110       CONTINUE
*
            F = F/H
*
            DO 120 I=M,IGH
  120       A(I,J) = A(I,J)-F*ORT(I)
*
  130    CONTINUE
*     FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
         DO 160 I=1,IGH
            F = 0.D0
            DO 140 JJ=M,IGH
               J = MP-JJ
               F = F+ORT(J)*A(I,J)
  140       CONTINUE
*
            F = F/H
*
            DO 150 J=M,IGH
  150       A(I,J) = A(I,J)-F*ORT(J)
*
  160    CONTINUE
*
         ORT(M) = SCALE*ORT(M)
         A(M,M-1) = SCALE*G
  180 CONTINUE
*
  200 RETURN
      END
*
      SUBROUTINE ETYT(NM,N,LOW,IGH,A,ORT,Z)
*     *************************************
*
*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTRANS,
*     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
*
*     THIS SUBROUTINE ACCUMULATES THE ORTHOGONAL SIMILARITY
*     TRANSFORMATIONS USED IN THE REDUCTION OF A REAL GENERAL
*     MATRIX TO UPPER HESSENBERG FORM BY  ETY.
*
*     ON INPUT-
*
*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
*          DIMENSION STATEMENT,
*
*        N IS THE ORDER OF THE MATRIX,
*
*        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
*          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
*          SET LOW=1, IGH=N,
*
*        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
*          FORMATIONS USED IN THE REDUCTION BY  ORTHES
*          IN ITS STRICT LOWER TRIANGLE,
*
*          ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANS-
*          FORMATIONS USED IN THE REDUCTION BY  ETY.
*          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
*
*     ON OUTPUT-
*
*        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
*          REDUCTION BY  ETY,
*
*        ORT HAS BEEN ALTERED.
*
*     FORTRAN ROUTINE BY B. S. GARBOW.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(NM,IGH),ORT(IGH),Z(NM,N)
*
*     INITIALIZE Z TO IDENTITY MATRIX
      DO 80 I=1,N
      DO 60 J=1,N
   60 Z(I,J) = 0.D0
      Z(I,I) = 1.D0
   80 CONTINUE
*
      KL = IGH-LOW-1
      IF(KL.LT.1) GOTO 200
      DO 140 MM=1,KL
         MP = IGH-MM
         IF(A(MP,MP-1).EQ.0.D0) GOTO 140
         MP1 = MP+1
*
         DO 100 I=MP1,IGH
  100    ORT(I) = A(I,MP-1)
*
         DO 130 J=MP,IGH
            G = 0.D0
*
            DO 110 I=MP,IGH
  110       G = G+ORT(I)*Z(I,J)
*     DIVISOR BELOW IS NEGATIVE OF H FORMED IN ORTHES.
*     DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW
            G = (G/ORT(MP))/A(MP,MP-1)
*
            DO 120 I=MP,IGH
  120       Z(I,J) = Z(I,J)+G*ORT(I)
*
  130    CONTINUE
*
  140 CONTINUE
*
  200 RETURN
      END
*
      SUBROUTINE ETY2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
*     ********************************************
*
*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR2,
*     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
*
*     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
*     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
*     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND
*     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
*     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
*     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
*
*     ON INPUT-
*
*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
*          DIMENSION STATEMENT,
*
*        N IS THE ORDER OF THE MATRIX,
*
*        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
*          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
*          SET LOW=1, IGH=N,
*
*        H CONTAINS THE UPPER HESSENBERG MATRIX,
*
*        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
*          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
*          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
*          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
*          IDENTITY MATRIX.
*
*     ON OUTPUT-
*
*        H HAS BEEN DESTROYED,
*
*        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
*          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
*          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
*          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
*          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
*          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
*          FOR INDICES IERR+1,...,N,
*
*        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
*          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
*          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
*          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
*          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
*          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
*          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND,
*
*        IERR IS SET TO
*          ZERO       FOR NORMAL RETURN,
*          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
*                     DETERMINED AFTER 30 ITERATIONS.
*
*     ARITHMETIC IS DOUBLE PRECISION. COMPLEX DIVISION
*     IS SIMULATED BY ROUTIN ETDIV.
*
*     FORTRAN ROUTINE BY B. S. GARBOW.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER EN,ENM2
      DOUBLE PRECISION H(NM,N),WR(N),WI(N),Z(NM,N),NORM,MACHEP
      LOGICAL NOTLAS
*
*     MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
*     THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
      MACHEP = 1.D-17
*     MACHEP = R1MACH(4)
*
      IERR = 0
      NORM = 0.D0
      K = 1
*     STORE ROOTS ISOLATED BY BALANC AND COMPUTE MATRIX NORM
      DO 50 I=1,N
*
         DO 40 J=K,N
   40    NORM = NORM+ABS(H(I,J))
*
         K = I
         IF(I.GE.LOW.AND.I.LE.IGH) GOTO 50
         WR(I) = H(I,I)
         WI(I) = 0.D0
   50 CONTINUE
*
      EN = IGH
      T = 0.D0
*     SEARCH FOR NEXT EIGENVALUES
   60 IF(EN.LT.LOW) GOTO 340
      ITS = 0
      NA = EN-1
      ENM2 = NA-1
*     LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
   70 DO 80 LL=LOW,EN
         L = EN+LOW-LL
         IF(L.EQ.LOW) GOTO 100
         S = ABS(H(L-1,L-1))+ABS(H(L,L))
         IF(S.EQ.0.D0) S = NORM
         IF(ABS(H(L,L-1)).LE.MACHEP*S) GOTO 100
   80 CONTINUE
*     FORM SHIFT
  100 X = H(EN,EN)
      IF(L.EQ.EN) GOTO 270
      Y = H(NA,NA)
      W = H(EN,NA)*H(NA,EN)
      IF(L.EQ.NA) GOTO 280
      IF(ITS.EQ.30) GOTO 1000
      IF(ITS.NE.10.AND.ITS.NE.20) GOTO 130
*     FORM EXCEPTIONAL SHIFT
      T = T+X
*
      DO 120 I=LOW,EN
  120 H(I,I) = H(I,I)-X
*
      S = ABS(H(EN,NA))+ABS(H(NA,ENM2))
      X = 0.75D0*S
      Y = X
      W = -0.4375D0*S*S
  130 ITS = ITS+1
*     LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS.
      DO 140 MM=L,ENM2
         M = ENM2+L-MM
         ZZ = H(M,M)
         R = X-ZZ
         S = Y-ZZ
         P = (R*S-W)/H(M+1,M)+H(M,M+1)
         Q = H(M+1,M+1)-ZZ-R-S
         R = H(M+2,M+1)
         S = ABS(P)+ABS(Q)+ABS(R)
         P = P/S
         Q = Q/S
         R = R/S
         IF(M.EQ.L) GOTO 150
         IF(ABS(H(M,M-1))*(ABS(Q)+ABS(R)).LE.
     *      MACHEP*ABS(P)*(ABS(H(M-1,M-1))+ABS(ZZ)+ABS(H(M+1,M+1))))
     *      GOTO 150
  140 CONTINUE
*
  150 MP2 = M+2
*
      DO 160 I=MP2,EN
         H(I,I-2) = 0.D0
         IF(I.EQ.MP2) GOTO 160
         H(I,I-3) = 0.D0
  160 CONTINUE
*     DOUBLE QR STEP INVOLVING ROWS L TO EN AND COLUMNS M TO EN
      DO 260 K=M,NA
         NOTLAS = K .NE. NA
         IF(K.EQ.M) GOTO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.D0
         IF(NOTLAS) R = H(K+2,K-1)
         X = ABS(P)+ABS(Q)+ABS(R)
         IF(X.EQ.0.D0) GOTO 260
         P = P/X
         Q = Q/X
         R = R/X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF(K.EQ.M) GOTO 180
         H(K,K-1) = -S*X
         GOTO 190
  180    IF(L.NE.M) H(K,K-1) = -H(K,K-1)
  190    P = P+S
         X = P/S
         Y = Q/S
         ZZ = R/S
         Q = Q/P
         R = R/P
*     ROW MODIFICATION
         DO 210 J=K,N
            P = H(K,J)+Q*H(K+1,J)
            IF(.NOT.NOTLAS) GOTO 200
            P = P+R*H(K+2,J)
            H(K+2,J) = H(K+2,J)-P*ZZ
  200       H(K+1,J) = H(K+1,J)-P*Y
            H(K,J) = H(K,J)-P*X
  210    CONTINUE
*
         J = MIN0(EN,K+3)
*     COLUMN MODIFICATION
         DO 230 I=1,J
            P = X*H(I,K)+Y*H(I,K+1)
            IF(.NOT.NOTLAS) GOTO 220
            P = P+ZZ*H(I,K+2)
            H(I,K+2) = H(I,K+2)-P*R
  220       H(I,K+1) = H(I,K+1)-P*Q
            H(I,K) = H(I,K)-P
  230    CONTINUE
*     ACCUMULATE TRANSFORMATIONS
         DO 250 I=LOW,IGH
            P = X*Z(I,K)+Y*Z(I,K+1)
            IF(.NOT.NOTLAS) GOTO 240
            P = P+ZZ*Z(I,K+2)
            Z(I,K+2) = Z(I,K+2)-P*R
  240       Z(I,K+1) = Z(I,K+1)-P*Q
            Z(I,K) = Z(I,K)-P
  250    CONTINUE
*
  260 CONTINUE
*
      GOTO 70
*     ONE ROOT FOUND
  270 H(EN,EN) = X+T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.D0
      EN = NA
      GOTO 60
*     TWO ROOTS FOUND
  280 P = (Y-X)/2.D0
      Q = P*P+W
      ZZ = SQRT(ABS(Q))
      H(EN,EN) = X+T
      X = H(EN,EN)
      H(NA,NA) = Y+T
      IF(Q.LT.0.D0) GOTO 320
*     REAL PAIR
      ZZ = P+SIGN(ZZ,P)
      WR(NA) = X+ZZ
      WR(EN) = WR(NA)
      IF(ZZ.NE.0.D0) WR(EN) = X-W/ZZ
      WI(NA) = 0.D0
      WI(EN) = 0.D0
      X = H(EN,NA)
      S = ABS(X)+ABS(ZZ)
      P = X/S
      Q = ZZ/S
      R = SQRT(P*P+Q*Q)
      P = P/R
      Q = Q/R
*     ROW MODIFICATION
      DO 290 J=NA,N
         ZZ = H(NA,J)
         H(NA,J) = Q*ZZ+P*H(EN,J)
         H(EN,J) = Q*H(EN,J)-P*ZZ
  290 CONTINUE
*     COLUMN MODIFICATION
      DO 300 I=1,EN
         ZZ = H(I,NA)
         H(I,NA) = Q*ZZ+P*H(I,EN)
         H(I,EN) = Q*H(I,EN)-P*ZZ
  300 CONTINUE
*     ACCUMULATE TRANSFORMATIONS
      DO 310 I=LOW,IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q*ZZ+P*Z(I,EN)
         Z(I,EN) = Q*Z(I,EN)-P*ZZ
  310 CONTINUE
*
      GOTO 330
*     COMPLEX PAIR
  320 WR(NA) = X+P
      WR(EN) = X+P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GOTO 60
*     ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND VECTORS OF UPPER TRIANGULAR FORM
  340 IF(NORM.EQ.0.D0) GOTO 1001
      DO 800 NN=1,N
         EN = N+1-NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN-1
         IF(Q) 710, 600, 800
*     REAL VECTOR
  600    M = EN
         H(EN,EN) = 1.D0
         IF(NA.EQ.0) GOTO 800
         DO 700 II=1,NA
            I = EN-II
            W = H(I,I)-P
            R = H(I,EN)
            IF(M.GT.NA) GOTO 620
*
            DO 610 J=M,NA
  610       R = R+H(I,J)*H(J,EN)
*
  620       IF(WI(I).GE.0.D0) GOTO 630
            ZZ = W
            S = R
            GOTO 700
  630       M = I
            IF(WI(I).NE.0.D0) GOTO 640
            T = W
            IF(W.EQ.0.D0) T = MACHEP*NORM
            H(I,EN) = -R/T
            GOTO 700
*     SOLVE REAL EQUATIONS
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
            T = (X*S-ZZ*R)/Q
            H(I,EN) = T
            IF(ABS(X).LE.ABS(ZZ)) GOTO 650
            H(I+1,EN) = (-R-W*T)/X
            GOTO 700
  650       H(I+1,EN) = (-S-Y*T)/ZZ
  700    CONTINUE
*     END REAL VECTOR
         GOTO 800
*     COMPLEX VECTOR
  710    M = NA
*     LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT EIGENVECTOR MATRIX IS
*     TRIANGULAR
         IF(ABS(H(EN,NA)).LE.ABS(H(NA,EN))) GOTO 720
         H(NA,NA) = Q/H(EN,NA)
         H(NA,EN) = -(H(EN,EN)-P)/H(EN,NA)
         GOTO 730
* 720    Z3 = CMPLX(0.D0,-H(NA,EN))/CMPLX(H(NA,NA)-P,Q)
*        H(NA,NA) = REAL(Z3)
*        H(NA,EN) = AIMAG(Z3)
  720    CALL ETDIV(Z3R,Z3I,0.D0,-H(NA,EN),H(NA,NA)-P,Q)
         H(NA,NA) = Z3R
         H(NA,EN) = Z3I
  730    H(EN,NA) = 0.D0
         H(EN,EN) = 1.D0
         ENM2 = NA-1
         IF(ENM2.EQ.0) GOTO 800
         DO 790 II=1,ENM2
            I = NA-II
            W = H(I,I)-P
            RA = 0.D0
            SA = H(I,EN)
*
            DO 760 J=M,NA
               RA = RA+H(I,J)*H(J,NA)
               SA = SA+H(I,J)*H(J,EN)
  760       CONTINUE
*
            IF(WI(I).GE.0.D0) GOTO 770
            ZZ = W
            R = RA
            S = SA
            GOTO 790
  770       M = I
            IF(WI(I).NE.0.D0) GOTO 780
*           Z3 = CMPLX(-RA,-SA)/CMPLX(W,Q)
*           H(I,NA) = REAL(Z3)
*           H(I,EN) = AIMAG(Z3)
            CALL ETDIV(Z3R,Z3I,-RA,-SA,W,Q)
            H(I,NA) = Z3R
            H(I,EN) = Z3I
            GOTO 790
*     SOLVE COMPLEX EQUATIONS
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
            VI = (WR(I)-P)*2.0*Q
            IF(VR.EQ.0.D0.AND.VI.EQ.0.D0)
     *         VR = MACHEP*NORM*(ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(ZZ))
*           Z3 = CMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA)/CMPLX(VR,VI)
*           H(I,NA) = REAL(Z3)
*           H(I,EN) = AIMAG(Z3)
            CALL ETDIV(Z3R,Z3I,X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI)
            H(I,NA) = Z3R
            H(I,EN) = Z3I
            IF(ABS(X).LE.ABS(ZZ)+ABS(Q)) GOTO 785
            H(I+1,NA) = (-RA-W*H(I,NA)+Q*H(I,EN))/X
            H(I+1,EN) = (-SA-W*H(I,EN)-Q*H(I,NA))/X
            GOTO 790
* 785       Z3 = CMPLX(-R-Y*H(I,NA),-S-Y*H(I,EN))/CMPLX(ZZ,Q)
*           H(I+1,NA) = REAL(Z3)
*           H(I+1,EN) = AIMAG(Z3)
  785       CALL ETDIV(Z3R,Z3I,-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q)
            H(I+1,NA) = Z3R
            H(I+1,EN) = Z3I
  790    CONTINUE
*     END COMPLEX VECTOR
  800 CONTINUE
*     END BACK SUBSTITUTION
*     VECTORS OF ISOLATED ROOTS
      DO 840 I=1,N
         IF(I.GE.LOW.AND.I.LE.IGH) GOTO 840
         DO 820 J=I,N
  820    Z(I,J) = H(I,J)
  840 CONTINUE
*     MULTIPLY BY TRANSFORMATION MATRIX TO GIVE VECTORS OF ORIGINAL FULL MATRIX
*
      DO 880 JJ=LOW,N
         J = N+LOW-JJ
         M = MIN0(J,IGH)
         DO 880 I=LOW,IGH
            ZZ = 0.D0
            DO 860 K=LOW,M
  860       ZZ = ZZ+Z(I,K)*H(K,J)
            Z(I,J) = ZZ
  880 CONTINUE
*
      GOTO 1001
*     SET ERROR -- NO CONVERGENCE TO AN EIGENVALUE AFTER 30 ITERATIONS
 1000 IERR = EN
 1001 RETURN
      END
*
      SUBROUTINE ETDIV(A,B,C,D,E,F)
*     *****************************
*
*     HIGH-PRECISION DIVISION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER FLIP
      FLIP = 0
      CC = C
      DD = D
      EE = E
      FF = F
      IF(ABS(F).GE.ABS(E)) THEN
        EE = F
        FF = E
        CC = D
        DD = C
        FLIP = 1
      ENDIF
      S = 1.D0/EE
      T = 1.D0/(EE+FF*(FF*S))
      IF(ABS(FF).GE.ABS(S)) THEN
        TEMP = FF
        FF = S
        S = TEMP
      ENDIF
      IF(ABS(DD).GE.ABS(S)) THEN
        A = T*(CC+S*(DD*FF))
      ELSEIF(ABS(DD).GE.ABS(FF)) THEN
        A = T*(CC+DD*(S*FF))
      ELSE
        A = T*(CC+FF*(S*DD))
      ENDIF
      IF(ABS(CC).GE.ABS(S)) THEN
        B = T*(DD-S*(CC*FF))
      ELSEIF(ABS(CC).GE.ABS(FF)) THEN
        B = T*(DD-CC*(S*FF))
      ELSE
        B = T*(DD-FF*(S*CC))
      ENDIF
      IF(FLIP.NE.0) B = -B
      RETURN
      END
*
      SUBROUTINE LSLINE(IMX,IMY,INN,INCA,INCB)
*     ****************************************
*
*     THIS SUBROUTINE COMPUTES THE LEAST SQUARE FIT LINE Y=A*X+B
*     FOR N PAIRS OF (X(),Y()).
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INN).NE.NRE) CALL FOXNTY(INN)
      N = NINT(CC(NBEG(INN)))
*
      IF(NTYP(IMX).GT.0.OR.NTYP(IMY).GT.0) THEN
         PRINT*,'$$$ ERROR IN LSLINE, 1ST OR 2ND ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(N.LE.0) THEN
         PRINT*,'$$$ ERROR IN LSLINE, NO DATA. (NUMBER OF PAIRS) < 1.'
         CALL FOXDEB
      ELSEIF(N.EQ.1) THEN
         PRINT*,'$$$ ERROR IN LSLINE, THERE IS ONLY ONE DATA POINT.'
         CALL FOXDEB
      ENDIF
*
      SX  = 0.D0
      SX2 = 0.D0
      SY  = 0.D0
      SXY = 0.D0
*
      DO 10 I=1,N
      IF(NTYP(IMX+I).NE.NRE) CALL FOXNTY(IMX+I)
      IF(NTYP(IMY+I).NE.NRE) CALL FOXNTY(IMY+I)
      X = CC(NBEG(IMX+I))
      Y = CC(NBEG(IMY+I))
      SX  = SX +X
      SX2 = SX2+X*X
      SY  = SY +Y
      SXY = SXY+X*Y
 10   CONTINUE
*
      DEN = N*SX2-SX*SX
      IF(DEN.EQ.0.D0) THEN
         PRINT*,'$$$ ERROR IN LSLINE, BAD DATA. NUMBER OF DATA POINTS ='
     *         ,N
         DO 20 I=1,N
         PRINT*,I,CC(NBEG(IMX+I)),CC(NBEG(IMY+I))
 20      CONTINUE
         CALL FOXDEB
      ENDIF
*
      CC(NBEG(INCA)) = (N*SXY-SX*SY)/DEN
      CC(NBEG(INCB)) = (SX2*SY-SX*SXY)/DEN
      NTYP(INCA) = NRE
      NTYP(INCB) = NRE
*
      RETURN
      END
*
      SUBROUTINE DAPEP(INA,INB,INV,INC)
*     *********************************
*
*     THIS SUBROUTINE PROVIDES AN INTERFACE TO DAPKP
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION JJ(LNV)
*
      IF(NTYP(INV).NE.NRE) CALL FOXNTY(INV)
      NV  = NINT(CC(NBEG(INV)))
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      CCC = CC(NBEG(INB))
*
      DO 10 J=1,LNV
  10  JJ(J) = 0
*
  20  CONTINUE
      CCC = (CCC/10.D0)
      CNEW = INT(CCC+1.D-2)
      J = NINT(10.D0*(CCC-CNEW))
      IF(J.NE.0) JJ(J) = JJ(J)+1
      CCC = CNEW
      IF(CCC.GT..5) GOTO 20
*
      CALL DAPKP(INA,JJ,NV,INC)
*
      RETURN
      END
*
      SUBROUTINE DAPKP(INA,JJ,NV,INC)
*     *******************************
*
*     THIS SUBROUTINE DETERMINES THE SUB DA VECTOR CONTAINING ALL TERMS WHOSE
*     FIRST NV EXPONENTS ARE THOSE IN JJ AND STORES IT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      INTEGER JJ(LNV),JJF(LNV)
*
      IF(NTYP(INA).EQ.NDA) THEN
         ITA = 1
      ELSEIF(NTYP(INA).EQ.NCD) THEN
         ITA = 2
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      DO 10 J=NV+1,NVMAX
  10  JJ(J) = 0
*
      CALL DADEC(JJ,ID1,ID2)
      IF(ID1.LT.0.OR.ID2.LT.0) THEN
         NTYP(INC) = NTYP(INA)
         NEND(INC) = NBEG(INC)-1
         RETURN
      ELSE
         DO 20 J=1,NV
         IF(JJ(J).GT.NOCUT) THEN
            NTYP(INC) = NTYP(INA)
            NEND(INC) = NBEG(INC)-1
            RETURN
         ENDIF
  20     CONTINUE
      ENDIF
*
      IC = NBEG(INC)-ITA
*
      DO 100 IA=NBEG(INA),NEND(INA),ITA
      NCIA = NC(IA)
      IF(IEO(NCIA).GT.NOCUT) GOTO 100
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJF)
      DO 90 J=1,NV
      IF (JJ(J).NE.JJF(J)) GOTO 100
  90  CONTINUE
      IC = IC+ITA
      CC(IC) = CC(IA)
      NC(IC) = IA1(IE1(NCIA)-ID1)+IA2(IE2(NCIA)-ID2)
      IF(ITA.EQ.2) THEN
         CC(IC+1) = CC(IA+1)
         NC(IC+1) = 0
      ENDIF
 100  CONTINUE
*
      NTYP(INC) = NTYP(INA)
      NEND(INC) = IC+ITA-1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAPKP')
*
      RETURN
      END
*
      SUBROUTINE DADIU(IIV,INA,INC)
*     *****************************
*
*     THIS SUBROUTINE DIVIDES BY THE UNIT IV AND STORES THE RESULT IN C.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      COMMON /MONDA/ ID1(LNV),ID2(LNV)
*
      IF(NTYP(IIV).NE.NRE) CALL FOXNTY(IIV)
      IF(NTYP(INA).EQ.NRE) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
         RETURN
      ELSEIF(NTYP(INA).EQ.NDA) THEN
         ITA = 1
      ELSEIF(NTYP(INA).EQ.NCD) THEN
         ITA = 2
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      IV = NINT(CC(NBEG(IIV)))
      IF(IEW(IV).GT.NOMAX) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
         RETURN
      ENDIF
*
      IC = NBEG(INC)-ITA
*
      DO 100 IA=NBEG(INA),NEND(INA),ITA
      NCIA = NC(IA)
      IFAC = IJJ(IV,NCIA)
      IF(IFAC.EQ.0) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
         RETURN
      ENDIF
      IC = IC+ITA
      CC(IC) = CC(IA)
      NC(IC) = IA1(IE1(NCIA)-ID1(IV))+IA2(IE2(NCIA)-ID2(IV))
      IF(ITA.EQ.2) THEN
        CC(IC+1) = CC(IA+1)
        NC(IC+1) = 0
      ENDIF
 100  CONTINUE
*
      NTYP(INC) = NTYP(INA)
      NEND(INC) = IC+ITA-1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DADIU')
*
      RETURN
      END
*
      SUBROUTINE DADMU(IIVD,IIVM,INA,INC)
*     ***********************************
*
*     THIS SUBROUTINE DIVIDES BY THE UNIT IVD THEN MULTIPLIES BY
*     THE UNIT IVM AND STORES THE RESULT IN C.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      COMMON /MONDA/ ID1(LNV),ID2(LNV)
*
      IF(NTYP(IIVD).NE.NRE) CALL FOXNTY(IIVD)
      IF(NTYP(IIVM).NE.NRE) CALL FOXNTY(IIVM)
      IF(NTYP(INA).EQ.NRE) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
         RETURN
      ELSEIF(NTYP(INA).EQ.NDA) THEN
         ITA = 1
      ELSEIF(NTYP(INA).EQ.NCD) THEN
         ITA = 2
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      IVD = NINT(CC(NBEG(IIVD)))
      IVM = NINT(CC(NBEG(IIVM)))
*
      IF(IEW(IVD).NE.IEW(IVM)) THEN
         PRINT*,'$$$ ERROR IN DADMU, THE DA WEIGHTS OF DIVIDING AND'
         PRINT*,'    MULTIPLYING VARIABLES HAVE TO AGREE'
         CALL FOXDEB
      ENDIF
*
      IF(IEW(IVD).GT.NOMAX) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
         RETURN
      ENDIF
*
      IC = NBEG(INC)-ITA
*
      DO 100 IA=NBEG(INA),NEND(INA),ITA
      NCIA = NC(IA)
      IFAC = IJJ(IVD,NCIA)
      IF(IFAC.EQ.0) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
         RETURN
      ENDIF
      IC = IC+ITA
      CC(IC) = CC(IA)
      NC(IC) = IA1(IE1(NCIA)+ID1(IVM)-ID1(IVD))+
     *         IA2(IE2(NCIA)+ID2(IVM)-ID2(IVD))
      IF(ITA.EQ.2) THEN
        CC(IC+1) = CC(IA+1)
        NC(IC+1) = 0
      ENDIF
 100  CONTINUE
*
      NTYP(INC) = NTYP(INA)
      NEND(INC) = IC+ITA-1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DADMU')
*
      RETURN
      END
*
      SUBROUTINE DADER(IIV,INA,INC)
*     *****************************
*
*     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
*     OF THE VECTOR A AND STORES THE RESULT IN C.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      COMMON /MONDA/ ID1(LNV),ID2(LNV)
*
      IF(NTYP(IIV).NE.NRE) CALL FOXNTY(IIV)
      IV = ABS(NINT(CC(NBEG(IIV))))
*
      IF(IV.LT.1) THEN
         PRINT*,'$$$ ERROR IN DADER, VARIABLE NUMBER IS 0'
         CALL FOXDEB
      ENDIF
*
      NTY = NTYP(INA)
*
      IW = IEW(IV)
*
      IF(NTY.EQ.NRE.OR.NTY.EQ.NCM.OR.IW.GT.NOMAX.OR.IV.GT.NVMAX) THEN
         NTYP(INC) = NDA
         IF(NTY.EQ.NCM.OR.NTY.EQ.NCD) NTYP(INC) = NCD
         NEND(INC) = NBEG(INC)-1
         CC(NBEG(INC)) = 0.D0
         RETURN
      ELSEIF(NTY.EQ.NDA) THEN
         ITA = 1
      ELSEIF(NTY.EQ.NCD) THEN
         ITA = 2
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      IC = NBEG(INC)-ITA
*
      DO 100 IA=NBEG(INA),NEND(INA),ITA
      NCIA = NC(IA)
      IFAC = IJJ(IV,NCIA)/IW
      IF(IFAC.EQ.0) GOTO 100
      IC = IC+ITA
      CC(IC) = CC(IA)*IFAC
      NC(IC) = IA1(IE1(NCIA)-ID1(IV))+IA2(IE2(NCIA)-ID2(IV))
      IF(ITA.EQ.2) THEN
        CC(IC+1) = CC(IA+1)*IFAC
        NC(IC+1) = 0
      ENDIF
 100  CONTINUE
*
      NTYP(INC) = NTY
      NEND(INC) = IC+ITA-1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DADER')
*
      RETURN
      END
*
      SUBROUTINE DAINT(IIV,INA,INC)
*     *****************************
*
*     THIS SUBROUTINE COMPUTES THE INTEGRAL WITH RESPECT TO VARIABLE I
*     OF THE VECTOR A AND STORES THE RESULT IN C.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      COMMON /MONDA/ ID1(LNV),ID2(LNV)
      INTEGER IMDA,LINT
      DOUBLE PRECISION CEX1(LNV,LNO+1),CEX2(LNV,LNO+1),A1(LNV),A2(LNV)
      COMMON /COMINT/ CEX1,CEX2,A1,A2,COBN1,COBN2,IMDA,LINT
      INTEGER JJ(LNV)
*
      IF(NTYP(IIV).NE.NRE) CALL FOXNTY(IIV)
      IV = ABS(NINT(CC(NBEG(IIV))))
*
      IF(IV.LT.1) THEN
         PRINT*,'$$$ ERROR IN DAINT, VARIABLE NUMBER IS 0'
         CALL FOXDEB
      ENDIF
*
      IW = IEW(IV)
*
      IF(NTYP(INA).EQ.NRE) THEN
         NC(NBEG(INC)) = IA1(ID1(IV))+IA2(ID2(IV))
         CC(NBEG(INC)) = CC(NBEG(INA))
         NTYP(INC) = NDA
         NEND(INC) = NBEG(INC)
         RETURN
      ELSEIF(NTYP(INA).EQ.NCM) THEN
         NC(NBEG(INC)  ) = IA1(ID1(IV))+IA2(ID2(IV))
         NC(NBEG(INC)+1) = 0
         CC(NBEG(INC)  ) = CC(NBEG(INA)  )
         CC(NBEG(INC)+1) = CC(NBEG(INA)+1)
         NTYP(INC) = NCD
         NEND(INC) = NBEG(INC)+1
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAINT')
         RETURN
      ELSEIF(NTYP(INA).EQ.NDA) THEN
         ITA = 1
      ELSEIF(NTYP(INA).EQ.NCD) THEN
         ITA = 2
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
C     IF(LINT.NE.-1) THEN
         IC = NBEG(INC)-ITA
*
         DO 100 IA=NBEG(INA),NEND(INA),ITA
         NCIA = NC(IA)
         IF(IEO(NCIA).GE.NOCUT) GOTO 100
         IFAC = IJJ(IV,NCIA)/IW+1
         IC = IC+ITA
         CC(IC) = CC(IA)/IFAC
         NC(IC) = IA1(IE1(NCIA)+ID1(IV))+IA2(IE2(NCIA)+ID2(IV))
         IF(ITA.EQ.2) THEN
            CC(IC+1) = CC(IA+1)/IFAC
            NC(IC+1) = 0
         ENDIF
 100     CONTINUE
C     ELSE
C        IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
*
C        IF(IW.NE.1) THEN
C           PRINT*,'$$$ ERROR IN TMINT, '//
C    *          'INTEGRAL WITH WEIGHTED DA IS NOT SUPPORTED'
C           CALL FOXDEB
C        ENDIF
*
C        TMT = 0.D0
C        COBN1 = 0.D0
C        COBN2 = 0.D0
C        IC = NBEG(INC)-1
*
C        DO 200 IA=NBEG(INA),NEND(INA)
C        NCIA = NC(IA)
C        IAO = IEO(NCIA)
C        IF(IAO.GT.NOCUT) GOTO 200
C        IFAC = IJJ(IV,NCIA)/IW+1
C        IF(IAO.EQ.NOCUT) THEN
C           CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
C           JJ(IV) = JJ(IV)+1
C           CALL INOBI(JJ,CB1,CB2)
C           CALL INMREO(CB1,CB2,CC(IA),CC1,CC2)
C           CALL INDREO(CC1,CC2,DBFLOAT(IFAC),CB1,CB2,IAF)
C           CALL INAINO(COBN1,COBN2,CB1,CB2,COBN1,COBN2)
C        ELSE
C           IC = IC+1
C           CC(IC) = CC(IA)/IFAC
C           TMT = TMT+ABS(CC(IC))
C           NC(IC) = IA1(IE1(NCIA)+ID1(IV))+IA2(IE2(NCIA)+ID2(IV))
C        ENDIF
C200     CONTINUE
C     ENDIF
*
      NTYP(INC) = NTYP(INA)
      NEND(INC) = IC+ITA-1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAINT')
*
      RETURN
      END
*
      SUBROUTINE DADIDA(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE COMPUTES THE DERIVATIVE OR THE INTEGRAL OF INA
*     WITH RESPECT TO VARIABLE I (INB) AND STORES THE RESULT IN INC.
*     IF INA IS DA OR RE, INC IS DA.  IF INA IS CD OR CM, INC IS CD.
*     IF THE NUMBER GIVEN BY INB IS POSITIVE, IT COMPUTES THE DERIVATIVE.
*     IF                            NEGATIVE,             THE INTEGRAL.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      IV = NINT(CC(NBEG(INB)))
*
      IF(IV.EQ.0) THEN
         PRINT*,'$$$ ERROR IN DADIDA, VARIABLE NUMBER IS 0'
         CALL FOXDEB
      ELSEIF(IV.GT.0) THEN
         CALL DADER(INB,INA,INC)
      ELSE
         CC(NBEG(INB)) = -CC(NBEG(INB))
         CALL DAINT(INB,INA,INC)
         CC(NBEG(INB)) = -CC(NBEG(INB))
      ENDIF
*
      RETURN
      END
*%%
      SUBROUTINE DASGN(INA,INSN,IISN,INC)
*     ***********************************
*
*     THIS SUBROUTINE FLIPS SIGNS OF COEFFICIENTS OF INA TO MAKE
*     THE FIRST NSN LINEAR COEFFS POSITIVE. THE ARRAY ISN WITH LENGTH NSN
*     CONTAINS THE SIGN OF EACH INDEPENDENT VARIABLE FOR FLIPPING.
*     THE RESULT IS INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      COMMON /INBCOM/ IJB(LEA)
      DOUBLE PRECISION CL(LNV)
      INTEGER ISN(LNV),JJ(LNV)
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INSN).NE.NRE) CALL FOXNTY(INSN)
      NSN = NINT(CC(NBEG(INSN)))
*
      IF(NSN.GT.NVMAX) THEN
         PRINT*,'$$$ ERROR IN DASGN, '//
     *          '2ND ARGUMENT > DA DIMENSION NVMAX',NVMAX
         CALL FOXDEB
      ENDIF
*
      CALL DACLIN(INA,CL)
*
      ISNT = 1
      DO 10 I=1,NSN
      IF(CL(I).LT.0) THEN
         ISN(I) = -1
         ISNT = -1
      ELSE
         ISN(I) = 1
      ENDIF
      NTYP(IISN+I) = NRE
      CC(NBEG(IISN+I)) = DBFLOAT(ISN(I))
 10   CONTINUE
*
      IF(ISNT.EQ.1) THEN
         CALL DACOP(INA,INC)
      ELSE
         IC = NBEG(INC)-1
         DO 30 IA=NBEG(INA),NEND(INA)
         IC = IC+1
         NN = NC(IA)
         NC(IC) = NN
         CC(IC) = CC(IA)
         IF(IJB(NN).GE.NSN) GOTO 30
*
         CALL DAENC(IE1(NN),IE2(NN),JJ)
         ICT = 1
         DO 20 I=1,NSN
         IF(ISN(I).EQ.-1.AND.MOD(JJ(I),2).EQ.1) ICT = -ICT
 20      CONTINUE
         IF(ICT.EQ.-1) CC(IC) = -CC(IC)
 30      CONTINUE
*
         NTYP(INC) = NDA
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DASGN')
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DACLIN(INA,CL)
*     *************************
*
*     THIS SUBROUTINE EXTRACTS THE LINEAR COEFFICIENTS CL OF INA (DA OR TM).
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----NON-ABORTING ARITHMETIC -----------------------------------------------
      INTEGER LARI
      COMMON /ARIST/ LARI
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CL(LNV)
      INTEGER JJ(LNV)
*
C     IF(NTYP(INA).NE.NDA.AND.NTYP(INA).NE.NTM) CALL FOXNTY(INA)
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
*
      DO 10 I=1,NVMAX
      CL(I) = 0.D0
 10   CONTINUE
*
      NOCOLD = NOCUT
      NOCUT  = 1
      CALL FOXALL(INA0,1,NMMAX)
*
      IF(NTYP(INA).EQ.NDA) THEN
         CALL DACOP(INA,INA0)
C     ELSEIF(NTYP(INA).EQ.NTM) THEN
C        NVTM = NC(NEND(INA)-1)
C        NETM = NEND(INA)
C        NEDA = NBEG(INA)+NC(NETM)-1
C        IF(NVTM.LT.0) THEN
C           IF(LARI.EQ.1) THEN
C              PRINT*,'$$$ WARNING IN DACLIN, TM ARITHMETIC FAILURE'
C           ELSEIF(LARI.EQ.0) THEN
C              PRINT*,'$$$ ERROR IN DACLIN, TM ARITHMETIC FAILURE'
C              CALL FOXDEB
C           ENDIF
C        ENDIF
*
C        NTYP(INA) = NDA
C        NEND(INA) = NEDA
C        CALL DACOP(INA,INA0)
C        NTYP(INA) = NTM
C        NEND(INA) = NETM
      ENDIF
*
      DO 100 II=NBEG(INA0),NEND(INA0)
      NCIA = NC(II)
      IF(IEO(NCIA).EQ.1) THEN
         CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
         DO 20 I=1,NVMAX
         IF(JJ(I).EQ.1) THEN
            CL(I) = CC(II)
            GOTO 30
         ENDIF
 20      CONTINUE
 30      CONTINUE
      ENDIF
 100  CONTINUE
*
      CALL FOXDAL(INA0,1)
      NOCUT = NOCOLD
*
      RETURN
      END
*
      SUBROUTINE DACLIW(INA,INCL,IMCL)
*     ********************************
*
*     THIS SUBROUTINE EXTRACTS INCL OF THE WEIGHTED DA "LINEAR" COEFFICIENTS OF
*     INA (DA OR TM) AND STORES THEM IN AN ARRAY IMCL.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----NON-ABORTING ARITHMETIC -----------------------------------------------
      INTEGER LARI
      COMMON /ARIST/ LARI
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CL(LNV)
      INTEGER JJ(LNV)
*
      IF(NTYP(INCL).NE.NRE) CALL FOXNTY(INCL)
      NCL = NINT(CC(NBEG(INCL)))
*
      IF(NCL.GT.NVMAX) THEN
         PRINT*,'$$$ ERROR IN DACLIW, '//
     *          '2ND ARGUMENT > DA DIMENSION NVMAX',NVMAX
         CALL FOXDEB
      ELSEIF(NTYP(IMCL).GT.0) THEN
         PRINT*,'$$$ ERROR IN DACLIW, 3RD ARGUMENT MUST BE AN ARRAY.'
         CALL FOXDEB
      ENDIF
*
      IF(LEW.EQ.0) THEN
         CALL DACLIN(INA,CL)
      ELSE
C        IF(NTYP(INA).NE.NDA.AND.NTYP(INA).NE.NTM) CALL FOXNTY(INA)
         IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
*
         IEWM = 0
         DO 10 I=1,NVMAX
         JJ(I) = 0
         IEWM = MAX(IEWM,IEW(I))
 10      CONTINUE
*
         NOCOLD = NOCUT
         NOCUT  = MIN(IEWM,NOCOLD)
         CALL FOXALL(INA0,1,NMMAX)
*
         IF(NTYP(INA).EQ.NDA) THEN
            CALL DACOP(INA,INA0)
C        ELSEIF(NTYP(INA).EQ.NTM) THEN
C           NVTM = NC(NEND(INA)-1)
C           NETM = NEND(INA)
C           NEDA = NBEG(INA)+NC(NETM)-1
C           IF(NVTM.LT.0) THEN
C              IF(LARI.EQ.1) THEN
C                 PRINT*,'$$$ WARNING IN DACLIW, TM ARITHMETIC FAILURE'
C              ELSEIF(LARI.EQ.0) THEN
C                 PRINT*,'$$$ ERROR IN DACLIW, TM ARITHMETIC FAILURE'
C                 CALL FOXDEB
C              ENDIF
C           ENDIF
*
C           NTYP(INA) = NDA
C           NEND(INA) = NEDA
C           CALL DACOP(INA,INA0)
C           NTYP(INA) = NTM
C           NEND(INA) = NETM
         ENDIF
*
         DO 20 I=1,NCL
         CJJ = 0.D0
         IF(IEW(I).LE.NOCUT) THEN
            JJ(I) = IEW(I)
            CALL DADEC(JJ,IC1,IC2)
            IF(IC1.GE.0.AND.IC2.GE.0) THEN
               NCA = IA1(IC1)+IA2(IC2)
               CALL DAPEK(INA0,NCA,CJJ)
            ENDIF
            JJ(I) = 0
         ENDIF
         CL(I) = CJJ
 20      CONTINUE
*
         CALL FOXDAL(INA0,1)
         NOCUT = NOCOLD
      ENDIF
*
      DO 30 I=1,NCL
      NTYP(IMCL+I) = NRE
      CC(NBEG(IMCL+I)) = CL(I)
 30   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE DACQDR(INA,CH,CL,C0)
*     *******************************
*
*     THIS SUBROUTINE EXTRACTS THE COEFFICIENTS UP TO ORDER 2 OF INA (DA OR TM)
*     IN THE FORM OF X*CH*X*1/2+CL*X+C0, WHERE CH IS THE HESSIAN.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----NON-ABORTING ARITHMETIC -----------------------------------------------
      INTEGER LARI
      COMMON /ARIST/ LARI
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CH(LNV,LNV),CL(LNV)
      INTEGER JJ(LNV)
*
C     IF(NTYP(INA).NE.NDA.AND.NTYP(INA).NE.NTM) CALL FOXNTY(INA)
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
*
      DO 10 J=1,NVMAX
      DO 10 I=1,NVMAX
      CH(I,J) = 0.D0
 10   CONTINUE
*
      DO 20 I=1,NVMAX
      CL(I) = 0.D0
 20   CONTINUE
*
      C0 = 0.D0
*
      NOCOLD = NOCUT
      NOCUT  = MIN(NOCUT,2)
      CALL FOXALL(INA0,1,NMMAX)
*
      IF(NTYP(INA).EQ.NDA) THEN
         CALL DACOP(INA,INA0)
C     ELSEIF(NTYP(INA).EQ.NTM) THEN
C        NVTM = NC(NEND(INA)-1)
C        NETM = NEND(INA)
C        NEDA = NBEG(INA)+NC(NETM)-1
C        IF(NVTM.LT.0) THEN
C           IF(LARI.EQ.1) THEN
C              PRINT*,'$$$ WARNING IN DACQDR, TM ARITHMETIC FAILURE'
C           ELSEIF(LARI.EQ.0) THEN
C              PRINT*,'$$$ ERROR IN DACQDR, TM ARITHMETIC FAILURE'
C              CALL FOXDEB
C           ENDIF
C        ENDIF
*
C        NTYP(INA) = NDA
C        NEND(INA) = NEDA
C        CALL DACOP(INA,INA0)
C        NTYP(INA) = NTM
C        NEND(INA) = NETM
      ENDIF
*
      DO 100 II=NBEG(INA0),NEND(INA0)
      NCIA = NC(II)
      NOIA = IEO(NCIA)
      IF(NOIA.EQ.2) THEN
         CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
         I1 = 0
         DO 30 I=1,NVMAX
         IF(JJ(I).EQ.1) THEN
            IF(I1.EQ.0) THEN
               I1 = I
            ELSE
               CH(I1,I) = CC(II)
               CH(I,I1) = CC(II)
               GOTO 40
            ENDIF
         ELSEIF(JJ(I).EQ.2) THEN
            CH(I,I) = 2*CC(II)
            GOTO 40
         ENDIF
 30      CONTINUE
 40      CONTINUE
      ELSEIF(NOIA.EQ.1) THEN
         CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
         DO 50 I=1,NVMAX
         IF(JJ(I).EQ.1) THEN
            CL(I) = CC(II)
            GOTO 60
         ENDIF
 50      CONTINUE
 60      CONTINUE
      ELSEIF(NOIA.EQ.0) THEN
         C0 = CC(II)
      ENDIF
 100  CONTINUE
*
      CALL FOXDAL(INA0,1)
      NOCUT = NOCOLD
*
      RETURN
      END
*
      SUBROUTINE DACQLC(INA,INNA,ICH,ICL,IC0)
*     ***************************************
*
*     THIS SUBROUTINE EXTRACTS THE COEFFICIENTS UP TO WEIGHTED DA ORDER 2 OF
*     INA (DA OR TM) IN THE FORM OF X*CH*X*1/2+CL*X+C0, WHERE CH IS THE HESSIAN.
*     THE RESULT IS STORED IN A 2D ARRAY ICH, A 1D ARRAY ICL AND A REAL IC0.
*     THE ARRAY SIZE NNA MUST BE GIVEN VIA INNA ( NNA >= NVMAX ).
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----NON-ABORTING ARITHMETIC -----------------------------------------------
      INTEGER LARI
      COMMON /ARIST/ LARI
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CH(LNV,LNV),CL(LNV)
      INTEGER JJ(LNV)
*
      IF(NTYP(INNA).NE.NRE) CALL FOXNTY(INNA)
*
      NNA = NINT(CC(NBEG(INNA)))
      IF(NNA.LT.NVMAX) THEN
         PRINT*,'$$$ ERROR IN DACQLC,'//
     *      ' THE ARRAY SIZE (2ND ARGUMENT) IS LESS THAN NVMAX.'
         PRINT*,'    NVMAX = ',NVMAX
         CALL FOXDEB
      ELSEIF(NTYP(ICH).GT.0.OR.NTYP(ICL).GT.0) THEN
         PRINT*,
     *   '$$$ ERROR IN DACQLC, 3RD AND 4TH ARGUMENTS MUST BE ARRAYS.'
         CALL FOXDEB
      ELSEIF(ICH+1.EQ.ICL+1) THEN
         PRINT*,'$$$ WARNING IN DACQLC, '//
     *          'VARIABLES FOR THE 3RD AND 4TH ARGUMENTS MAY BE SAME.'
      ENDIF
*
      DO 20 I=1,NNA
      JCH = ICH+I
      DO 10 J=1,NNA
      NTYP(JCH) = NRE
      CC(NBEG(JCH)) = 0.D0
      JCH = JCH+NNA
 10   CONTINUE
      NTYP(ICL+I) = NRE
      CC(NBEG(ICL+I)) = 0.D0
 20   CONTINUE
*
      IF(LEW.EQ.0) THEN
         CALL DACQDR(INA,CH,CL,C0)
      ELSE
C        IF(NTYP(INA).NE.NDA.AND.NTYP(INA).NE.NTM) CALL FOXNTY(INA)
         IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
*
         IEWM = 0
         DO 30 I=1,NVMAX
         JJ(I) = 0
         IEWM = MAX(IEWM,IEW(I))
 30      CONTINUE
*
         NOCOLD = NOCUT
         NOCUT  = MIN(2*IEWM,NOCOLD)
         CALL FOXALL(INA0,1,NMMAX)
*
         IF(NTYP(INA).EQ.NDA) THEN
            CALL DACOP(INA,INA0)
C        ELSEIF(NTYP(INA).EQ.NTM) THEN
C           NVTM = NC(NEND(INA)-1)
C           NETM = NEND(INA)
C           NEDA = NBEG(INA)+NC(NETM)-1
C           IF(NVTM.LT.0) THEN
C              IF(LARI.EQ.1) THEN
C                 PRINT*,'$$$ WARNING IN DACQLC, TM ARITHMETIC FAILURE'
C              ELSEIF(LARI.EQ.0) THEN
C                 PRINT*,'$$$ ERROR IN DACQLC, TM ARITHMETIC FAILURE'
C                 CALL FOXDEB
C              ENDIF
C           ENDIF
*
C           NTYP(INA) = NDA
C           NEND(INA) = NEDA
C           CALL DACOP(INA,INA0)
C           NTYP(INA) = NTM
C           NEND(INA) = NETM
         ENDIF
*
         C0 = RE(INA0)
*
         DO 40 I=1,NVMAX
*
         CJJ = 0.D0
         IF(IEW(I).LE.NOCUT) THEN
            JJ(I) = IEW(I)
            CALL DADEC(JJ,IC1,IC2)
            IF(IC1.GE.0.AND.IC2.GE.0) THEN
               NCA = IA1(IC1)+IA2(IC2)
               CALL DAPEK(INA0,NCA,CJJ)
            ENDIF
            JJ(I) = 0
         ENDIF
         CL(I) = CJJ
*
         CJJ = 0.D0
         IF(2*IEW(I).LE.NOCUT) THEN
            JJ(I) = 2*IEW(I)
            CALL DADEC(JJ,IC1,IC2)
            IF(IC1.GE.0.AND.IC2.GE.0) THEN
               NCA = IA1(IC1)+IA2(IC2)
               CALL DAPEK(INA0,NCA,CJJ)
            ENDIF
            JJ(I) = 0
         ENDIF
         CH(I,I) = 2*CJJ
*
         DO 40 J=I+1,NVMAX
         CJJ = 0.D0
         IF(IEW(I)+IEW(J).LE.NOCUT) THEN
            JJ(I) = IEW(I)
            JJ(J) = IEW(J)
            CALL DADEC(JJ,IC1,IC2)
            IF(IC1.GE.0.AND.IC2.GE.0) THEN
               NCA = IA1(IC1)+IA2(IC2)
               CALL DAPEK(INA0,NCA,CJJ)
            ENDIF
            JJ(I) = 0
            JJ(J) = 0
         ENDIF
         CH(I,J) = CJJ
         CH(J,I) = CJJ
*
 40      CONTINUE
*
         CALL FOXDAL(INA0,1)
         NOCUT = NOCOLD
      ENDIF
*
      DO 60 I=1,NVMAX
      JCH = ICH+I
      DO 50 J=1,NVMAX
      CC(NBEG(JCH)) = CH(I,J)
      JCH = JCH+NNA
 50   CONTINUE
      CC(NBEG(ICL+I)) = CL(I)
 60   CONTINUE
*
      NTYP(IC0) = NRE
      CC(NBEG(IC0)) = C0
*
      RETURN
      END
*
*%%
      SUBROUTINE DAPOI2(INA,INB,INC,N)
*     ********************************
*
*     THIS SUBROUTINE COMPUTES THE POISSON BRACKET OF THE VECTORS A AND
*     B AND STORES THE RESULT IN C. N IS THE DEGREE OF FREEDOM OF THE SYSTEM.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      INTEGER IS(3)
*
      IF((INA.EQ.INC).OR.(INB.EQ.INC)) CALL DANFI(INA,INB,INC)
*
      CALL FOXALL(IS,3,NMMAX)
      CALL DACON(INC,0.D0)
*
      DO 100 I=1,N
*
      CALL DADER(2*I-1,INA,IS(1))
      CALL DADER(2*I,  INB,IS(2))
      CALL DAMDA(IS(1),IS(2),IS(3))
      CALL DAADA(INC,IS(3),IS(1))
      CALL DACOP(IS(1),INC)
*
      CALL DADER(2*I,  INA,IS(1))
      CALL DADER(2*I-1,INB,IS(2))
      CALL DAMDA(IS(1),IS(2),IS(3))
      CALL DASDA(INC,IS(3),IS(1))
      CALL DACOP(IS(1),INC)
*
 100  CONTINUE
*
      CALL FOXDAL(IS,3)
*
      RETURN
      END
*%%
      SUBROUTINE DAPLU(INA,IIV,INB,INC)
*     *********************************
*
*     THIS SUBROUTINE REPLACES VARIABLE IV IN INA WITH THE CONSTANT IN INB.
*     IF WEIGHTED DA IS USED FOR THE VARIABLE IV, THE WEIGHTED VARIABLE IV IS
*     REPLACED WITH THE CONSTANT IN INB.
*     (PARTIAL PLUG)   THE OPERATION IS COMPATIBLE WITH POLVAL.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION JJ(LNV),ID1(0:LNO),ID2(0:LNO),CR(0:LNO),CI(0:LNO)
*
      IF(NTYP(IIV).NE.NRE) CALL FOXNTY(IIV)
      IV = NINT(CC(NBEG(IIV)))
      IF(IV.LT.1.OR.IV.GT.NVMAX) THEN
         PRINT*,'$$$ ERROR IN DAPLU, VARIABLE IS OUT OF RANGE'
         PRINT*,'    VARIABLE INDEX, BOUND = ',IV,',',NVMAX
         CALL FOXDEB
      ENDIF
*
      IF(NTYP(INA).EQ.NDA) THEN
         IF(IEW(IV).GT.NOMAX) THEN
            CALL DACOP(INA,INC)
            RETURN
         ENDIF
         ITA = 1
      ELSEIF(NTYP(INA).EQ.NCD) THEN
         IF(IEW(IV).GT.NOMAX) THEN
            CALL CDCOP(INA,INC)
            RETURN
         ENDIF
         ITA = 2
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      IF(NTYP(INB).EQ.NRE) THEN
         ITB = 1
      ELSEIF(NTYP(INB).EQ.NCM) THEN
         ITB = 2
      ELSE
         CALL FOXNTY(INB)
      ENDIF
*
      CALL CDCLR
      ITC    = MAX(ITA,ITB)
      CR (0) = 1.D0
      CI (0) = 0.D0
      ID1(0) = 0
      ID2(0) = 0
      CR (1) = CC(NBEG(INB))
      CI (1) = 0.D0
      IF(ITB.EQ.2) CI(1) = CC(NBEG(INB)+1)
*
      DO 10 J=1,NVMAX
  10  JJ(J) = 0
*
      DO 20 J=1,NOMAX/IEW(IV)
      CR(J) = CR(J-1)*CR(1)-CI(J-1)*CI(1)
      CI(J) = CI(J-1)*CR(1)+CR(J-1)*CI(1)
      JJ(IV) = J*IEW(IV)
      CALL DADEC(JJ,II1,II2)
      ID1(J) = II1
      ID2(J) = II2
  20  CONTINUE
*
      DO 100 IA=NBEG(INA),NEND(INA),ITA
*
      NCIA = NC(IA)
      IE   = IJJ(IV,NCIA)/IEW(IV)
      CFAR = CR(IE)
      CFAI = CI(IE)
      NNC = ITC*(IA1(IE1(NCIA)-ID1(IE))+IA2(IE2(NCIA)-ID2(IE)))+1-ITC
      CDA(NNC  ) = CDA(NNC  )+CC(IA)*CFAR
      IF(ITC.EQ.1) GOTO 100
      CCIAI      = 0.D0
      IF(ITA.EQ.2) CCIAI = CC(IA+1)
      CDA(NNC  ) = CDA(NNC  )              -CCIAI*CFAI
      CDA(NNC+1) = CDA(NNC+1)+CC(IA)*CFAI+CCIAI*CFAR
*
 100  CONTINUE
*
      IF(ITC.EQ.1) THEN
         CALL DAPAC(INC)
      ELSE
         CALL CDPAC(INC)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DASCL(INA,IIV,INB,INC)
*     *********************************
*
*     THIS SUBROUTINE SCALES THE VARIABLES X_IV IN THE DA INA BY B*X_IV.
*     THE RESULT IS A DA IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION BR(0:LNO)
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(IIV).NE.NRE) CALL FOXNTY(IIV)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IV = NINT(CC(NBEG(IIV)))
      IF(IV.LT.1.OR.IV.GT.NVMAX) THEN
         PRINT*,'$$$ ERROR IN DASCL, VARIABLE IS OUT OF RANGE'
         PRINT*,'    VARIABLE INDEX, BOUND = ',IV,',',NVMAX
         CALL FOXDEB
      ELSEIF(IEW(IV).GT.NOMAX) THEN
         CALL DACOP(INA,INC)
         RETURN
      ENDIF
*
      B = CC(NBEG(INB))
      BR(0)  = 1.D0
      DO 10 I=1,NOMAX/IEW(IV)
      BR(I)  = BR(I-1)*B
 10   CONTINUE
*
      IC = NBEG(INC)-1
      DO 20 IA=NBEG(INA),NEND(INA)
      NCIA = NC(IA)
      IE   = IJJ(IV,NCIA)/IEW(IV)
      CCC = CC(IA)*BR(IE)
      IF(ABS(CCC).LT.EPS) GOTO 20
      IC = IC+1
      CC(IC) = CCC
      NC(IC) = NCIA
 20   CONTINUE
*
      NTYP(INC) = NDA
      NEND(INC) = IC
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DASCL')
*
      RETURN
      END
*
      SUBROUTINE DABINO
*     *****************
*
*     THIS SUBROUTINE COMPUTES THE BINOMIAL COEFFICIENTS AND INITIAL DATA
*     FOR DATRN AND TMTRNF
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----DATA FOR DATRN, TMTRNF ------------------------------------------------
      PARAMETER(LBC=2000)
      INTEGER IBINO(0:LNO,0:LNO),ICBIN(LBC),ITRD1(0:LNO),ITRD2(0:LNO),
     *        ITREC(0:LNO),ITREA(0:LNO),JJTR(LNV),ICALLT
      DOUBLE PRECISION TRAR(0:LNO),TRCR(0:LNO),TRINUP,CBIN(LBC)
      COMMON /TRNCOM/ TRAR,TRCR,TRINUP,CBIN,IBINO,ICBIN,
     *       ITRD1,ITRD2,ITREC,ITREA,JJTR,ICALLT
*----------------------------------------------------------------------------
*-----ROUNDING --------------------------------------------------------------
      INTEGER INTINI
      DOUBLE PRECISION RINUP(2),RINDN(2),FINSR
      COMMON /INCOM/ RINUP,RINDN,FINSR,INTINI
*----------------------------------------------------------------------------
*
      ICALLT = 1
*
      DO 10 N=0,LNO
      IBINO(N,0) = 1
      IBINO(N,N) = 1
 10   CONTINUE
*
      ICB = 0
      DO 40 N=2,LNO
      IBINO(N,1) = N
      CFAC = DBFLOAT(N)
      ICT = 0
      DO 20 I=2,N/2
      CFAC = (CFAC*DBFLOAT(N-I+1))/DBFLOAT(I)
      IFAC = NINT(CFAC)
      IF(CFAC.LE.2147483647.D0.AND.ABS(CFAC-IFAC).EQ.0.D0) THEN
         IBINO(N,I) = IFAC
      ELSE
         ICB = ICB+1
         ICT = ICT+1
         IF(ICB.GT.LBC) THEN
            PRINT*,'!!! MEMORY EXHAUSTION IN DABINO, INCREASE LBC'
            CALL FOXSTL
         ENDIF
         IBINO(N,I) = -ICB
         CBIN(ICB) = CFAC
         ICBIN(ICB) = 2*ICT
      ENDIF
 20   CONTINUE
      DO 30 I=N/2+1,N-1
      IBINO(N,I) = IBINO(N,N-I)
 30   CONTINUE
 40   CONTINUE
*
*     DIAGNOSTICS
*
C     WRITE(6,'(A,I4/A)') ' IN DABINO: LNO = ',LNO,
C    *   '   N   I  IBINO(N,I)     CBIN(-IBINO(N,I)) ICBIN(-IBINO(N,I))'
      DO 50 N=0,LNO
C     WRITE(6,'()')
      DO 50 I=0,N
      IF(IBINO(N,I).EQ.0) THEN
         PRINT*,'@@@ ERROR IN DABINO, IBINO DATA WRONG.'
         CALL FOXDEB
C     ELSEIF(IBINO(N,I).LT.0) THEN
C        WRITE(6,'(2I4,I12,G27.16E3,I8)')
C    *      N,I,IBINO(N,I),CBIN(-IBINO(N,I)),ICBIN(-IBINO(N,I))
C     ELSE
C        WRITE(6,'(2I4,I12)') N,I,IBINO(N,I)
      ENDIF
 50   CONTINUE
*
      TRAR(0)  = 1.D0
      TRCR(0)  = 1.D0
      ITRD1(0) = 0
      ITRD2(0) = 0
      ITREC(0) = 0
      ITREA(0) = 1
      ITREC(1) = 1
      ITREA(1) = 1
      DO 60 I=2,LNO
      ITREC(I) = I+1
      ITREA(I) = I+1
 60   CONTINUE
*
      DO 70 I=1,LNV
      JJTR(I) = 0
 70   CONTINUE
*
C     TRINUP = 1.D0+(RINUP(1)-1.D0)*50
*
      RETURN
      END
*
      SUBROUTINE DATRN(INA,IMA,IMC,INM1,INM2,INC)
*     *******************************************
*
*     THIS SUBROUTINE REPLACES THE VARIABLES X(I) IN THE DA INA WITH
*     A(I)*X(I)+C(I) FOR I=M1,...,M2. THE RESULT IS A DA IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----DATA FOR DATRN, TMTRNF ------------------------------------------------
      PARAMETER(LBC=2000)
      INTEGER IBINO(0:LNO,0:LNO),ICBIN(LBC),ITRD1(0:LNO),ITRD2(0:LNO),
     *        ITREC(0:LNO),ITREA(0:LNO),JJTR(LNV),ICALLT
      DOUBLE PRECISION TRAR(0:LNO),TRCR(0:LNO),TRINUP,CBIN(LBC)
      COMMON /TRNCOM/ TRAR,TRCR,TRINUP,CBIN,IBINO,ICBIN,
     *       ITRD1,ITRD2,ITREC,ITREA,JJTR,ICALLT
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INM1).NE.NRE) CALL FOXNTY(INM1)
      IF(NTYP(INM2).NE.NRE) CALL FOXNTY(INM2)
      MM1 = NINT(CC(NBEG(INM1)))
      MM2 = NINT(CC(NBEG(INM2)))
      M1 = MAX(1,MIN(MM1,MM2))
      M2 = MIN(MAX(MM1,MM2),NVMAX)
*
      IF(NTYP(IMA).GT.0.OR.NTYP(IMC).GT.0) THEN
         PRINT*,'$$$ ERROR IN DATRN, 2ND OR 3RD ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ENDIF
*
      IF(ICALLT.EQ.0) CALL DABINO
*
      CALL FOXALL(ISC,1,NMMAX)
      IS = NBEG(ISC)
      DO 20 IA=NBEG(INA),NEND(INA)
      CC(IS) = CC(IA)
      NC(IS) = NC(IA)
      IS = IS+1
 20   CONTINUE
      NEND(ISC) = IS-1
      NTYP(ISC) = NDA
*
*     PROCESSING EACH INDEPENDENT VARIABLE
*     ************************************
*
      DO 60 IV=M1,M2
      IF(IEW(IV).GT.NOMAX) GOTO 60
*
*     REPLACING THE INDEPENDENT VARIABLE X(IV) BY A(IV)*X(IV)+C(IV)
*     *************************************************************
*
      IF(NTYP(IMA+IV).NE.NRE) CALL FOXNTY(IMA+IV)
      IF(NTYP(IMC+IV).NE.NRE) CALL FOXNTY(IMC+IV)
      AKON = CC(NBEG(IMA+IV))
      CKON = CC(NBEG(IMC+IV))
      IF(AKON.EQ.1.D0.AND.CKON.EQ.0.D0) GOTO 60
*
      CALL DACLR
*
      DO 30 I=1,NOMAX/IEW(IV)
      TRAR(I)  = TRAR(I-1)*AKON
      TRCR(I)  = TRCR(I-1)*CKON
      JJTR(IV) = I*IEW(IV)
      CALL DADEC(JJTR,II1,II2)
      ITRD1(I) = II1
      ITRD2(I) = II2
 30   CONTINUE
*
      DO 50 IS=NBEG(ISC),NEND(ISC)
      NCIS = NC(IS)
      CCIS = CC(IS)
      IE   = IJJ(IV,NCIS)/IEW(IV)
*
      DO 40 I=0,IE
      IF(IE.EQ.0) THEN
         CK = CCIS
      ELSEIF(IE.EQ.1) THEN
         IF(I.EQ.0) THEN
            CK = CCIS*TRAR(1)
         ELSE
            CK = CCIS*TRCR(1)
         ENDIF
      ELSE
         IF(IBINO(IE,I).LT.0) THEN
            CK = CBIN(-IBINO(IE,I))*CCIS*TRAR(IE-I)*TRCR(I)
         ELSE
            CK = IBINO(IE,I)*CCIS*TRAR(IE-I)*TRCR(I)
         ENDIF
      ENDIF
      IF(CK.EQ.0.D0) GOTO 40
      NNC = IA1(IE1(NCIS)-ITRD1(I))+IA2(IE2(NCIS)-ITRD2(I))
      CDA(NNC) = CDA(NNC)+CK
 40   CONTINUE
 50   CONTINUE
*
      JJTR(IV) = 0
      CALL DAPAC(ISC)
*
 60   CONTINUE
*
      IC = NBEG(INC)
      DO 70 IS=NBEG(ISC),NEND(ISC)
      CC(IC) = CC(IS)
      NC(IC) = NC(IS)
      IC = IC+1
 70   CONTINUE
*
      NTYP(INC) = NDA
      NEND(INC) = IC-1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DATRN')
*
      CALL FOXDAL(ISC,1)
*
      RETURN
      END
*
      SUBROUTINE DAPOI(INA,INB,INC,IND)
*     *********************************
*
*     THIS SUBROUTINE PERFORMS A POISSON BRACKET THE DA VECTORS A AND B.
*     THE RESULT IS STORED IN INC. 2*ND IS THE DIMENSIONALITY OF PHASE SPACE.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION IBNO(0:LNO,LNV),IENO(0:LNO,LNV)
      COMMON /MONDA/ ID1(LNV),ID2(LNV)
      INTEGER JJ(LNV)
*
      IF(LEW.EQ.1) THEN
         PRINT*,'$$$ ERROR IN DAPOI, '//
     *          'WEIGHTED DA (DANOTW) IS NOT SUPPORTED'
         CALL FOXDEB
      ENDIF
*
      IF((NTYP(INA).EQ.NRE).OR.(NTYP(INB).EQ.NRE)) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
         RETURN
      ENDIF
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NDA) CALL FOXNTY(INB)
      IF(NTYP(IND).NE.NRE) CALL FOXNTY(IND)
      ND = NINT(CC(NBEG(IND)))
*
      CALL DACLR
*
      JMEM = IMEM
*
      DO 20 J=1,2*ND
      DO 20 I=0,NOCUT
      IBNO(I,J) = JMEM+1
      IENO(I,J) = JMEM
  20  JMEM = JMEM+NMMAX
*
      CALL MEMCHK(JMEM)
*
*     DIFFERENTIATING B BY ORDER
*     **************************
*
      DO 40 IB=NBEG(INB),NEND(INB)
*
      NCIB = NC(IB)
      NOIB = IEO(NCIB)-1
*
      JFAC = 1
      CALL DAENC(IE1(NCIB),IE2(NCIB),JJ)
      DO 35 J=1,2*ND
      JFAC = -JFAC
      JDFA = JJ(J)
      IF(JDFA.EQ.0) GOTO 35
      JS = J-JFAC
      IPOS = IENO(NOIB,JS)+1
      IENO(NOIB,JS) = IPOS
*
      CC(IPOS) = JFAC*JDFA*CC(IB)
      NC(IPOS) = IA1(IE1(NCIB)-ID1(J))+IA2(IE2(NCIB)-ID2(J))
  35  CONTINUE
  40  CONTINUE
*
*     PERFORMING MULTIPLICATIONS
*     **************************
*
      DO 100 IA=NBEG(INA),NEND(INA)
*
      NCIA = NC(IA)
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
      DO 100 J=1,2*ND
*
      JDFA = JJ(J)
      IF(JDFA.EQ.0) GOTO 100
      CCIA = JDFA*CC(IA)
      I1IA = IE1(NCIA)-ID1(J)
      I2IA = IE2(NCIA)-ID2(J)
*
      DO 90 NOIB = 0,NOCUT-IEO(NCIA)+1
      DO 90 IB = IBNO(NOIB,J),IENO(NOIB,J)
*
      NCIB = NC(IB)
      ICC = IA2(I2IA+IE2(NCIB))+IA1(I1IA+IE1(NCIB))
      CDA(ICC) = CDA(ICC)+CCIA*CC(IB)
*
  90  CONTINUE
 100  CONTINUE
*
      CALL DAPAC(INC)
*
      RETURN
      END
*
      SUBROUTINE DAGMD(INA,INB,INC,IND)
*     *********************************
*
*     THIS SUBROUTINE COMPUTES GRAD(INA)*INB, WHERE INA IS A DA VECTOR AND
*     INB IS AN ARRAY OF DA VECTORS.
*     THE RESULT IS STORED IN INC. ND IS THE LENGTH OF INB.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION IBNO(0:LNO,LNV),IENO(0:LNO,LNV),NB(LNV)
      COMMON /MONDA/ ID1(LNV),ID2(LNV)
      INTEGER JJ(LNV)
*
      IF((NTYP(INA).EQ.NRE)) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
         RETURN
      ENDIF
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(IND).NE.NRE) CALL FOXNTY(IND)
      JND = NINT(CC(NBEG(IND)))
      IF(JND.GT.LNV) THEN
         PRINT*,'$$$ ERROR IN DAGMD, ND > LNV. INCREASE LNV > ',JND-1
         CALL FOXDEB
      ELSEIF(NTYP(INB).GT.0) THEN
         PRINT*,'$$$ ERROR IN DAGMD, SECOND ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ENDIF
*
      DO 15 I=1,JND
      NB(I) = INB+I
  15  IF(NTYP(NB(I)).NE.NDA) CALL FOXNTY(NB(I))
*
      CALL DACLR
      JMEM = IMEM
*
      DO 20 J=1,JND
      DO 20 I=0,NOCUT
      IBNO(I,J) = JMEM+1
      IENO(I,J) = JMEM
  20  JMEM = JMEM+NMMAX
*
      CALL MEMCHK(JMEM)
*
*     DIFFERENTIATING A BY ORDER
*     **************************
*
      DO 40 IA=NBEG(INA),NEND(INA)
*
      NCIA = NC(IA)
      NOIA = IEO(NCIA)
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
*
      DO 35 J=1,JND
      IF(JJ(J).EQ.0) GOTO 35
      JDFA = JJ(J)/IEW(J)
      NOIAJ = NOIA-IEW(J)
      IPOS = IENO(NOIAJ,J)+1
      IENO(NOIAJ,J) = IPOS
*
      CC(IPOS) = JDFA*CC(IA)
      NC(IPOS) = IA1(IE1(NCIA)-ID1(J))+IA2(IE2(NCIA)-ID2(J))
  35  CONTINUE
  40  CONTINUE
*
*     PERFORMING MULTIPLICATIONS
*     **************************
*
      DO 100 J=1,JND
      DO 100 IB=NBEG(NB(J)),NEND(NB(J))
*
      NCIB = NC(IB)
*
      CCIB = CC(IB)
      I1IB = IE1(NCIB)
      I2IB = IE2(NCIB)
*
      DO 90 NOIA = 0,NOCUT-IEO(NCIB)
      DO 90 IA = IBNO(NOIA,J),IENO(NOIA,J)
*
      NCIA = NC(IA)
      ICC = IA2(I2IB+IE2(NCIA))+IA1(I1IB+IE1(NCIA))
      CDA(ICC) = CDA(ICC)+CCIB*CC(IA)
*
  90  CONTINUE
 100  CONTINUE
*
      CALL DAPAC(INC)
*
      RETURN
      END
*
      SUBROUTINE CDGMD(INA,INB,INC,IND)
*     *********************************
*
*     THIS SUBROUTINE COMPUTES GRAD(INA)*INB, WHERE INA IS A CD VECTOR AND
*     INB IS AN ARRAY OF CD VECTORS.
*     THE RESULT IS STORED IN INC. ND IS THE LENGTH OF INB.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION IBNO(0:LNO,LNV),IENO(0:LNO,LNV),NB(LNV)
      COMMON /MONDA/ ID1(LNV),ID2(LNV)
      INTEGER JJ(LNV)
*
      IF((NTYP(INA).EQ.NRE)) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
         RETURN
      ENDIF
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IF(NTYP(IND).NE.NRE) CALL FOXNTY(IND)
      JND = NINT(CC(NBEG(IND)))
      IF(JND.GT.LNV) THEN
         PRINT*,'$$$ ERROR IN CDGMD, ND > LNV. INCREASE LNV > ',JND-1
         CALL FOXDEB
      ELSEIF(NTYP(INB).GT.0) THEN
         PRINT*,'$$$ ERROR IN CDGMD, SECOND ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ENDIF
*
      DO 15 I=1,JND
      NB(I) = INB+I
  15  IF(NTYP(NB(I)).NE.NCD) CALL FOXNTY(NB(I))
*
      CALL CDCLR
      JMEM = IMEM
*
      DO 20 J=1,JND
      DO 20 I=0,NOCUT
      IBNO(I,J) = JMEM+1
      IENO(I,J) = JMEM-1
  20  JMEM = JMEM+2*NMMAX
*
      CALL MEMCHK(JMEM)
*
*     DIFFERENTIATING A BY ORDER
*     **************************
*
      DO 40 IA=NBEG(INA),NEND(INA),2
*
      NCIA = NC(IA)
      NOIA = IEO(NCIA)
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
*
      DO 35 J=1,JND
      IF(JJ(J).EQ.0) GOTO 35
      JDFA = JJ(J)/IEW(J)
      NOIAJ = NOIA-IEW(J)
      IPOS = IENO(NOIAJ,J)+2
      IENO(NOIAJ,J) = IPOS
*
      CC(IPOS  ) = JDFA*CC(IA  )
      CC(IPOS+1) = JDFA*CC(IA+1)
      NC(IPOS  ) = IA1(IE1(NCIA)-ID1(J))+IA2(IE2(NCIA)-ID2(J))
      NC(IPOS+1) = 0
  35  CONTINUE
  40  CONTINUE
*
*     PERFORMING MULTIPLICATIONS
*     **************************
*
      DO 100 J=1,JND
      DO 100 IB=NBEG(NB(J)),NEND(NB(J)),2
*
      NCIB = NC(IB)
      I1IB = IE1(NCIB)
      I2IB = IE2(NCIB)
      CCIR = CC(IB  )
      CCII = CC(IB+1)
*
      DO 90 NOIA = 0,NOCUT-IEO(NCIB)
      DO 90 IA = IBNO(NOIA,J),IENO(NOIA,J),2
*
      NCIA = NC(IA)
      ICC = 2*(IA2(I2IB+IE2(NCIA))+IA1(I1IB+IE1(NCIA)))-1
      CDA(ICC  ) = CDA(ICC  )+CCIR*CC(IA  )-CCII*CC(IA+1)
      CDA(ICC+1) = CDA(ICC+1)+CCIR*CC(IA+1)+CCII*CC(IA  )
*
  90  CONTINUE
 100  CONTINUE
*
      CALL CDPAC(INC)
*
      RETURN
      END
*
      SUBROUTINE CDGMD2(INA,INB,INC,IND)
*     *********************************
*
*     THIS SUBROUTINE COMPUTES GRAD(INA)*INB, WHERE INA IS A CD VECTOR AND
*     INB IS AN ARRAY OF CD VECTORS.
*     THE RESULT IS STORED IN INC. ND IS THE LENGTH OF INB.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION INBR(LNV),INBI(LNV),INCR(LNV),INCI(LNV)
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IF(NTYP(IND).NE.NRE) CALL FOXNTY(IND)
      IF(NTYP(INB).GT.0) THEN
         PRINT*,'$$$ ERROR IN CDGMD2, SECOND ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ENDIF
      JND = NINT(CC(NBEG(IND)))
*
      CALL FOXALL(INAR,1,  NMMAX)
      CALL FOXALL(INAI,1,  NMMAX)
      CALL FOXALL(ISC1,1,  NMMAX)
      CALL FOXALL(ISC2,1,2*NMMAX)
      CALL FOXALL(IIMA,1,2      )
      DO 10 I=1,JND
      CALL FOXALL(INBR(I),1,NMMAX)
      CALL FOXALL(INBI(I),1,NMMAX)
      CALL FOXALL(INCR(I),1,NMMAX)
      CALL FOXALL(INCI(I),1,NMMAX)
  10  CONTINUE
*
      CC(NBEG(IIMA)  ) = 0.D0
      CC(NBEG(IIMA)+1) = 1.D0
      NEND(IIMA) = NBEG(IIMA)+1
      NTYP(IIMA) = NCM
*
      CALL CDRE(INA,INAR)
      CALL CDIM(INA,INAI)
      DO 20 I=1,JND
      CALL CDRE(INB+I,INBR(I))
      CALL CDIM(INB+I,INBI(I))
  20  CONTINUE
*
      CALL DAGMD(INAR,INBR(1)-1,ISC1,IND)
      CALL DAGMD(INAI,INBI(1)-1,ISC2,IND)
      CALL DASDA(ISC1,ISC2,INCR)
      CALL DAGMD(INAR,INBI(1)-1,ISC1,IND)
      CALL DAGMD(INAI,INBR(1)-1,ISC2,IND)
      CALL DAADA(ISC1,ISC2,INCI)
*
      CALL DAMCM(INCI,IIMA,ISC1)
      CALL DAACD(INCR,ISC1,INC)
*
      CALL FOXDAL(IIMA,1  )
      CALL FOXDAL(ISC2,1  )
      CALL FOXDAL(ISC1,1  )
      DO 30 I=1,JND
      CALL FOXDAL(INCI(I),1)
      CALL FOXDAL(INCR(I),1)
      CALL FOXDAL(INBI(I),1)
      CALL FOXDAL(INBR(I),1)
  30  CONTINUE
      CALL FOXDAL(INAI,1  )
      CALL FOXDAL(INAR,1  )
*
      RETURN
      END
*
      SUBROUTINE DALEX(INA,INB,INC,IND)
*     *********************************
*
*     THIS SUBROUTINE PERFORMS A LIE EXPONENTIATION EXP(:A:) ON B. THE RESULT
*     IS STORED IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      EPS = EPS/10.D0
      CALL FOXALL(IPB,1,NMMAX)
      CALL FOXALL(ISC,1,NMMAX)
      CALL FOXALL(ICC,1,1)
      CALL FOXALL(ICF,1,1)
*
      NTYP(ICC) = 1
      NTYP(ICF) = 1
      NBICC = NBEG(ICC)
      NBICF = NBEG(ICF)
*
      CALL DAPOI(INA,INB,IPB,IND)
      CALL DAADA(INB,IPB,INC)
*
      CFAC = 1.D0
*
  10  CFAC = CFAC+1.D0
      CC(NBICF) = 1.D0/CFAC
      CALL DAMRE(IPB,ICF,ISC)
      CALL DAPOI(INA,ISC,IPB,IND)
      CALL DAADA(IPB,INC,ISC)
      CALL DACOP(ISC,INC)
      CALL DANOR(IPB,ICC)
      ILEN = NEND(IPB)-NBEG(IPB)+1
      IF(CFAC.GT.200) THEN
         PRINT*,'$$$ ERROR IN DALEX, NO CONVERGENCE AFTER 200 STEPS'
         CALL FOXDEB
      ELSEIF(CC(NBICC).GT.EPS) THEN
         GOTO 10
      ENDIF
*
*     PRINT*,'REQUIRED STEPS IN DALEX: ',NINT(CFAC)
*
      CALL FOXDAL(ICF,1)
      CALL FOXDAL(ICC,1)
      CALL FOXDAL(ISC,1)
      CALL FOXDAL(IPB,1)
      EPS = EPS*10.D0
*
      RETURN
      END
*
      SUBROUTINE DAFLO(INA,INB,INC,IND)
*     *********************************
*
*     THIS SUBROUTINE COMPUTES THE FLOW OF X' = F(X,T)
*     INB: INITIAL CONDITION, INA: ARRAY OF RIGHT HAND SIDES.
*     THE RESULT IS STORED IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      EPS = EPS/10.D0
      CALL FOXALL(IFL,1,NMMAX)
      CALL FOXALL(ISC,1,NMMAX)
      CALL FOXALL(ICC,1,1)
      CALL FOXALL(ICF,1,1)
*
      NTYP(ICC) = 1
      NTYP(ICF) = 1
      NBICC = NBEG(ICC)
      NBICF = NBEG(ICF)
*
      CALL DAGMD(INB,INA,IFL,IND)
      CALL DAADA(INB,IFL,INC)
*
      CFAC = 1.D0
*
  10  CFAC = CFAC+1.D0
      CC(NBICF) = 1.D0/CFAC
      CALL DAMRE(IFL,ICF,ISC)
      CALL DAGMD(ISC,INA,IFL,IND)
      CALL DAADA(IFL,INC,ISC)
      CALL DACOP(ISC,INC)
      CALL DANOR(IFL,ICC)
      ILEN = NEND(IFL)-NBEG(IFL)+1
      IF(CFAC.GT.200) THEN
         PRINT*,'$$$ ERROR IN DAFLO, NO CONVERGENCE AFTER 200 STEPS'
         CALL FOXDEB
      ELSEIF(CC(NBICC).GT.EPS) THEN
         GOTO 10
      ENDIF
*
      CALL FOXDAL(ICF,1)
      CALL FOXDAL(ICC,1)
      CALL FOXDAL(ISC,1)
      CALL FOXDAL(IFL,1)
      EPS = EPS*10.D0
*
      RETURN
      END
*
      SUBROUTINE CDFLO(INA,INB,INC,IND)
*     *********************************
*
*     THIS SUBROUTINE COMPUTES THE FLOW OF X' = F(X,T)
*     INB: INITIAL CONDITION, INA: ARRAY OF RIGHT HAND SIDES.
*     THE RESULT IS STORED IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      EPS = EPS/10.D0
      CALL FOXALL(IFL,1,2*NMMAX)
      CALL FOXALL(ISC,1,2*NMMAX)
      CALL FOXALL(ICC,1,1)
      CALL FOXALL(ICF,1,1)
*
      NTYP(ICC) = 1
      NTYP(ICF) = 1
      NBICC = NBEG(ICC)
      NBICF = NBEG(ICF)
*
      CALL CDGMD(INB,INA,IFL,IND)
      CALL CDACD(INB,IFL,INC)
*
      CFAC = 1.D0
*
  10  CFAC = CFAC+1.D0
      CC(NBICF) = 1.D0/CFAC
      CALL CDMRE(IFL,ICF,ISC)
      CALL CDGMD(ISC,INA,IFL,IND)
      CALL CDACD(IFL,INC,ISC)
      CALL CDCOP(ISC,INC)
      CALL CDNOR(IFL,ICC)
      ILEN = NEND(IFL)-NBEG(IFL)+1
      IF(CFAC.GT.200) THEN
         PRINT*,'$$$ ERROR IN CDFLO, NO CONVERGENCE AFTER 200 STEPS'
         CALL FOXDEB
      ELSEIF(CC(NBICC).GT.EPS) THEN
         GOTO 10
      ENDIF
*
      CALL FOXDAL(ICF,1)
      CALL FOXDAL(ICC,1)
      CALL FOXDAL(ISC,1)
      CALL FOXDAL(IFL,1)
      EPS = EPS*10.D0
*
      RETURN
      END
*
      SUBROUTINE DALEX2(INB,INA,INC,ND)
*     *********************************
*
*     THIS SUBROUTINE PERFORMS A LIE EXPONENTIATION OF THE FIRST VECTOR
*     ON THE SECOND VECTOR.
*     THE RESULT IS STORED IN INC. 2*ND IS THE DIMENSIONALITY OF PHASE SPACE.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION IBNO(0:LNO,LNV),IENO(0:LNO,LNV)
      COMMON /MONDA/ ID1(LNV),ID2(LNV)
      INTEGER JJ(LNV)
*
      IF(LEW.EQ.1) THEN
         PRINT*,'$$$ ERROR IN DALEX2, '//
     *          'WEIGHTED DA (DANOTW) IS NOT SUPPORTED'
         CALL FOXDEB
      ENDIF
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NDA) CALL FOXNTY(INB)
*
      CALL DACLR
*
      JMEM = IMEM
*
      DO 20 J=1,2*ND
      DO 20 I=0,NOCUT
      IBNO(I,J) = JMEM+1
      IENO(I,J) = JMEM
  20  JMEM = JMEM+NMMAX
*
      CALL MEMCHK(JMEM)
*
*     DIFFERENTIATING B BY ORDER
*     **************************
*
      DO 40 IB=NBEG(INB),NEND(INB)
*
      NCIB = NC(IB)
      NOIB = IEO(NCIB)-1
      CALL DAENC(IE1(NCIB),IE2(NCIB),JJ)
*
      JFAC = 1
      DO 35 J=1,2*ND
      JFAC = -JFAC
      JDFA = JJ(J)
      IF(JDFA.EQ.0) GOTO 35
      JS = J-JFAC
      IPOS = IENO(NOIB,JS)+1
      IENO(NOIB,JS) = IPOS
*
      CC(IPOS) = JFAC*JDFA*CC(IB)
      NC(IPOS) = IA1(IE1(NCIB)-ID1(J))+IA2(IE2(NCIB)-ID2(J))
  35  CONTINUE
  40  CONTINUE
*
*
*     PERFORMING MULTIPLICATIONS
*     **************************
*
      DO 100 IA=NBEG(INA),NEND(INA)
*
      NCIA = NC(IA)
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
      DO 100 J=1,2*ND
*
      JDFA = JJ(J)
      IF(JDFA.EQ.0) GOTO 100
      CCIA = JDFA*CC(IA)
      I1IA = IE1(NCIA)-ID1(J)
      I2IA = IE2(NCIA)-ID2(J)
*
      DO 90 NOIB = 0,NOCUT-IEO(NCIA)+1
      DO 90 IB = IBNO(NOIB,J),IENO(NOIB,J)
*
      NCIB = NC(IB)
      ICC = IA2(I2IA+IE2(NCIB))+IA1(I1IA+IE1(NCIB))
      CDA(ICC) = CDA(ICC)+CCIA*CC(IB)
*
  90  CONTINUE
 100  CONTINUE
*
      CALL DAPAC(INC)
*
      RETURN
      END
*
*%%
      SUBROUTINE DAPRI(INA,IUNIT)
*     ****************************
*
*     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CHARACTER LIN*80
      INTEGER JJ(LNV)
*
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
*
      IOUT = 0
      TMS = 0.D0
*
      DO 100 IOA = 0,NOMAX
      DO 100 II=NBEG(INA),NEND(INA)
      NCIA = NC(II)
      IF(IEO(NCIA).NE.IOA) GOTO 100
      CCC = CC(II)
      IF(ABS(CCC).LT.EPS) THEN
         TMS = TMS+ABS(CCC)
         GOTO 100
      ENDIF
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
*
      IOUT = IOUT+1
      IF(IOUT.EQ.1) THEN
      WRITE(IUNIT,'(A)')'     I  COEFFICIENT            ORDER EXPONENTS'
      ENDIF
*
      WRITE(IUNIT,'(I6,2X,G22.16,I4,2X,8(2I2,1X))')
     *   IOUT,CCC,IOA,(JJ(III),III=1,MIN(NVMAX,16))
      IF(NVMAX.GT.16) THEN
         DO 10 IIC=17,NVMAX,16
         WRITE(IUNIT,'(36X,8(2I2,1X))')
     *                (JJ(III),III=IIC,MIN(NVMAX,IIC+15))
 10      CONTINUE
      ENDIF
*
 100  CONTINUE
*
      IF(IOUT.EQ.0) THEN
         WRITE(IUNIT,'(5X,A)') 'ALL COMPONENTS ZERO '
         WRITE(IUNIT,'(5X,A)') '------------------- '
         RETURN
      ENDIF
*
      J = 5-INT(LOG(FLOAT(IOUT))/LOG(10.D0))
      DO 110 I=1,J
 110  LIN(I:I) = ' '
      N = 35+INT(FLOAT(MIN(NVMAX,16))*2.5E0+.5E0)
      DO 120 I=J+1,N
 120  LIN(I:I) = '-'
      WRITE(IUNIT,'(A)') LIN(1:N)
*
      RETURN
      END
*%%
*
      SUBROUTINE DAPRV(INA,NNA,NMP,NNP,INUNIT)
*     ****************************************
*
*     THIS SUBROUTINE PRINTS THE NNA DA VECTORS INA SIDE BY SIDE
*     TO UNIT IUNIT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CHARACTER FORM*47,HEXA*16,COUT*20
      DIMENSION MA(LNV),NN(LNV),CCN(LNV)
      INTEGER JJ(LNV)
      DATA HEXA / '0123456789ABCDEF' /
*
      IF(NTYP(NNA).NE.NRE) CALL FOXNTY(NNA)
      IF(NTYP(NNP).NE.NRE) CALL FOXNTY(NNP)
      IF(NTYP(NMP).NE.NRE) CALL FOXNTY(NMP)
      IF(NTYP(INUNIT).NE.NRE) CALL FOXNTY(INUNIT)
      IF(INA.LT.1.OR.INA.GT.IVAR) THEN
         PRINT*,'$$$ ERROR IN DAPRV, VARIABLE NOT ALLOCATED: ',INA
         CALL FOXDEB
      ELSEIF(NVMAX.GT.20) THEN
         PRINT*,'$$$ ERROR IN DAPRV, NVMAX>20 IS NOT SUPPORTED'
         CALL FOXDEB
      ENDIF
*
      NA    = NINT(CC(NBEG(NNA)))
      NP    = NINT(CC(NBEG(NNP)))
      MP    = NINT(CC(NBEG(NMP)))
      IUNITO = NINT(CC(NBEG(INUNIT)))
      CALL PRDIRB(IUNITO,IUNIT)
*
      GOTO(1,2,3,4,5), NA
*
      PRINT*,'$$$ ERROR IN DAPRV, '//
     *       'ILLEGAL NUMBER OF VECTORS (2ND ARGUMENT): ',NA
      PRINT*,'    BOUND = 5'
      CALL FOXDEB
*
   1  FORM = '(1X, G20.13,1X,A)'
      GOTO 6
   2  FORM = '(1X,2G20.13,1X,A)'
      GOTO 6
   3  FORM = '(1X,3G20.13,1X,A)'
      GOTO 6
   4  FORM = '(1X,4G16. 9,1X,A)'
      GOTO 6
   5  FORM = '(1X,5G14. 7,1X,A)'
      GOTO 6
   6  CONTINUE
*
      DO 10 I=1,NA
      MA(I) = INA+I
      IF(NTYP(INA+I).NE.NDA) CALL FOXNTY(INA+I)
  10  CONTINUE
*
      DO 100 IOA = 0,NOMAX
*
      DO 50 N=1,NA
      DO 30 J=NBEG(MA(N)),NEND(MA(N))
      IF(IEO(NC(J)).EQ.IOA) GOTO 40
  30  CONTINUE
  40  NN(N) = J
  50  CONTINUE
*
  60  CONTINUE
      NNC = LEA+1
      DO 70 N=1,NA
      IF(NN(N).GT.NEND(MA(N))) GOTO 70
      NNC = MIN(NNC,NC(NN(N)))
  70  CONTINUE
      IF(NNC.EQ.LEA+1) GOTO 100
*
      DO 90 N=1,NA
      IF(NN(N).GT.NEND(MA(N))) THEN
         CCN(N) = 0.D0
      ELSEIF(NC(NN(N)).EQ.NNC) THEN
         CCN(N) = CC(NN(N))
         DO 80 J=NN(N)+1,NEND(MA(N))
         IF(IEO(NC(J)).EQ.IOA) GOTO 85
  80     CONTINUE
  85     NN(N) = J
      ELSE
         CCN(N) = 0.D0
      ENDIF
  90  CONTINUE
*
      CALL DAENC(IE1(NNC),IE2(NNC),JJ)
*
      DO 95 III=1,NP
      IF(JJ(III).GT.15) THEN
         PRINT*,'$$$ ERROR IN DAPRV, EXPONENTS TOO LARGE'
         CALL FOXDEB
      ELSE
         IOUT = JJ(III)+1
         COUT(III:III) = HEXA(IOUT:IOUT)
      ENDIF
  95  CONTINUE
      DO 96 III=NP+1,MP
  96  COUT(III:III) = '0'
      DO 97 III=MP+1,NVMAX-NP+MP
      IOUT = JJ(III-MP+NP)+1
  97  COUT(III:III) = HEXA(IOUT:IOUT)
      WRITE(IUNIT,FORM) (CCN(N),N=1,NA),COUT(1:NVMAX-NP+MP)
      GOTO 60
*
 100  CONTINUE
*
      WRITE(IUNIT,'(1X,78A)') ('-',I=1,78)
*
      CALL PRDIRE(IUNITO,IUNIT)
*
      RETURN
      END
*
      SUBROUTINE DAREV(INA,NNA,NNM,NNI,IIUNIT)
*     ****************************************
*
*     THIS SUBROUTINE READS THE NNA DA VECTORS INA SIDE BY SIDE
*     FROM UNIT IUNIT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CHARACTER FORM*47,HEXA*16,COUT*20
      DIMENSION C(LEA,5),NN(LEA),J(LNV)
      DATA HEXA / '0123456789ABCDEF' /
*
      IF(NTYP(NNA).NE.NRE) CALL FOXNTY(NNA)
      IF(NTYP(NNI).NE.NRE) CALL FOXNTY(NNI)
      IF(NTYP(NNM).NE.NRE) CALL FOXNTY(NNM)
      IF(NTYP(IIUNIT).NE.NRE) CALL FOXNTY(IIUNIT)
      IF(INA.LT.1.OR.INA.GT.IVAR) THEN
         PRINT*,'$$$ ERROR IN DAREV, VARIABLE NOT ALLOCATED: ',INA
         CALL FOXDEB
      ENDIF
*
      NA    = NINT(CC(NBEG(NNA)))
      NI    = NINT(CC(NBEG(NNI)))
      NM    = NINT(CC(NBEG(NNM)))
      IUNIT = NINT(CC(NBEG(IIUNIT)))
*
      ILIN = 0
      IF(NA.NE.5) THEN
         PRINT*,'$$$ ERROR IN DAREV, NA DIFFERS FROM 5'
         CALL FOXDEB
      ELSEIF(NVMAX.GT.20) THEN
         PRINT*,'$$$ ERROR IN DAREV, NVMAX>20 IS NOT SUPPORTED'
         CALL FOXDEB
      ENDIF
      FORM = '(1X,5G14. 7,1X,A)'
*
  50  ILIN = ILIN+1
      READ(IUNIT,FORM,ERR=100) (C(ILIN,K),K=1,NA),COUT
      IO = 0
      DO 55 I=1,NVMAX+NM-NI
      IF(I.LE.NI) THEN
         J(I) = MAX(0,INDEX(HEXA,COUT(I:I))-1)
         IO = IO+J(I)
      ELSEIF(I.GT.NI.AND.I.LE.NM) THEN
         IF(INDEX(HEXA,COUT(I:I))-1.NE.0) THEN
            ILIN = ILIN-1
            GOTO 50
         ENDIF
      ELSE
         J(I+NI-NM) = MAX(0,INDEX(HEXA,COUT(I:I))-1)
         IO = IO+J(I+NI-NM)
      ENDIF
  55  CONTINUE
      DO 56 I=NVMAX+NM-NI+1,20
      IF(COUT(I:I).NE.' '.AND.COUT(I:I).NE.'0') THEN
         ILIN = ILIN-1
         GOTO 50
      ENDIF
  56  CONTINUE
      IF(IO.GT.NOCUT) THEN
         ILIN = ILIN-1
         GOTO 50
      ENDIF
      CALL DADEC(J,I1,I2)
      IF(I1.LT.0.OR.I2.LT.0) THEN
         ILIN = ILIN-1
         GOTO 50
      ENDIF
      NN(ILIN) = IA1(I1)+IA2(I2)
      IF(NN(ILIN).GT.NMMAX.OR.NN(ILIN).LT.1) THEN
         PRINT*,'@@@ ERROR IN DAREV, NN(ILIN) = ',NN(ILIN)
      ELSEIF(IO.NE.IEO(NN(ILIN))) THEN
         PRINT*,'@@@ ERROR IN DAREV, IO CHECK FAILED'
      ENDIF
      GOTO 50
*
  100 CONTINUE
      ILIN = ILIN-1
      DO 120 K=1,NA
      CALL DACLR
      DO 110 I=1,ILIN
  110 CDA(NN(I)) = C(I,K)
  120 CALL DAPAC(INA+K)
*
      RETURN
      END
*
      SUBROUTINE DAREA(IIUNIT,INA,INVA)
*     *********************************
*
*     THIS SUBROUTINE READS THE DA VECTOR INA FROM UNIT IUNIT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CHARACTER ALIN1*80,ALIN2*80,FORM1*60,FORM2*60
      DIMENSION J(LNV+100)
*
      IF(NTYP(IIUNIT).NE.NRE) CALL FOXNTY(IIUNIT)
      IF(NTYP(INVA).NE.NRE) CALL FOXNTY(INVA)
      IUNIT = NINT(CC(NBEG(IIUNIT)))
      NVA   = NINT(CC(NBEG(INVA)))
      IF(NVA.GT.LNV) THEN
         PRINT*,'$$$ ERROR IN DAREA, 3RD ARGUMENT > LNV'
         CALL FOXDEB
      ENDIF
*
      CALL DACLR
      DO 10 I=1,LNV
 10   J(I) = 0
*
      ITM = 0
      TMT = 0.D0
      TMS = 0.D0
*
 20   READ(IUNIT,'(A)',END=200) ALIN1
      IF(INDEX(ALIN1,'ALL COMPONENTS ZERO').NE.0) THEN
         READ(IUNIT,'(A)') ALIN1
         NTYP(INA) = NDA
         NEND(INA) = NBEG(INA)-1
         RETURN
      ENDIF
      IF(INDEX(ALIN1,'     I  COEFFICIENT ').EQ.0) GOTO 20
*
      FORM1 = '(I6,2X,G22.16,I4,2X,8(2I2,1X))'
      FORM2 = '(36X,8(2I2,1X))'
*
      READ(IUNIT,'(A)',END=300) ALIN1
      ILIN = 0
      NLIN = 1
      NVLIN = (ILAST(ALIN1,1,80)-35)*2/5
      READ(ALIN1,FORM1,ERR=300) II,C,IO,(J(I),I=1,NVLIN)
      NVAMAX = NVLIN
      IF(NVLIN.LT.16) GOTO 100
*
 30   READ(IUNIT,'(A)') ALIN2
      IF(ALIN2(1:34).EQ.'                                  ') THEN
         NLIN = NLIN+1
         NVLIN = (ILAST(ALIN2,1,80)-35)*2/5
         IF(NVAMAX+NVLIN.GT.LNV+100) GOTO 300
         READ(ALIN2,FORM2,ERR=300) (J(I),I=NVAMAX+1,NVAMAX+NVLIN)
         NVAMAX = NVAMAX+NVLIN
      ELSE
         ILIN = -1
         NVLIN = 0
         ALIN1 = ALIN2
      ENDIF
      IF(NVLIN.EQ.16) GOTO 30
*
 100  CONTINUE
      IIO = 0
      DO 40 I=NVA+1,NVAMAX
 40   IIO = IIO+J(I)
      IF(IIO.NE.0) GOTO 60
      IIO = 0
      DO 50 I=1,NVA
 50   IIO = IIO+J(I)
      IF(IIO.GT.NOCUT) GOTO 60
      CALL DADEC(J,II1,II2)
      IF(II1.LT.0.OR.II2.LT.0) GOTO 60
      IC = IA1(II1)+IA2(II2)
      CDA(IC) = C
      TMT = TMT+ABS(C)
*
 60   CONTINUE
      IF(ILIN.EQ.0) READ(IUNIT,'(A)') ALIN1
      IF(INDEX(ALIN1,'----------').NE.0) GOTO 200
*
      ILIN = 0
      NVLIN = (ILAST(ALIN1,1,80)-35)*2/5
      READ(ALIN1,FORM1,ERR=300) II,C,IO,(J(I),I=1,NVLIN)
      NVAMAX = NVLIN
*
      IF(NLIN.GE.2) THEN
         DO 70 IL=2,NLIN
         READ(IUNIT,'(A)',END=300) ALIN2
         NVLIN = (ILAST(ALIN2,1,80)-35)*2/5
         READ(ALIN2,FORM2,ERR=300) (J(I),I=NVAMAX+1,NVAMAX+NVLIN)
         NVAMAX = NVAMAX+NVLIN
 70      CONTINUE
      ENDIF
      GOTO 100
*
 200  CALL DAPAC(INA)
      RETURN
*
 300  CONTINUE
      PRINT*,'$$$ ERROR IN DAREA, DATA IS WRONG.'
      CALL FOXDEB
*
      END
*%%
*
      SUBROUTINE DAPEW(INUNIT,INA,INIA,INOR)
*     **************************************
*
*     THIS SUBROUTINE PRINTS THE PART OF DA VECTOR INA THAT HAS ORDER INOR
*     IN THE INIA-TH COLUMN.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CHARACTER LIN*80
      INTEGER JJ(LNV)
*
      IF(NTYP(INA).NE.NDA) THEN
         CALL FOXNTY(INA)
      ELSEIF(NTYP(INUNIT).NE.NRE) THEN
         CALL FOXNTY(INUNIT)
      ELSEIF(NTYP(INIA).NE.NRE) THEN
         CALL FOXNTY(INIA)
      ELSEIF(NTYP(INOR).NE.NRE) THEN
         CALL FOXNTY(INOR)
      ENDIF
*
      IIA   = NINT(CC(NBEG(INIA)))
      IPOS  = NINT(CC(NBEG(INOR)))
      IUNITO = NINT(CC(NBEG(INUNIT)))
      CALL PRDIRB(IUNITO,IUNIT)
*
*     WRITE(IUNIT,'(A)') ' '
*     WRITE(IUNIT,'(2A)')  ' DA VECTOR'
*     WRITE(IUNIT,'(A)')   ' *********'
*
      WRITE(IUNIT,'(8X,A5,I4,A11,I5)') 'ORDER',IPOS,'  IN COLUMN',IIA
      IOUT = 0
      DO 100 IOA = IPOS,MIN(NOCUT+IPOS,NOMAX)
      DO 100 II=NBEG(INA),NEND(INA)
      NCIA = NC(II)
      IF(IEO(NCIA).NE.IOA) GOTO 100
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
      IF(JJ(IIA).NE.IPOS) GOTO 100
      IF(ABS(CC(II)).LT.EPS) GOTO 100
      IOUT = IOUT+1
      IF(IOUT.EQ.1) THEN
         WRITE(IUNIT,'(5X,A)') 'I  COEFFICIENT          ORDER EXPONENTS'
      ENDIF
      WRITE(IUNIT,'(I6,2X,G22.16,I4,2X,8(2I2,1X),20(/34X,8(2I2,1X)))')
     * IOUT,CC(II),IOA,(JJ(III),III=1,IIA-1),0,
     * (JJ(III),III=IIA+1,NVMAX)
*
 100  CONTINUE
*
      IF(IOUT.EQ.0) THEN
         WRITE(IUNIT,'(5X,A)') 'ALL COMPONENTS ZERO '
         WRITE(IUNIT,'(5X,A)') '------------------- '
*
         CALL PRDIRE(IUNITO,IUNIT)
         RETURN
      ENDIF
*
      J = 5-INT(LOG(FLOAT(IOUT))/LOG(10.D0))
      DO 110 I=1,J
 110  LIN(I:I) = ' '
      N = 35+INT(FLOAT(NVMAX)*2.5E0+.5E0)
      DO 120 I=J+1,N
 120  LIN(I:I) = '-'
      WRITE(IUNIT,'(A)') LIN(1:N)
*
      CALL PRDIRE(IUNITO,IUNIT)
*
      RETURN
      END
*
      SUBROUTINE DARAN(INA,ICM)
*     *************************
*
*     THIS SUBROUTINE FILLS THE DA VECTOR A WITH RANDOM ENTRIES.
*     FOR CM > 0, THE VECTOR IS FILLED WITH REALS,
*     FOR CM < 0, THE VECTOR IS FILLED WITH SINGLE DIGIT INTEGERS
*     ABS(CM) IS THE FILLING FACTOR
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      SAVE IRAN
      DATA IRAN /1234321/
*
      IF(NTYP(ICM).NE.NRE) CALL FOXNTY(ICM)
      CM = CC(NBEG(ICM))
*
      DO 100 I=1,NMMAX
      IF(CM.GE.0.D0) THEN
         CDA(I) = 2*(BRAN(IRAN)-.5)
         IF(ABS(CDA(I)).GT.CM) CDA(I) = 0.D0
      ELSEIF(CM.LT.0.D0) THEN
         CDA(I) = INT(20*BRAN(IRAN)-10)
         IF(ABS(CDA(I)).GT.-10.D0*CM) CDA(I) = 0.D0
      ENDIF
 100  CONTINUE
*
      CALL DAPAC(INA)
*
      RETURN
      END
*
      SUBROUTINE RERAN(INA)
*     *********************
*
*     THIS SUBROUTINE FILLS THE REAL VECTOR A WITH A RANDOM NUMBER BETWEEN -1,1
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      SAVE IRAN
      DATA IRAN /1234321/
*
      CC(NBEG(INA)) = 2.D0*(BRAN(IRAN)-.5D0)
      NEND(INA) = NBEG(INA)-1
      NTYP(INA) = NRE
*
      RETURN
      END
*
      DOUBLE PRECISION FUNCTION DBFLOAT(I)
*     ************************************
*
*     PRECISE FLOAT TO DOUBLE PRECISION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(0:100)
*
      DATA F / 0.D0,
     *    1.D0, 2.D0, 3.D0, 4.D0, 5.D0, 6.D0, 7.D0, 8.D0, 9.D0,10.D0,
     *   11.D0,12.D0,13.D0,14.D0,15.D0,16.D0,17.D0,18.D0,19.D0,20.D0,
     *   21.D0,22.D0,23.D0,24.D0,25.D0,26.D0,27.D0,28.D0,29.D0,30.D0,
     *   31.D0,32.D0,33.D0,34.D0,35.D0,36.D0,37.D0,38.D0,39.D0,40.D0,
     *   41.D0,42.D0,43.D0,44.D0,45.D0,46.D0,47.D0,48.D0,49.D0,50.D0,
     *   51.D0,52.D0,53.D0,54.D0,55.D0,56.D0,57.D0,58.D0,59.D0,60.D0,
     *   61.D0,62.D0,63.D0,64.D0,65.D0,66.D0,67.D0,68.D0,69.D0,70.D0,
     *   71.D0,72.D0,73.D0,74.D0,75.D0,76.D0,77.D0,78.D0,79.D0,80.D0,
     *   81.D0,82.D0,83.D0,84.D0,85.D0,86.D0,87.D0,88.D0,89.D0,90.D0,
     *   91.D0,92.D0,93.D0,94.D0,95.D0,96.D0,97.D0,98.D0,99.D0,100.D0/
*
      IF(I.LE.100) THEN
         IF(I.GE.0) THEN
            DBFLOAT = F(I)
            RETURN
         ELSEIF(I.GE.-100) THEN
            DBFLOAT = -F(-I)
            RETURN
         ELSE
            DBFLOAT = DBLE(FLOAT(I))
            RETURN
         ENDIF
      ELSE
         DBFLOAT = DBLE(FLOAT(I))
         RETURN
      ENDIF
*
      END
*
      DOUBLE PRECISION FUNCTION BRAN(IRAN)
*     ************************************
*
*     RANDOM NUMBER GENERATOR
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      SAVE XRAN
      DATA XRAN / .234D0 /
*
      XRAN = XRAN+10.D0+IRAN/1.D6
      IF(XRAN.GT.1.D4) XRAN = XRAN-9999*ABS(COS(XRAN))
      BRAN = ABS(SIN(XRAN))
      BRAN = 1E6*BRAN
      BRAN = BRAN-INT(BRAN)
*
      RETURN
      END
*
*%%
      SUBROUTINE DANUM(NO,NV,NUMDA)
*     *****************************
*
*     THIS SUBROUTINE COMPUTES THE NUMBER OF MONOMIALS UP TO
*     ORDER NO AND NUMBER OF VARIABLES NV
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
*
      NUMDA = 1
      MM = MAX(NV,NO)
*
      DO 5 I=1,MIN(NV,NO)
  5   NUMDA = (NUMDA*(MM+I))/I
*
      RETURN
      END
*
      SUBROUTINE DANUMD(NO,NV,CNUMDA)
*     *******************************
*
*     THIS SUBROUTINE COMPUTES THE NUMBER OF MONOMIALS UP TO
*     ORDER NO AND NUMBER OF VARIABLES NV IN DOUBLE PRECISION.
*     NO AND NV ARE ASSUMED TO BE POSITIVE AND LESS THAN 2147483647.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
*
      CNUMDA = 1.D0
      CMM = DBLE(MAX(NV,NO))
*
      DO 5 I=1,MIN(NV,NO)
  5   CNUMDA = (CNUMDA*(CMM+DBLE(I)))/DBLE(I)
*
      RETURN
      END
*%%
      SUBROUTINE DAFSET(INA)
*     **********************
*
*     THIS SUBROUTINE SETS THE DA FILTERING MODE.
*     A TEMPLATE DA IS SUPPLIED BY INA.
*     IF INA IS 0 OR DAINI IS CALLED AHEAD, THE MODE IS TURNED OFF.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).EQ.NRE) THEN
         IF(NINT(CC(NBEG(INA))).EQ.0) LFLT = 0
      ELSEIF(NTYP(INA).EQ.NDA) THEN
         LFLT = 1
         NFLT = 0
         DO 10 I=NBEG(INA),NEND(INA)
         NFLT = NFLT+1
         NCFLT(NFLT) = NC(I)
 10      CONTINUE
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DAFILT(INA,INB)
*     **************************
*
*     THIS SUBROUTINE FILTERS THE DA OR CD VECTOR A TO B
*     THROUGH THE DA FILTERING TEMPLATE SUPPLIED BY DAFSET.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IAMAX = NEND(INA)
      IA = NBEG(INA)
      NA = NC(IA)
*
      IF(NTYP(INA).EQ.NDA) THEN
         IB = NBEG(INB)-1
         IF(IA.GT.IAMAX) GOTO 150
*
         DO 140 I=1,NFLT
         NCF = NCFLT(I)
 110     IF(NA-NCF) 120,130,140
*
 120     CONTINUE
         IF(IA.EQ.IAMAX) GOTO 150
         IA = IA+1
         NA = NC(IA)
         GOTO 110
*
 130     CONTINUE
         IB = IB+1
         CC(IB) = CC(IA)
         NC(IB) = NA
         IF(IA.EQ.IAMAX) GOTO 150
         IA = IA+1
         NA = NC(IA)
 140     CONTINUE
*
 150     CONTINUE
         NTYP(INB) = NDA
         NEND(INB) = IB
*
      ELSEIF(NTYP(INA).EQ.NCD) THEN
         IB = NBEG(INB)-2
         IF(IA.GT.IAMAX) GOTO 250
*
         DO 240 I=1,NFLT
         NCF = NCFLT(I)
 210     IF(NA-NCF) 220,230,240
*
 220     CONTINUE
         IF(IA.EQ.IAMAX) GOTO 250
         IA = IA+2
         NA = NC(IA)
         GOTO 210
*
 230     CONTINUE
         IB = IB+2
         CC(IB  ) = CC(IA)
         CC(IB+1) = CC(IA+1)
         NC(IB  ) = NA
         NC(IB+1) = 0
         IF(IA.EQ.IAMAX) GOTO 250
         IA = IA+2
         NA = NC(IA)
 240     CONTINUE
*
 250     CONTINUE
         NTYP(INB) = NCD
         NEND(INB) = IB+1
*
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('DAFILT')
*
      RETURN
      END
*
      SUBROUTINE DACLR
*     ****************
*
*     THIS SUBROUTINE SETS ALL THE SCRATCH SPACE IN CDA TO ZERO
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(LFLT.EQ.0) THEN
         DO 10 I=1,NMMAX
  10     CDA(I) = 0.D0
      ELSE
         DO 20 I=1,NFLT
  20     CDA(NCFLT(I)) = 0.D0
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DAPAC(INC)
*     *********************
*
*     THIS SUBROUTINE PACKS THE INFORMATION IN THE SCRATCH VECTOR CDA
*     INTO THE VECTOR INC. INVERSE IS DAUNP.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IC = NBEG(INC)-1
*
      IF(LFLT.EQ.0) THEN
         IF(ITM.NE.1) THEN
            DO 100 I=1,NMMAX
            CCC = CDA(I)
            IF(ABS(CCC).LT.EPS) GOTO 100
            IC = IC+1
            CC(IC) = CCC
            NC(IC) = I
 100        CONTINUE
         ELSEIF(ITM.EQ.1) THEN
            DO 110 I=1,NMMAX
            CCC = CDA(I)
            IF(ABS(CCC).LT.EPS) THEN
               TMS = TMS+ABS(CCC)
               GOTO 110
            ENDIF
            IC = IC+1
            CC(IC) = CCC
            NC(IC) = I
 110        CONTINUE
         ENDIF
      ELSE
         IF(ITM.NE.1) THEN
            DO 200 I=1,NFLT
            NCF = NCFLT(I)
            CCC = CDA(NCF)
            IF(ABS(CCC).LT.EPS) GOTO 200
            IC = IC+1
            CC(IC) = CCC
            NC(IC) = NCF
 200        CONTINUE
         ELSEIF(ITM.EQ.1) THEN
            DO 210 I=1,NFLT
            NCF = NCFLT(I)
            CCC = CDA(NCF)
            IF(ABS(CCC).LT.EPS) THEN
               TMS = TMS+ABS(CCC)
               GOTO 210
            ENDIF
            IC = IC+1
            CC(IC) = CCC
            NC(IC) = NCF
 210        CONTINUE
         ENDIF
      ENDIF
*
      NTYP(INC) = NDA
      NEND(INC) = IC
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAPAC')
*
      RETURN
      END
*%%
      SUBROUTINE DAUNP(INC)
*     *********************
*
*     THIS SUBROUTINE UNPACKS THE VECTOR C.
*     INVERSE IS DAPAC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL DACLR
*
      DO 10 I=NBEG(INC),NEND(INC)
      IC = NC(I)
      CDA(IC) = CC(I)
  10  CONTINUE
*
      RETURN
      END
*%%
      INTEGER FUNCTION IJJ(J,N)
*     *************************
*
*     THIS FUNCTION EXTRACTS THE EXPONENT OF THE J-TH VARIABLE IN THE N-TH TERM
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(J.GT.NVMAX) THEN
         IJJ = 0
         RETURN
      ELSEIF(J.LE.IESP) THEN
         IE = IE1(N)
         I1 = 1
      ELSE
         IE = IE2(N)
         I1 = IESP+1
      ENDIF
*
      IF(LEW.EQ.0) THEN
         IBASE = NOMAX+1
         IJJ=MOD(IE/IED(J),IBASE)
      ELSE
         DO 10 I=I1,J
         IBASE = NOMAX/IEW(I)+1
         X  = IE/DBFLOAT(IBASE)
         IE = INT(X+EPSMAC)
  10     IJJ = NINT(IBASE*(X-IE))*IEW(I)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DADEC(JJ,IC1,IC2)
*     ****************************
*
*     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES I1,I2.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      INTEGER JJ(LNV)
*
      IC1 = 0
      IC2 = 0
*
      IF(LEW.EQ.0) THEN
         IBASE = NOMAX+1
*
         DO 10 I=NVMAX,IESP+1,-1
 10      IC2 = IC2*IBASE+JJ(I)
*
         DO 20 I=IESP,1,-1
 20      IC1 = IC1*IBASE+JJ(I)
      ELSE
         DO 30 I=NVMAX,IESP+1,-1
         IF(MOD(JJ(I),IEW(I)).NE.0) THEN
            IC2 = -1
            GOTO 40
         ENDIF
 30      IC2 = IC2*(NOMAX/IEW(I)+1)+JJ(I)/IEW(I)
 40      CONTINUE
*
         DO 50 I=IESP,1,-1
         IF(MOD(JJ(I),IEW(I)).NE.0) THEN
            IC1 = -1
            GOTO 60
         ENDIF
 50      IC1 = IC1*(NOMAX/IEW(I)+1)+JJ(I)/IEW(I)
 60      CONTINUE
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DAENC(IC1,IC2,JJ)
*     ****************************
*
*     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES I1,I2.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      INTEGER JJ(LNV)
*
      IF(LEW.EQ.0) THEN
         IBASE = NOMAX+1
*
         IC = IC1
         DO 10 I=1,IESP
         X  = IC/DBFLOAT(IBASE)
         IC = INT(X+EPSMAC)
  10     JJ(I) = NINT(IBASE*(X-IC))
*
         IC = IC2
         DO 20 I=IESP+1,NVMAX
         X  = IC/DBFLOAT(IBASE)
         IC = INT(X+EPSMAC)
  20     JJ(I) = NINT(IBASE*(X-IC))
      ELSE
         IC = IC1
         DO 30 I=1,IESP
         IBASE = NOMAX/IEW(I)+1
         X  = IC/DBFLOAT(IBASE)
         IC = INT(X+EPSMAC)
  30     JJ(I) = NINT(IBASE*(X-IC))*IEW(I)
*
         IC = IC2
         DO 40 I=IESP+1,NVMAX
         IBASE = NOMAX/IEW(I)+1
         X  = IC/DBFLOAT(IBASE)
         IC = INT(X+EPSMAC)
  40     JJ(I) = NINT(IBASE*(X-IC))*IEW(I)
      ENDIF
*
      DO 50 I=NVMAX+1,LNV
  50  JJ(I) = 0
*
      RETURN
      END
*
      SUBROUTINE DAENCW(IC1,IC2,JJ)
*     *****************************
*
*     THIS SUBROUTINE ENCODES THE WEIGHT DIVIDED EXPONENTS IN JJ FROM
*     THEIR DA CODES I1,I2.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      INTEGER JJ(LNV)
*
      IC = IC1
      DO 10 I=1,IESP
      IBASE = NOMAX/IEW(I)+1
      X  = IC/DBFLOAT(IBASE)
      IC = INT(X+EPSMAC)
  10  JJ(I) = NINT(IBASE*(X-IC))
*
      IC = IC2
      DO 20 I=IESP+1,NVMAX
      IBASE = NOMAX/IEW(I)+1
      X  = IC/DBFLOAT(IBASE)
      IC = INT(X+EPSMAC)
  20  JJ(I) = NINT(IBASE*(X-IC))
*
      DO 30 I=NVMAX+1,LNV
  30  JJ(I) = 0
*
      RETURN
      END
*%%
      SUBROUTINE DAVAR(INA,CKON,I)
*     ****************************
*
*     THIS SUBROUTINE DECLARES THE DA VECTOR
*     AS THE INDEPENDENT VARIABLE NUMBER I WITH VALUE CKON
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      INTEGER JJ(LNV)
*
      IPOA = NBEG(INA)
*
      IF(IEW(I).GT.NOMAX) THEN
         NEND(INA) = IPOA-1
         NTYP(INA) = NDA
         TMT = 0.D0
         TMS = 0.D0
         RETURN
      ENDIF
*
      DO 15 J=1,LNV
  15  JJ(J) = 0
      JJ(I) = IEW(I)
*
      CALL DADEC(JJ,IC1,IC2)
*
      NCA = IA1(IC1)+IA2(IC2)
*
      TMT = 1.D0
      TMS = 0.D0
*
      IF(ABS(CKON).GT.EPS) THEN
         TMT = TMT+ABS(CKON)+MAX(ABS(CKON),1.D0)
         NEND(INA) = IPOA+1
         CC(IPOA) = CKON
         NC(IPOA) = 1
*
         CC(IPOA+1) = 1.D0
         NC(IPOA+1) = NCA
         IF(NMAX(INA).LT.NBEG(INA)+2) CALL FOXERV('DAVAR')
         NTYP(INA) = NDA
      ELSE
         TMS = ABS(CKON)
         NEND(INA) = IPOA
         CC(IPOA) = 1.D0
         NC(IPOA) = NCA
         NTYP(INA) = NDA
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DACON(INA,CKON)
*     **************************
*
*     THIS SUBROUTINE SETS THE VECTOR C TO THE CONSTANT CKON
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IPOA = NBEG(INA)
      NTYP(INA) = NDA
      NEND(INA) = IPOA
      CC(IPOA) = CKON
      NC(IPOA) = 1
      IF(ABS(CKON).LT.EPS) NEND(INA) = IPOA-1
*
      RETURN
      END
*
      SUBROUTINE DAPEE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE DETERMINES THE COEFFICIENT OF THE DA, CD, TM BY
*     SPECIFYING THE EXPONENTS BY TRANSPORT NOTATION NUMBER IN INB
*
*     THIS SUBROUTINE PROVIDES AN INTERFACE TO DAPEA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION JJ(LNV)
*
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      CCC = CC(NBEG(INB))
*
      IF(CCC+1.D-8.LT.0.D0) THEN
         PRINT*,'$$$ ERROR IN DAPEE, ID IS NEGATIVE: ',CCC
         CALL FOXDEB
      ENDIF
*
      DO 10 J=1,LNV
  10  JJ(J) = 0
*
  20  CONTINUE
      CCC = (CCC/10.D0)
      CNEW = INT(CCC+1.D-2)
      J = NINT(10.D0*(CCC-CNEW))
      IF(J.GT.0.AND.J.LE.LNV) JJ(J) = JJ(J)+1
      CCC = CNEW
      IF(CCC.GT..5D0) GOTO 20
*
      CALL FOXAAL(ISJJ,LNV,1)
      CALL FOXALL(ISR,1,1)
*
      NTYP(ISR) = NRE
      CC(NBEG(ISR)) = LNV
*
      DO 30 J=1,LNV
      NTYP(ISJJ+J) = NRE
      CC(NBEG(ISJJ+J)) = JJ(J)
  30  CONTINUE
*
      CALL DAPEA(INA,ISJJ,ISR,INC)
*
      CALL FOXDAL(ISR,1)
      CALL FOXADA(ISJJ,LNV)
*
      RETURN
      END
*
      SUBROUTINE DAPEA(INA,IMB,INB,INC)
*     *********************************
*
*     THIS SUBROUTINE DETERMINES THE COEFFICIENT OF THE DA, CD, TM BY
*     SPECIFYING THE EXPONENTS IN AN ARRAY JJ FROM IMB WITH LENGTH INB
*
*     THIS SUBROUTINE PROVIDES AN INTERFACE TO DAPEK ALLOWING HIGHER ORDER
*     AND HIGHER DIMENSIONAL CASE.        (CF DAPEE)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION JJ(LNV)
*
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      NV = NINT(CC(NBEG(INB)))
*
      IF(NV.LE.0) THEN
         PRINT*,'$$$ ERROR IN DAPEA, 3RD ARGUMENT IS NOT POSITIVE'
         CALL FOXDEB
      ENDIF
*
      NV = MIN(NV,LNV)
*
      IO = 0
      DO 10 J=1,LNV
  10  JJ(J) = 0
*
      DO 20 J=1,NV
      IF(NTYP(IMB+J).NE.NRE) CALL FOXNTY(IMB+J)
      JJ(J) = NINT(CC(NBEG(IMB+J)))
      IF(JJ(J).LT.0) THEN
         PRINT*,'$$$ ERROR IN DAPEA, ARRAY ELEMENT (',J,') IS NAGATIVE'
         CALL FOXDEB
      ENDIF
      IO = IO+JJ(J)
  20  CONTINUE
*
      IF(IO.GT.NOMAX) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
      ELSEIF(NTYP(INA).EQ.NRE) THEN
         IF(IO.EQ.0) THEN
            CC(NBEG(INC)) = CC(NBEG(INA))
            NTYP(INC) = NRE
         ELSE
            CC(NBEG(INC)) = 0.D0
            NTYP(INC) = NRE
         ENDIF
C     ELSEIF(NTYP(INA).EQ.NDA.OR.NTYP(INA).EQ.NTM) THEN
      ELSEIF(NTYP(INA).EQ.NDA) THEN
         CALL DADEC(JJ,IC1,IC2)
         IF(IC1.LT.0.OR.IC2.LT.0) THEN
            CJJ = 0.D0
         ELSE
            NCA = IA1(IC1)+IA2(IC2)
            CALL DAPEK(INA,NCA,CJJ)
         ENDIF
*
         NTYP(INC) = NRE
         CC(NBEG(INC)) = CJJ
      ELSEIF(NTYP(INA).EQ.NCD) THEN
         CALL FOXALL(IAR,1,2*NMMAX)
         CALL FOXALL(IAI,1,2*NMMAX)
         CALL CDRE(INA,IAR)
         CALL CDIM(INA,IAI)
*
         CALL DADEC(JJ,IC1,IC2)
         IF(IC1.LT.0.OR.IC2.LT.0) THEN
            CJR = 0.D0
            CJI = 0.D0
         ELSE
            NCA = IA1(IC1)+IA2(IC2)
            CALL DAPEK(IAR,NCA,CJR)
            CALL DAPEK(IAI,NCA,CJI)
         ENDIF
*
         NTYP(INC) = NCM
         CC(NBEG(INC)  ) = CJR
         CC(NBEG(INC)+1) = CJI
         NEND(INC) = NBEG(INC)+1
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAPEA')
         CALL FOXDAL(IAI,1)
         CALL FOXDAL(IAR,1)
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DAPIC(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE DETERMINES THE COEFFICIENT OF THE DA, CD, TM BY
*     SPECIFYING THE EXPONENTS IN A VECTOR INB
*
*     THIS SUBROUTINE PROVIDES AN INTERFACE TO DAPEK.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION JJ(LNV)
*
      IF(NTYP(INB).NE.NRE.AND.NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IB = NBEG(INB)
      LB = NEND(INB)-IB+1
*
      IO = 0
      DO 10 J=1,LNV
 10   JJ(J) = 0
*
      DO 20 J=1,MIN(LB,LNV)
      JJ(J) = NINT(CC(IB))
      IF(JJ(J).LT.0) THEN
         PRINT*,'$$$ ERROR IN DAPIC, NEGATIVE EXPONENT IS GIVEN'
         CALL FOXDEB
      ENDIF
      IO = IO+JJ(J)
      IB = IB+1
 20   CONTINUE
*
      IF(IO.GT.NOMAX) THEN
         NTYP(INC) = NRE
         CC(NBEG(INC)) = 0.D0
      ELSEIF(NTYP(INA).EQ.NRE) THEN
         IF(IO.EQ.0) THEN
            CC(NBEG(INC)) = CC(NBEG(INA))
            NTYP(INC) = NRE
         ELSE
            CC(NBEG(INC)) = 0.D0
            NTYP(INC) = NRE
         ENDIF
C     ELSEIF(NTYP(INA).EQ.NDA.OR.NTYP(INA).EQ.NTM) THEN
      ELSEIF(NTYP(INA).EQ.NDA) THEN
         CALL DADEC(JJ,IC1,IC2)
         IF(IC1.LT.0.OR.IC2.LT.0) THEN
            CJJ = 0.D0
         ELSE
            NCA = IA1(IC1)+IA2(IC2)
            CALL DAPEK(INA,NCA,CJJ)
         ENDIF
*
         NTYP(INC) = NRE
         CC(NBEG(INC)) = CJJ
      ELSEIF(NTYP(INA).EQ.NCD) THEN
         CALL FOXALL(IAR,1,2*NMMAX)
         CALL FOXALL(IAI,1,2*NMMAX)
         CALL CDRE(INA,IAR)
         CALL CDIM(INA,IAI)
*
         CALL DADEC(JJ,IC1,IC2)
         IF(IC1.LT.0.OR.IC2.LT.0) THEN
            CJR = 0.D0
            CJI = 0.D0
         ELSE
            NCA = IA1(IC1)+IA2(IC2)
            CALL DAPEK(IAR,NCA,CJR)
            CALL DAPEK(IAI,NCA,CJI)
         ENDIF
*
         NTYP(INC) = NCM
         CC(NBEG(INC)  ) = CJR
         CC(NBEG(INC)+1) = CJI
         NEND(INC) = NBEG(INC)+1
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('DAPIC')
         CALL FOXDAL(IAI,1)
         CALL FOXDAL(IAR,1)
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DAPEK(INA,IC,CJJ)
*     ****************************
*
*     THIS SUBROUTINE DETERMINES THE COEFFICIENT OF THE DA INDEX NUMBER IC
*     AND RETURNS IT IN CJJ
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTY = NTYP(INA)
*
C     IF(NTY.EQ.NTM) THEN
C        NETM = NEND(INA)
C        NEDA = NBEG(INA)+NC(NETM)-1
C        NTYP(INA) = NDA
C        NEND(INA) = NEDA
C     ELSEIF(NTY.NE.NDA) THEN
      IF(NTY.NE.NDA) THEN
         CALL FOXNTY(INA)
      ENDIF
*
*     DETERMINE IF MONOMIAL IS INSIDE FIRST AND LAST MONOMIALS OF A
*     *************************************************************
*
      IU = NBEG(INA)
      IZ = NEND(INA)
      ICU = NC(IU)
      ICZ = NC(IZ)
*
      IF(NBEG(INA).EQ.NEND(INA)+1) THEN
         CJJ = 0.D0
         GOTO 100
      ELSEIF(IC.EQ.ICU) THEN
         CJJ = CC(IU)
         GOTO 100
      ELSEIF(IC.EQ.ICZ) THEN
         CJJ = CC(IZ)
         GOTO 100
      ELSEIF(IC.LT.ICU.OR.IC.GT.ICZ) THEN
         CJJ = 0.D0
         GOTO 100
      ENDIF
*
*     SEARCHING PROPER MONOMIAL
*     *************************
*
 10   CONTINUE
      IF(IZ-IU.LE.1) THEN
         CJJ = 0.D0
         GOTO 100
      ENDIF
      I = (IU+IZ)/2
*
      IF(NC(I)-IC) 20,30,40
 20   IU = I
      GOTO 10
 30   CJJ = CC(I)
      GOTO 100
 40   IZ = I
      GOTO 10
*
 100  CONTINUE
C     IF(NTY.EQ.NTM) THEN
C        NTYP(INA) = NTM
C        NEND(INA) = NETM
C     ENDIF
*
      RETURN
      END
*
      SUBROUTINE DACODE(IDAPR,INCL,IMCE)
*     **********************************
*
*     THIS SUBROUTINE DECODES THE DA INTERNAL MONOMIAL NUMBERS TO THE EXPONENTS.
*     THE INPUT DA SETUP IDAPR (VE) IS CHECKED AGAINST THE CURRENT DAINI SETUP.
*     THE EXPONENTS ARE RETURNED AS THE VE ARRAY IMCE WITH LENGTH INCL.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION JJ(LNV)
*
      CALL DACHK(IDAPR,'DACODE')
*
      IF(NTYP(INCL).NE.NRE) CALL FOXNTY(INCL)
      NCL = NINT(CC(NBEG(INCL)))
*
      IF(NTYP(IMCE).GT.0) THEN
         PRINT*,'$$$ ERROR IN DACODE, 3RD ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(NCL.LT.NMMAX) THEN
         PRINT*,'$$$ ERROR IN DACODE, '//
     *          'ARRAY LENGTH (2ND ARGUMENT) NOT ENOUGH. NEED',NMMAX
         CALL FOXDEB
      ELSEIF(NMAX(IMCE+1)-NBEG(IMCE+1)+1.LT.NVMAX) THEN
         CALL FOXERV('DACODE')
      ENDIF
*
      DO 20 I=1,NMMAX
      CALL DAENC(IE1(I),IE2(I),JJ)
      NTYP(IMCE+I) = NRE
      CC(NBEG(IMCE+I)) = JJ(1)
      IF(NVMAX.GE.2) THEN
         NTYP(IMCE+I) = NVE
         NEND(IMCE+I) = NBEG(IMCE+I)+NVMAX-1
         DO 10 J=2,NVMAX
         CC(NBEG(IMCE+I)+J-1) = JJ(J)
 10      CONTINUE
      ENDIF
 20   CONTINUE
*
      RETURN
      END
*%%
      BLOCKDATA POLDT
*     ***************
*
*     THIS BLOCKDATA CONTAINS AND INITIALIZES ALL THE PARAMETERS AND
*     COMMON BLOCK USED BY THE POLYNOMIAL EVALUATION ROUTINES.
*     (CF POLVAL,MTREEF)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*     PARAMETERS:
*
*     LPS: MAXIMUM NUMBER OF TREE SETS;   THE LPS-TH MAY BE RESERVED
*                                         FOR INTERNAL USE
*     LPL: MAXIMUM NUMBER OF POLYNOMIALS
*     LPM: MAXIMUM LENGTH OF TREE NODES+1
*
*-----POLYNOMIAL EVALUATION -------------------------------------------------
      PARAMETER(LPS=5,LPL=16,LPM=10000)
      INTEGER IL(LPM,LPS),IV(LPM,LPS),NNP(LPS),LSW,LSWS(LPS),IPOL
      DOUBLE PRECISION CF(LPM,LPL,LPS)
      COMMON /POLDAT/ CF,IL,IV,NNP,LSW,LSWS,IPOL
*----------------------------------------------------------------------------
*
*     CF:   COEFFICIENTS OF POLYNOMIAL AT EACH TREE NODE
*     IL:   ORDER OF EACH MONOMIAL AT EACH TREE NODE
*     IV:   VARIABLE NUMBER OF THE IL-TH ORDER OF EACH MONOMIAL
*           AT EACH TREE NODE
*     NNP:  TOTAL NUMBER OF TREE NODES
*     LSW:  LAST USED TREE SET NUMBER
*     LSWS: FLAG OF EACH TREE SET TO IDENTIFY IF THE TREE WAS MADE BEFORE
*     IPOL: SWITCH OF POLYNOMIAL EVALUATION ALGORITHM;
*           0: EXPANDED, 1: HORNER
*
      DATA LSW    / 0      /
      DATA LSWS   / 0,0,0,0,0 /
      DATA IPOL   / 0      /
*
      END
*
      SUBROUTINE POLVAL(ILT,IMP,INP,IMA,INA,IMR,INR)
*     **********************************************
*
*     THIS SUBROUTINE EVALUATES POLYNOMIALS, AND WORKS AS CONCATENATOR.
*
*     IPOL: SWITCH ALGORITHM; 0: EXPANDED, 1: HORNER
*
*     LTR:  SWITCH SET NUMBER OF POLYNOMIALS
*       1 -   5   (i) : CALCULATES TREE AND STORES IT IN i-TH ENTRY,
*                       AND EVALUATES POLYNOMIALS
*     101 - 105 (10i) : CALCULATES TREE AND STORES IT IN i-TH ENTRY,
*                       NO POLYNOMIAL EVALUATION
*      -1 -  -5  (-i) : EVALUATES POLYNOMIALS WITH PREVIOUSLY COMPUTED
*                       TREE IN i-TH ENTRY
*         0           : EVALUATES POLYNOMIALS WITH THE LAST ACCESSED TREE
*                       WHICH WAS PREVIOUSLY COMPUTED
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----POLYNOMIAL EVALUATION -------------------------------------------------
      PARAMETER(LPS=5,LPL=16,LPM=10000)
      INTEGER IL(LPM,LPS),IV(LPM,LPS),NNP(LPS),LSW,LSWS(LPS),IPOL
      DOUBLE PRECISION CF(LPM,LPL,LPS)
      COMMON /POLDAT/ CF,IL,IV,NNP,LSW,LSWS,IPOL
*----------------------------------------------------------------------------
*-----NON-ABORTING ARITHMETIC -----------------------------------------------
      INTEGER LARI
      COMMON /ARIST/ LARI
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CFF(LPM,LPL),RMNS(LPL,LNO+1),RMN(LNO+1)
      INTEGER MNS(LPL*(LNO+1)),MN(LNO+1),IARG(LNO+1),NVTMP(LNV)
*
*     CONSISTENCY CHECKS AND PREPARATION
*     **********************************
*
      IF(NTYP(ILT).NE.NRE) CALL FOXNTY(ILT)
      IF(NTYP(INP).NE.NRE) CALL FOXNTY(INP)
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INR).NE.NRE) CALL FOXNTY(INR)
      LTR  = NINT(CC(NBEG(ILT)))
      NPOL = NINT(CC(NBEG(INP)))
      NARG = NINT(CC(NBEG(INA)))
      NRES = NINT(CC(NBEG(INR)))
*
      IF(NPOL.NE.NRES) THEN
         PRINT*,'$$$ ERROR IN POLVAL, '//
     *          'NPOL(3RD ARGUMENT) # NRES(7TH ARGUMENT)'
         CALL FOXDEB
      ELSEIF(NPOL.GT.LPL) THEN
         PRINT*,'$$$ ERROR IN POLVAL, '//
     *          'NPOL(3RD ARGUMENT) > LPL. INCREASE LPL'
         CALL FOXDEB
      ELSEIF(NARG.GT.LNV) THEN
         PRINT*,'$$$ ERROR IN POLVAL, NARG(5TH ARGUMENT) > LNV=',LNV
         CALL FOXDEB
      ENDIF
*
      NTA = NTYP(IMA+1)
      DO 10 I=2,NARG
 10   IF(NTYP(IMA+I).NE.NTA) CALL FOXNTY(IMA+I)
*
      LENA = 2
      IF(NTA.EQ.NDA) THEN
         LENA = NMMAX
C     ELSEIF(NTA.EQ.NTM) THEN
C        LENA = LENTM
      ELSEIF(NTA.EQ.NCD) THEN
         LENA = 2*NMMAX
      ELSEIF(NTA.EQ.NVE) THEN
         LENA = NEND(IMA+1)-NBEG(IMA+1)+1
      ENDIF
*
      NTP = NTYP(IMP+1)
      DO 11 I=2,NPOL
 11   IF(NTYP(IMP+I).NE.NTP) CALL FOXNTY(IMP+I)
*
C     IF(NTP.EQ.NTM) THEN
C        DO 12 I=1,NPOL
C        NVTMP(I) = NC(NEND(IMP+I)-1)
C        IF(NVTMP(I).LT.0.AND.LARI.EQ.0) THEN
C           PRINT*,'$$$ ERROR IN POLVAL, TM ARITHMETIC FAILURE'
C           CALL FOXDEB
C        ENDIF
C12      CONTINUE
*
C        CALL TMCHKP(IMP,NPOL,IMA,NTA,NARG)
C        CALL FOXAAL(IBND,NPOL,2)
C        DO 13 I=1,NPOL
C13      CALL TMRBND(IMP+I,IBND+I)
C     ENDIF
*
      IF(LTR.GE.1) THEN
         ISW = MOD(LTR,100)
         CALL MTREEF(ISW,IMP,NPOL)
         LSW = ISW
         LSWS(LSW) = 1
C        WRITE(6,*) ' IN POLVAL: NNP(ISW) = ',NNP(ISW),
C    *              '  , ISW = ',ISW
C        WRITE(6,*) '  I      IL(I,ISW)   IV(I,ISW)       CF(I,1,ISW)'
C        DO 14 I=1,NNP(ISW)+1
C        WRITE(6,'(I3,2I12,G28.16)') I,IL(I,ISW),IV(I,ISW),CF(I,1,ISW)
C14      CONTINUE
         IF(LTR.GT.100) THEN
C           IF(NTP.EQ.NTM) CALL FOXADA(IBND,NPOL)
            RETURN
         ENDIF
      ELSEIF(LTR.EQ.0) THEN
         ISW = LSW
         IF(LSW.EQ.0) THEN
            PRINT*,'$$$ ERROR IN POLVAL, POLVAL IS NOT CALLED BEFORE.'
            CALL FOXDEB
         ENDIF
      ELSE
         ISW = -LTR
         IF(ISW.LE.0.OR.ISW.GT.LPS) THEN
            PRINT*,'$$$ ERROR IN POLVAL, SWITCH ISW IS OUT OF RANGE.'
            PRINT*,'    ISW, BOUND = ',ISW,',',LPS
            CALL FOXDEB
         ELSEIF(LSWS(ISW).EQ.0) THEN
            PRINT*,'$$$ ERROR IN POLVAL, POLVAL IS NOT CALLED BEFORE'
     *             //' WITH L = ',ISW
            CALL FOXDEB
         ENDIF
         LSW = ISW
      ENDIF
*
      DO 20 K=1,NPOL
      DO 20 I=1,NNP(ISW)
      IF(IL(I,ISW).LE.NOCUT.AND.IV(I,ISW).LE.NARG) THEN
         CFF(I,K) = CF(I,K,ISW)
      ELSE
         CFF(I,K) = 0.D0
      ENDIF
 20   CONTINUE
*
*     REAL NUMBERS
*     ************
*
      IF(NTA.EQ.NRE) THEN
         DO 100 K=1,NPOL
 100     RMNS(K,1) = CFF(1,K)
*
         IF(IPOL.EQ.0) THEN
            RMN(1) = 1.D0
            DO 120 I=2,NNP(ISW)
            IVI  = IV(I,ISW)
            ILI  = IL(I,ISW)
            ILIP = ILI+1
            IF(ILI.GT.NOCUT.OR.IVI.GT.NARG) THEN
               RMN(ILIP) = 0.D0
            ELSE
               RMN(ILIP) = RMN(ILI)*CC(NBEG(IMA+IVI))
               DO 110 K=1,NPOL
 110           RMNS(K,1) = RMNS(K,1)+RMN(ILIP)*CFF(I,K)
            ENDIF
 120        CONTINUE
         ELSE
            DO 170 I=2,NNP(ISW)
            IVI  = IV(I,ISW)
            ILI  = IL(I,ISW)
            ILIN = IL(I+1,ISW)
            IF(ILI.LT.ILIN) THEN
               IARG(ILI) = IVI
               DO 130 K=1,NPOL
 130           RMNS(K,ILIN) = CFF(I,K)
            ELSE
               IF(ILI.LE.NOCUT.AND.IVI.LE.NARG) THEN
                  DO 140 K=1,NPOL
                  RMNS(K,ILI) = RMNS(K,ILI)
     *                         +CFF(I,K)*CC(NBEG(IMA+IVI))
 140              CONTINUE
               ENDIF
               IF(ILI.GT.ILIN) THEN
                  DO 160 J=MIN(NOCUT,ILI-1),ILIN,-1
                  IF(IARG(J).LE.NARG) THEN
                     DO 150 K=1,NPOL
 150                 RMNS(K,J) = RMNS(K,J)
     *                          +RMNS(K,J+1)*CC(NBEG(IMA+IARG(J)))
                  ENDIF
 160              CONTINUE
               ENDIF
            ENDIF
 170        CONTINUE
         ENDIF
*
         DO 180 K=1,NPOL
         NTYP(IMR+K) = NRE
 180     CC(NBEG(IMR+K)) = RMNS(K,1)
*
C        IF(NTP.EQ.NTM) THEN
C           CALL TMVERP(IBND,IMR,NPOL,NTA,NVTMP,LENA)
C           CALL FOXADA(IBND,NPOL)
C        ENDIF
*
         RETURN
      ENDIF
*
*     ALL THE OTHER VARIABLE TYPES
*     ****************************
*
*     ALLOCATION OF LOCAL VARIABLES
*     *****************************
*
      CALL FOXALL(MN,NOMAX+1,LENA)
      CALL FOXALL(MNS,NPOL*(NOMAX+1),LENA)
      CALL FOXALL(IST0,1,LENA)
      CALL FOXALL(IST1,1,LENA)
      CALL FOXALL(IST,1,LENA)
      CALL FOXALL(ISR,1,1)
*
      NTYP(ISR) = NRE
      NBSR = NBEG(ISR)
*
*     INITIALIZATION
*     **************
*
      CC(NBSR) = 0.D0
*
      IF(NTA.EQ.NVE) THEN
         CALL VEMRE(IMA+1,ISR,IST0)
         DO 210 K=1,NPOL
         CC(NBSR) = CFF(1,K)
 210     CALL VEARE(IST0,ISR,MNS(K))
      ELSEIF(NTA.EQ.NDA) THEN
         CALL DAZERO(IMA+1,IST0)
         DO 220 K=1,NPOL
         CC(NBSR) = CFF(1,K)
 220     CALL DAARE(IST0,ISR,MNS(K))
C     ELSEIF(NTA.EQ.NTM) THEN
C        CALL TMZERO(IMA+1,IST0)
C        DO 230 K=1,NPOL
C        CC(NBSR) = CFF(1,K)
C230     CALL TMARE(IST0,ISR,MNS(K))
      ELSEIF(NTA.EQ.NCD) THEN
         CALL CDMRE(IMA+1,ISR,IST0)
         DO 240 K=1,NPOL
         CC(NBSR) = CFF(1,K)
 240     CALL CDARE(IST0,ISR,MNS(K))
      ELSEIF(NTA.EQ.NCM) THEN
         CALL CMMRE(IMA+1,ISR,IST0)
         DO 270 K=1,NPOL
         CC(NBSR) = CFF(1,K)
 270     CALL CMARE(IST0,ISR,MNS(K))
      ELSE
         PRINT*,'$$$ ERROR IN POLVAL, DATA TYPE OF 4TH ARGUMENT'//
     *          ' IS NOT SUPPORTED.'
         CALL FOXDEB
      ENDIF
*
*     IPOL = 0 :   INITIALIZATION AND POLYNOMIAL EVALUATION
*     *****************************************************
*
      IF(IPOL.EQ.0) THEN
*     ******************
         CC(NBSR) = 1.D0
*
         IF(NTA.EQ.NVE) THEN
            CALL VEARE(IST0,ISR,MN(1))
            DO 311 I=2,NNP(ISW)
            IVI  = IV(I,ISW)
            ILI  = IL(I,ISW)
            ILIP = ILI+1
            IF(ILI.GT.NOCUT.OR.IVI.GT.NARG) THEN
               CALL VECOP(IST0,MN(ILIP))
            ELSE
               CALL VEMVE(MN(ILI),IMA+IVI,MN(ILIP))
               DO 310 K=1,NPOL
               CFI = CFF(I,K)
               IF(CFI.EQ.0.D0) GOTO 310
               CC(NBSR) = CFI
               CALL VEMRE(MN(ILIP),ISR,IST)
               CALL VEAVE(MNS(K),IST,MNS(K))
 310           CONTINUE
            ENDIF
 311        CONTINUE
         ELSEIF(NTA.EQ.NDA) THEN
            CALL DAARE(IST0,ISR,MN(1))
            DO 321 I=2,NNP(ISW)
            IVI  = IV(I,ISW)
            ILI  = IL(I,ISW)
            ILIP = ILI+1
            IF(ILI.GT.NOCUT.OR.IVI.GT.NARG) THEN
               CALL DACOP(IST0,MN(ILIP))
            ELSE
               CALL DAMDA(MN(ILI),IMA+IVI,MN(ILIP))
               DO 320 K=1,NPOL
               CFI = CFF(I,K)
               IF(CFI.EQ.0.D0) GOTO 320
               CC(NBSR) = CFI
               CALL DAMRE(MN(ILIP),ISR,IST)
               CALL DAADA(MNS(K),IST,IST1)
               CALL DACOP(IST1,MNS(K))
 320           CONTINUE
            ENDIF
 321        CONTINUE
C        ELSEIF(NTA.EQ.NTM) THEN
C           CALL TMARE(IST0,ISR,MN(1))
C           DO 331 I=2,NNP(ISW)
C           IVI  = IV(I,ISW)
C           ILI  = IL(I,ISW)
C           ILIP = ILI+1
C           IF(ILI.GT.NOCUT.OR.IVI.GT.NARG) THEN
C              CALL TMCOP(IST0,MN(ILIP))
C           ELSE
C              CALL TMMTM(MN(ILI),IMA+IVI,MN(ILIP))
C              DO 330 K=1,NPOL
C              CFI = CFF(I,K)
C              IF(CFI.EQ.0.D0) GOTO 330
C              CC(NBSR) = CFI
C              CALL TMMRE(MN(ILIP),ISR,IST)
C              CALL TMATM(MNS(K),IST,IST1)
C              CALL TMCOP(IST1,MNS(K))
C330           CONTINUE
C           ENDIF
C331        CONTINUE
         ELSEIF(NTA.EQ.NCD) THEN
            CALL CDARE(IST0,ISR,MN(1))
            DO 341 I=2,NNP(ISW)
            IVI  = IV(I,ISW)
            ILI  = IL(I,ISW)
            ILIP = ILI+1
            IF(ILI.GT.NOCUT.OR.IVI.GT.NARG) THEN
               CALL CDCOP(IST0,MN(ILIP))
            ELSE
               CALL CDMCD(MN(ILI),IMA+IVI,MN(ILIP))
               DO 340 K=1,NPOL
               CFI = CFF(I,K)
               IF(CFI.EQ.0.D0) GOTO 340
               CC(NBSR) = CFI
               CALL CDMRE(MN(ILIP),ISR,IST)
               CALL CDACD(MNS(K),IST,IST1)
               CALL CDCOP(IST1,MNS(K))
 340           CONTINUE
            ENDIF
 341        CONTINUE
         ELSEIF(NTA.EQ.NCM) THEN
            CALL CMARE(IST0,ISR,MN(1))
            DO 371 I=2,NNP(ISW)
            IVI  = IV(I,ISW)
            ILI  = IL(I,ISW)
            ILIP = ILI+1
            IF(ILI.GT.NOCUT.OR.IVI.GT.NARG) THEN
               CALL CMCOP(IST0,MN(ILIP))
            ELSE
               CALL CMMCM(MN(ILI),IMA+IVI,MN(ILIP))
               DO 370 K=1,NPOL
               CFI = CFF(I,K)
               IF(CFI.EQ.0.D0) GOTO 370
               CC(NBSR) = CFI
               CALL CMMRE(MN(ILIP),ISR,IST)
               CALL CMACM(MNS(K),IST,MNS(K))
 370           CONTINUE
            ENDIF
 371        CONTINUE
         ENDIF
*
      ELSE
*     ****
*
*     IPOL = 1 :   POLYNOMIAL EVALUATION
*     **********************************
*
         DO 800 I=2,NNP(ISW)
*        ===================
*
*        I : NODE NUMBER OF TREE FOR A MONOMIAL
*
         IVI  = IV(I,ISW)
         ILI  = IL(I,ISW)
         ILIN = IL(I+1,ISW)
*
         IF(ILI.LT.ILIN) THEN
*        ====================
*        OPEN A PARENTHESIS FOR THE NEXT ORDER (ONE ORDER UP)
*
            IARG(ILI) = IVI
            JILIN = (ILIN-1)*NPOL
*
            IF(NTA.EQ.NVE) THEN
               DO 410 K=1,NPOL
               CC(NBSR) = CFF(I,K)
               JK = JILIN+K
 410           CALL VEARE(IST0,ISR,MNS(JK))
            ELSEIF(NTA.EQ.NDA) THEN
               DO 420 K=1,NPOL
               CC(NBSR) = CFF(I,K)
               JK = JILIN+K
 420           CALL DAARE(IST0,ISR,MNS(JK))
C           ELSEIF(NTA.EQ.NTM) THEN
C              DO 430 K=1,NPOL
C              CC(NBSR) = CFF(I,K)
C              JK = JILIN+K
C430           CALL TMARE(IST0,ISR,MNS(JK))
            ELSEIF(NTA.EQ.NCD) THEN
               DO 440 K=1,NPOL
               CC(NBSR) = CFF(I,K)
               JK = JILIN+K
 440           CALL CDARE(IST0,ISR,MNS(JK))
            ELSEIF(NTA.EQ.NCM) THEN
               DO 470 K=1,NPOL
               CC(NBSR) = CFF(I,K)
               JK = JILIN+K
 470           CALL CMARE(IST0,ISR,MNS(JK))
            ENDIF
*
         ELSE
*        ====
*        ADD THE MONOMIAL CONTRIBUTION IN THE PARENTHESIS
*        (AND CLOSE THE PARENTHESES IF NECESSARY)
*
            IF(ILI.LE.NOCUT.AND.IVI.LE.NARG) THEN
*
               JILI = (ILI-1)*NPOL
*
               IF(NTA.EQ.NVE) THEN
                  DO 510 K=1,NPOL
                  CFI = CFF(I,K)
                  IF(CFI.EQ.0.D0) GOTO 510
                  CC(NBSR) = CFI
                  JK = JILI+K
                  CALL VEMRE(IMA+IVI,ISR,IST)
                  CALL VEAVE(MNS(JK),IST,MNS(JK))
 510              CONTINUE
               ELSEIF(NTA.EQ.NDA) THEN
                  DO 520 K=1,NPOL
                  CFI = CFF(I,K)
                  IF(CFI.EQ.0.D0) GOTO 520
                  CC(NBSR) = CFI
                  JK = JILI+K
                  CALL DAMRE(IMA+IVI,ISR,IST)
                  CALL DAADA(MNS(JK),IST,IST1)
                  CALL DACOP(IST1,MNS(JK))
 520              CONTINUE
C              ELSEIF(NTA.EQ.NTM) THEN
C                 DO 530 K=1,NPOL
C                 CFI = CFF(I,K)
C                 IF(CFI.EQ.0.D0) GOTO 530
C                 CC(NBSR) = CFI
C                 JK = JILI+K
C                 CALL TMMRE(IMA+IVI,ISR,IST)
C                 CALL TMATM(MNS(JK),IST,IST1)
C                 CALL TMCOP(IST1,MNS(JK))
C530              CONTINUE
               ELSEIF(NTA.EQ.NCD) THEN
                  DO 540 K=1,NPOL
                  CFI = CFF(I,K)
                  IF(CFI.EQ.0.D0) GOTO 540
                  CC(NBSR) = CFI
                  JK = JILI+K
                  CALL CDMRE(IMA+IVI,ISR,IST)
                  CALL CDACD(MNS(JK),IST,IST1)
                  CALL CDCOP(IST1,MNS(JK))
 540              CONTINUE
               ELSEIF(NTA.EQ.NCM) THEN
                  DO 570 K=1,NPOL
                  CFI = CFF(I,K)
                  IF(CFI.EQ.0.D0) GOTO 570
                  CC(NBSR) = CFI
                  JK = JILI+K
                  CALL CMMRE(IMA+IVI,ISR,IST)
                  CALL CMACM(MNS(JK),IST,MNS(JK))
 570              CONTINUE
               ENDIF
*
            ENDIF
*
            IF(ILI.GT.ILIN) THEN
*           --------------------
*           CLOSE THE PARENTHESES DOWN TO THE ORDER OF THE NEXT MONOMIAL
*
               DO 700 J=MIN(NOCUT,ILI-1),ILIN,-1
*
               IF(IARG(J).LE.NARG) THEN
*
                  JJP = J*NPOL
                  JJ  = (J-1)*NPOL
*
                  IF(NTA.EQ.NVE) THEN
                     DO 610 K=1,NPOL
                     JPK = JJP+K
                     JK  = JJ +K
                     CALL VEMVE(MNS(JPK),IMA+IARG(J),IST)
 610                 CALL VEAVE(MNS(JK),IST,MNS(JK))
                  ELSEIF(NTA.EQ.NDA) THEN
                     DO 620 K=1,NPOL
                     JPK = JJP+K
                     JK  = JJ +K
                     CALL DAMDA(MNS(JPK),IMA+IARG(J),IST)
                     CALL DAADA(MNS(JK),IST,IST1)
 620                 CALL DACOP(IST1,MNS(JK))
C                 ELSEIF(NTA.EQ.NTM) THEN
C                    DO 630 K=1,NPOL
C                    JPK = JJP+K
C                    JK  = JJ +K
C                    CALL TMMTM(MNS(JPK),IMA+IARG(J),IST)
C                    CALL TMATM(MNS(JK),IST,IST1)
C630                 CALL TMCOP(IST1,MNS(JK))
                  ELSEIF(NTA.EQ.NCD) THEN
                     DO 640 K=1,NPOL
                     JPK = JJP+K
                     JK  = JJ +K
                     CALL CDMCD(MNS(JPK),IMA+IARG(J),IST)
                     CALL CDACD(MNS(JK),IST,IST1)
 640                 CALL CDCOP(IST1,MNS(JK))
                  ELSEIF(NTA.EQ.NCM) THEN
                     DO 670 K=1,NPOL
                     JPK = JJP+K
                     JK  = JJ +K
                     CALL CMMCM(MNS(JPK),IMA+IARG(J),IST)
 670                 CALL CMACM(MNS(JK),IST,MNS(JK))
                  ENDIF
*
               ENDIF
*
 700           CONTINUE
*
            ENDIF
*           -----
*
         ENDIF
*        =====
*
 800     CONTINUE
*===     ========
*
      ENDIF
*     *****
*
*     SET THE RETURNING RESULTS
*     *************************
*
      IF(NTA.EQ.NVE) THEN
         DO 910 K=1,NPOL
 910     CALL VECOP(MNS(K),IMR+K)
      ELSEIF(NTA.EQ.NDA) THEN
         DO 920 K=1,NPOL
 920     CALL DACOP(MNS(K),IMR+K)
C     ELSEIF(NTA.EQ.NTM) THEN
C        DO 930 K=1,NPOL
C930     CALL TMCOP(MNS(K),IMR+K)
      ELSEIF(NTA.EQ.NCD) THEN
         DO 940 K=1,NPOL
 940     CALL CDCOP(MNS(K),IMR+K)
      ELSEIF(NTA.EQ.NCM) THEN
         DO 970 K=1,NPOL
 970     CALL CMCOP(MNS(K),IMR+K)
      ENDIF
*
*     DEALLOCATION OF LOCAL VARIABLES
*     *******************************
*
      CALL FOXDAL(ISR,1)
      CALL FOXDAL(IST,1)
      CALL FOXDAL(IST1,1)
      CALL FOXDAL(IST0,1)
      CALL FOXDAL(MNS,NPOL*(NOMAX+1))
      CALL FOXDAL(MN,NOMAX+1)
*
C     IF(NTP.EQ.NTM) THEN
C        CALL TMVERP(IBND,IMR,NPOL,NTA,NVTMP,LENA)
C        CALL FOXADA(IBND,NPOL)
C     ENDIF
*
      RETURN
      END
*
      SUBROUTINE POLSET(IIPOL)
*     ************************
*
*     THIS SUBROUTINE SETS THE ALGORITHM NUMBER TO EVALUATE POLYNOMIALS
*     IN THE SUBROUTINE POLVAL
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----POLYNOMIAL EVALUATION -------------------------------------------------
      PARAMETER(LPS=5,LPL=16,LPM=10000)
      INTEGER IL(LPM,LPS),IV(LPM,LPS),NNP(LPS),LSW,LSWS(LPS),IPOL
      DOUBLE PRECISION CF(LPM,LPL,LPS)
      COMMON /POLDAT/ CF,IL,IV,NNP,LSW,LSWS,IPOL
*----------------------------------------------------------------------------
*
      IF(NTYP(IIPOL).NE.NRE) CALL FOXNTY(IIPOL)
      IPOL = NINT(CC(NBEG(IIPOL)))
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     REAL OPERATIONS                                                         *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE RECOP(INA,INC)
*     *************************
*
*     THIS SUBROUTINE STORES THE REAL VARIABLE INA IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = CC(NBEG(INA))
      RETURN
      END
*
      SUBROUTINE LRE(INA,INC)
*     ***********************
*
*     THIS SUBROUTINE RETURNS THE ALLOCATION LENGTH OF REAL IN INC.
*     THE INCOMING ARGUMENT INA IS VOID.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = 1.D0
*
      RETURN
      END
*
      SUBROUTINE REZERO(INA,INC)
*     **************************
*
*     THIS SUBROUTINE SETS INC TO REAL ZERO
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = 0.D0
      RETURN
      END
*
      SUBROUTINE RESQR(INA,INC)
*     *************************
*
*     THIS SUBROUTINE SQUARES THE REAL VARIABLE INA IN INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = CC(NBEG(INA))*CC(NBEG(INA))
      RETURN
      END
*
      SUBROUTINE REISRT(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE RECIPROCAL OF THE SQUARE ROOT OF REAL INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = 1.D0/SQRT(CC(NBEG(INA)))
      RETURN
      END
*
      SUBROUTINE REISR3(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE RECIPROCAL TO THE POWER 3/2 OF REAL INA.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = 1.D0/(CC(NBEG(INA))*SQRT(CC(NBEG(INA))))
      RETURN
      END
*
      SUBROUTINE REERF(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE ERROR FUNCTION ERF OF REAL INA
*
*     USE 7.1.26 PAGE 299 HANDBOOK MATH. FUNCTIONS.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      A1  = .254829592D0
      A2  = -.284496736D0
      A3  = 1.421413741D0
      A4  = -1.453152027D0
      A5  = 1.061405429D0
      P   = .3275911D0
      PI  = 4.D0*DATAN(1.D0)
      SPI = SQRT(PI)/2.D0
*
      A0  = CC(NBEG(INA))
      AA = A0
      IF(A0.LT.0.D0) AA = -A0
*
      IF(AA.LT.1.D0) THEN
        XN   = AA
        FAK  = 1.D0
        L    = 1
        X1   = L*XN/FAK
        DO 341 KI = 1,20
          L   = L *(-1)
          XN  = XN*AA*AA
          FAK = FAK*KI
          D   = 2*KI+1
          X1  = X1+L*XN/FAK/D
  341   CONTINUE
        E2  = X1/SPI
      ELSE
        E1 = DEXP(-AA*AA)
        T  = 1.D0/(1+P*AA)
        E2 = 1-T*(A1+T*(A2+T*(A3+T*(A4+T*A5))))*E1
      ENDIF
*
      IF(A0.LT.0.D0) E2 = -E2
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = E2
      PRINT*,'$$$ WARNING IN REERF, THIS ROUTINE IS UNDER DEVELOPMENT.'
      PRINT*,'    CONTACT US AT BERZ@MSU.EDU'
      RETURN
      END
*
      SUBROUTINE REPRI(INA,IUNIT)
*     ***************************
*
*     THIS SUBROUTINE PRINTS REAL INA TO UNIT IUNIT
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      WRITE(IUNIT,'(1X,G23.16E3)') CC(NBEG(INA))
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     STRING OPERATIONS                                                       *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE STCOP(INA,INB)
*     *************************
*
*     THIS SUBROUTINE COPIES STRING INA INTO INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      LA = NEND(INA)-IA
      LB = NMAX(INB)-IB
*
      LM = MIN(LA,LB)
*
      DO 10 I=0,LM
 10   NC(IB+I) = NC(IA+I)
*
      NEND(INB) = IB+LM
      NTYP(INB) = NST
*
      RETURN
      END
*
      SUBROUTINE LST(INA,INC)
*     ***********************
*
*     THIS SUBROUTINE RETURNS THE ALLOCATION LENGTH OF STRING IN INC.
*     THE REAL VARIABLE INA SUPPLIES THE LENGTH OF STRING.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      A0 = CC(NBEG(INA))
      IF(A0.LE.0.D0) THEN
         PRINT*,'$$$ WARNING IN LST, ARGUMENT MUST BE POSITIVE INTEGER.'
         CC(NBEG(INC)) = 1.D0
      ELSEIF(A0.LE.2147483647.D0) THEN
         NA0 = NINT(A0)
         IF(ABS(A0-NA0).GT.1.D-12)
     *      PRINT*,'$$$ WARNING IN LST, ARGUMENT MUST BE INTEGER.'
         CC(NBEG(INC)) = DBLE(NA0)
      ELSE
         PRINT*,'$$$ WARNING IN LST, ARGUMENT IS TOO LARGE.'
         CC(NBEG(INC)) = A0
      ENDIF
*
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE STUST(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MERGES TWO STRINGS INA AND INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(INB.EQ.INC) CALL DANFI(INA,INB,INC)
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NST) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)
      LA = NEND(INA)-IA
      LB = NEND(INB)-IB
      LC = NMAX(INC)-IC
*
      LM = MIN(LA+LB+1,LC)
      MI = MIN(LA,LM)
      IB = IB-MI-1
*
      DO 10 I=0,MI
 10   NC(IC+I) = NC(IA+I)
*
      DO 20 I=MI+1,LM
 20   NC(IC+I) = NC(IB+I)
*
      NEND(INC) = IC+LM
      NTYP(INC) = NST
*
      RETURN
      END
*
      SUBROUTINE SUBSTR(INA,INB,INC,IND)
*     **********************************
*
*     THIS SUBROUTINE PICKS THE SUBSTRING OF INA FROM INB TO INC AND
*     STORES IT IN IND
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      IF(NTYP(INC).NE.NRE) CALL FOXNTY(INC)
*
      IA = NBEG(INA)
      ID = NBEG(IND)
      LA = NEND(INA)-IA+1
      NNB = NINT(CC(NBEG(INB)))
      NNC = NINT(CC(NBEG(INC)))
      N1 = MAX(1,MIN(NNB,NNC))
      N2 = MIN(LA,MAX(NNB,NNC))
      IDL = ID+N2-N1
*
      IF(IDL.GT.NMAX(IND)) CALL FOXERV('SUBSTR')
*
      DO 10 I=IA+N1-1,IA+N2-1
      NC(ID) = NC(I)
 10   ID = ID+1
*
      NEND(IND) = IDL
      NTYP(IND) = NST
*
      RETURN
      END
*
      SUBROUTINE STSUB(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PICKS THE SUBSTRING OF INA SPECIFIED BY INB AND
*     STORES IT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE.AND.NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      LA = NEND(INA)-IA+1
      NNB1 = NINT(CC(NBEG(INB)))
      NNB2 = NNB1
      IF(NTYP(INB).EQ.NVE) NNB2 = NINT(CC(NBEG(INB)+1))
      N1 = MAX(1,MIN(NNB1,NNB2))
      N2 = MIN(LA,MAX(NNB1,NNB2))
      ICL = IC+N2-N1
*
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('STSUB')
*
      DO 10 I=IA+N1-1,IA+N2-1
      NC(IC) = NC(I)
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NST
*
      RETURN
      END
*
      SUBROUTINE STLST(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A STRING "<" OPERATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF((INA.EQ.INC).OR.(INB.EQ.INC)) CALL DANFI(INA,INB,INC)
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NST) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      LA = NEND(INA)-IA
      LB = NEND(INB)-IB
*
      LM = MIN(LA,LB)
*
      NTYP(INC) = NLO
      NC(NBEG(INC)) = 1
*
      DO 10 I=0,LM
      IF(NC(IA+I).GT.NC(IB+I)) THEN
         NC(NBEG(INC)) = 0
         RETURN
      ELSEIF(NC(IA+I).LT.NC(IB+I)) THEN
         RETURN
      ENDIF
 10   CONTINUE
*
      IF(LA.GE.LB) THEN
         NC(NBEG(INC)) = 0
         RETURN
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE STGST(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A STRING ">" OPERATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF((INA.EQ.INC).OR.(INB.EQ.INC)) CALL DANFI(INA,INB,INC)
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NST) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      LA = NEND(INA)-IA
      LB = NEND(INB)-IB
*
      LM = MIN(LA,LB)
*
      NTYP(INC) = NLO
      NC(NBEG(INC)) = 1
*
      DO 10 I=0,LM
      IF(NC(IA+I).GT.NC(IB+I)) THEN
         RETURN
      ELSEIF(NC(IA+I).LT.NC(IB+I)) THEN
         NC(NBEG(INC)) = 0
         RETURN
      ENDIF
 10   CONTINUE
*
      IF(LA.LE.LB) THEN
         NC(NBEG(INC)) = 0
         RETURN
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE STEST(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A STRING "=" OPERATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF((INA.EQ.INC).OR.(INB.EQ.INC)) CALL DANFI(INA,INB,INC)
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NST) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      LA = NEND(INA)-IA
      LB = NEND(INB)-IB
*
      NTYP(INC) = NLO
      NC(NBEG(INC)) = 1
*
      IF(LA.NE.LB) THEN
         NC(NBEG(INC)) = 0
         RETURN
      ENDIF
*
      DO 10 I=0,LA
      IF(NC(IA+I).NE.NC(IB+I)) THEN
         NC(NBEG(INC)) = 0
         RETURN
      ENDIF
 10   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE STNST(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A STRING "#" OPERATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF((INA.EQ.INC).OR.(INB.EQ.INC)) CALL DANFI(INA,INB,INC)
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NST) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      LA = NEND(INA)-IA
      LB = NEND(INB)-IB
*
      NTYP(INC) = NLO
      NC(NBEG(INC)) = 0
*
      IF(LA.NE.LB) THEN
         NC(NBEG(INC)) = 1
         RETURN
      ENDIF
*
      DO 10 I=0,LA
      IF(NC(IA+I).NE.NC(IB+I)) THEN
         NC(NBEG(INC)) = 1
         RETURN
      ENDIF
 10   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE STPRI(INA,IUNIT)
*     ***************************
*     THIS SUBROUTINE PRINTS A STRING TO UNIT IUNIT
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
*
      WRITE(IUNIT,'(4096A)') (CHAR(NC(J)),J=NBEG(INA),NEND(INA))
*
      RETURN
      END
*
      SUBROUTINE LOCST(INA,INC)
*     *************************
*
*     THIS SUBROUTINE WRITES THE LOGICAL INA TO THE STRING INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER TARGET*5
*
      IF(NTYP(INA).NE.NLO) CALL FOXNTY(INA)
*
      IF(NC(NBEG(INA)).EQ.0) THEN
         TARGET = 'FALSE'
      ELSE
         TARGET = 'TRUE '
      ENDIF
*
      LCH = 5
      IC = NBEG(INC)
      NEND(INC) = IC+LCH-1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('LOCST')
*
      DO 10 I=1,LCH
      NC(IC) = ICHAR(TARGET(I:I))
 10   IC = IC+1
*
      NTYP(INC) = NST
*
      RETURN
      END
*
      SUBROUTINE RECSTC(INA,INC)
*     **************************
*
*     THIS SUBROUTINE WRITES THE REAL NUMBER OR THE COMPLEX NUMBERS
*     INA WITH THE COSY DEFAULT FORMAT TO THE STRING INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER FORM*10
*
      NTA = NTYP(INA)
      IF(NTA.NE.NRE.AND.NTA.NE.NCM) CALL FOXNTY(INA)
*
      IF(NTA.EQ.NCM) THEN
         LFORM = 9
         FORM = '(G17.9E3)'
      ELSE
         LFORM = 10
         FORM = '(G23.16E3)'
      ENDIF
*
      CALL FOXALL(INB,1,LFORM)
      IB = NBEG(INB)
      NEND(INB) = IB+LFORM-1
      NTYP(INB) = NST
*
      DO 10 I=1,LFORM
      NC(IB) = ICHAR(FORM(I:I))
 10   IB = IB+1
*
      CALL RECST(INA,INB,INC)
      CALL FOXDAL(INB,1)
*
      RETURN
      END
*
      SUBROUTINE RECST(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE WRITES THE REAL NUMBER OR THE COMPLEX NUMBERS
*     INA WITH THE FORMAT INB TO THE STRING INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER CB,FORM*20,TARGET*200
*
      NTA = NTYP(INA)
*
      IF(NTA.NE.NRE.AND.NTA.NE.NCM) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NST) CALL FOXNTY(INB)
*
      FORM = '                    '
      TARGET = '                                                  '//
     *         '                                                  '//
     *         '                                                  '//
     *         '                                                  '
*
      IT = 0
      DO 10 I=NBEG(INB),NEND(INB)
      CB = CHAR(NC(I))
      IF(CB.NE.' ') THEN
         IT = IT+1
         FORM(IT:IT) = CB
      ENDIF
  10  CONTINUE
*
      IF(IT.LT.4.OR.IT.GT.10) THEN
         GOTO 40
      ELSEIF(FORM(1:1).NE.'('.OR.FORM(IT:IT).NE.')') THEN
         PRINT*,'$$$ ERROR IN RECST, FORMAT IS NOT ENCLOSED BY ( )'
         GOTO 50
      ELSEIF(INDEX('IiFfEeDdGg',FORM(2:2)).EQ.0.OR.
     *       INDEX('0123456789',FORM(3:3)).EQ.0) THEN
         GOTO 40
      ENDIF
*
      DO 20 I=4,IT-1
      IF(INDEX('Ee0123456789.',FORM(I:I)).EQ.0) GOTO 40
  20  CONTINUE
*
      IF(NTA.EQ.NRE) THEN
         IFORM = INDEX('Ii',FORM(2:2))
         CALL REWSTO(CC(NBEG(INA)),FORM,IFORM,TARGET)
      ELSEIF(NTA.EQ.NCM) THEN
         CALL CMCSTO(CC(NBEG(INA)),CC(NBEG(INA)+1),FORM,TARGET)
      ENDIF
*
      LC = MIN(ILAST(TARGET,1,200),NMAX(INC)-NBEG(INC)+1)
*
      DO 30 I=1,LC
  30  NC(NBEG(INC)+I-1) = ICHAR(TARGET(I:I))
*
      NTYP(INC) = NST
      NEND(INC) = NBEG(INC)+LC-1
*
      RETURN
*
  40  PRINT*,'$$$ ERROR IN RECST, GIVEN FORMAT IS NOT SUPPORTED'
  50  PRINT*,'    FORMAT = '//FORM
      CALL FOXDEB
*
      END
*
      SUBROUTINE REWSTO(A,FORM,IFORM,C)
*     *********************************
*
*     THIS SUBROUTINE WRITES THE REAL A TO THE STRING C WITH THE FORMAT FORM
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER FORM*20,C*100,FORMD*22,CDGT*2
*
      IF(IFORM.EQ.0) THEN
         WRITE(C,FORM,ERR=20) A
         RETURN
      ELSE
         IL = ILAST(FORM,1,20)
         LPOI = IL
         IF(INDEX(FORM(1:IL),'.').NE.0) LPOI = INDEX(FORM(1:IL),'.')
         IF(LPOI.EQ.5) THEN
            READ(FORM(3:LPOI-1),'(I2)',ERR=20) IDGT
         ELSEIF(LPOI.EQ.4) THEN
            READ(FORM(3:LPOI-1),'(I1)',ERR=20) IDGT
         ELSE
            PRINT*,'$$$ ERROR IN REWSTO, TOO MANY DIGITS ARE DEMANDED'
            GOTO 30
         ENDIF
C        FORTRAN INTEGER LIMIT IS BETWEEN -2147483648 AND 2147483647
         IF(ABS(A).LE.2147483647.D0) THEN
            WRITE(C,FORM,ERR=10) NINT(A)
            RETURN
         ENDIF
 10      FORMD = '                      '
         IF(ABS(A).LT.1.D15.AND.IDGT.GT.INT(LOG10(ABS(A)))+2) THEN
            WRITE(CDGT,'(I2)') 0
            FORMD = '(F'//FORM(3:LPOI-1)//'.'//CDGT//')'
         ELSE
            WRITE(CDGT,'(I2)') IDGT-6
            FORMD = '(D'//FORM(3:LPOI-1)//'.'//CDGT//')'
         ENDIF
         WRITE(C,FORMD,ERR=20) A
         RETURN
      ENDIF
*
 20   PRINT*,'$$$ ERROR IN REWSTO, '//
     *       'SYSTEM CANNOT ACCEPT THE GIVEN FORMAT'
 30   PRINT*,'    FORMAT = '//FORM
      CALL FOXDEB
*
      END
*
      SUBROUTINE STCRE(INA,INB)
*     *************************
*
*     THIS SUBROUTINE CONVERTS THE STRING INA TO THE REAL NUMBER INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CHARACTER C*80
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
*
      IF(NEND(INA)-NBEG(INA).GT.79) THEN
         CC(NBEG(INB)) = 0.D0
         NTYP(INB) = NRE
         RETURN
      ENDIF
*
      J = 0
      DO 10 I=NBEG(INA),NEND(INA)
      J = J+1
  10  C(J:J) = CHAR(NC(I))
*
      CALL VALCH(C,1,J,CC(NBEG(INB)),IER)
      NTYP(INB) = NRE
*
      RETURN
      END
*
      SUBROUTINE VALCH(A,I1,I2,VAL,IER)
*     *********************************
*
*     THIS FUNCTION EXTRACTS THE VALUE OF A NUMBER STORED IN A CHARACTER
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER A*(*),NUM*10,FORM*10,DIGITS*3
      DATA NUM / '1234567890' /
*
      IEND = ILAST(A,I1,I2)
      IDGT = 0
      IFRA = 0
      LPOI = 0
      LEXP = 0
      I = I1-1
 10   I = I+1
      IF(I.GT.IEND) GOTO 20
      IF(A(I:I).EQ.' ') GOTO 10
      IBEG = I
      IDGT = IEND-IBEG+1
      I = I-1
*
 15   I = I+1
      IF(I.GT.IEND) GOTO 20
      IF(INDEX(NUM,A(I:I)).NE.0) THEN
         IF(LEXP.EQ.0.AND.LPOI.EQ.1) IFRA = IFRA+1
      ELSEIF(A(I:I).EQ.'.') THEN
         LPOI = 1
      ELSEIF(INDEX('deDE',A(I:I)).NE.0) THEN
         LEXP = 1
      ELSEIF(A(I:I).EQ.' ') THEN
         GOTO 100
      ENDIF
      GOTO 15
*
 20   CONTINUE
      VAL  = 0.D0
      IER = 0
      IF(IDGT.GT.0) THEN
         IF(LEXP.EQ.1) THEN
            FORM(1:2) = '(E'
         ELSE
            FORM(1:2) = '(F'
         ENDIF
         DIGITS = '   '
         WRITE(DIGITS,'(I3)') IDGT
         FORM(3:6) = DIGITS//'.'
         DIGITS = '   '
         WRITE(DIGITS,'(I3)') IFRA
         FORM(7:10) = DIGITS//')'
         READ(A(IBEG:IEND),FORM,ERR=100) VAL
      ENDIF
      RETURN
*
 100  CONTINUE
      VAL = 0.D0
      IER = 1
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     LOGICAL OPERATIONS                                                      *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE LOCOP(INA,INB)
*     *************************
*
*     THIS SUBROUTINE COPIES LOGICAL INA INTO INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NLO) CALL FOXNTY(INA)
*
      NC(NBEG(INB)) = NC(NBEG(INA))
      NTYP(INB) = NLO
*
      RETURN
      END
*
      SUBROUTINE LLO(INA,INC)
*     ***********************
*
*     THIS SUBROUTINE RETURNS THE ALLOCATION LENGTH OF LOGICAL IN INC.
*     THE INCOMING ARGUMENT INA IS VOID.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = 1.D0
*
      RETURN
      END
*
      SUBROUTINE LOALO(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A LOGICAL OR OPERATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NLO) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NLO) CALL FOXNTY(INB)
*
      LC = NC(NBEG(INA))+NC(NBEG(INB))
      IF(LC.EQ.2) LC = 1
      NC(NBEG(INC)) = LC
      NTYP(INC) = NLO
*
      RETURN
      END
*
      SUBROUTINE LOMLO(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A LOGICAL AND OPERATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NLO) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NLO) CALL FOXNTY(INB)
*
      NC(NBEG(INC)) = NC(NBEG(INA))*NC(NBEG(INB))
      NTYP(INC) = NLO
*
      RETURN
      END
*
      SUBROUTINE RELRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A LOGICAL "<" OPERATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      NTYP(INC) = NLO
      IF(CC(NBEG(INA)).LT.CC(NBEG(INB))) THEN
         NC(NBEG(INC)) = 1
      ELSE
         NC(NBEG(INC)) = 0
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE REGRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A LOGICAL ">" OPERATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      NTYP(INC) = NLO
      IF(CC(NBEG(INA)).GT.CC(NBEG(INB))) THEN
         NC(NBEG(INC)) = 1
      ELSE
         NC(NBEG(INC)) = 0
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE REERE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A LOGICAL "=" OPERATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      NTYP(INC) = NLO
      IF(CC(NBEG(INA)).EQ.CC(NBEG(INB))) THEN
         NC(NBEG(INC)) = 1
      ELSE
         NC(NBEG(INC)) = 0
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE RENRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A LOGICAL "#" OPERATION
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      NTYP(INC) = NLO
      IF(CC(NBEG(INA)).NE.CC(NBEG(INB))) THEN
         NC(NBEG(INC)) = 1
      ELSE
         NC(NBEG(INC)) = 0
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE LTRUE(INC)
*     *********************
*
*     THIS SUBROUTINE SETS INC TO LOGICAL TRUE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NC(NBEG(INC)) = 1
      NTYP(INC) = NLO
*
      RETURN
      END
*
      SUBROUTINE LFALSE(INC)
*     **********************
*
*     THIS SUBROUTINE SETS INC TO LOGICAL FALSE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NC(NBEG(INC)) = 0
      NTYP(INC) = NLO
*
      RETURN
      END
*
      SUBROUTINE LO(INA,INC)
*     **********************
*
*     THIS SUBROUTINE SETS INC TO LOGICAL TRUE  IF INA IS 1,
*                          AND TO LOGICAL FALSE IF INA IS 0
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      A0 = CC(NBEG(INA))
*
      IF(ABS(A0).LT.1.D-12) THEN
         NC(NBEG(INC)) = 0
      ELSEIF(ABS(A0-1.D0).LT.1.D-12) THEN
         NC(NBEG(INC)) = 1
      ELSE
         PRINT*,'$$$ ERROR IN LO, '//
     *              'ARGUMENT HAS TO BE 1 (TRUE) OR 0 (FALSE).'
         CALL FOXDEB
      ENDIF
*
      NTYP(INC) = NLO
*
      RETURN
      END
*
      SUBROUTINE LONOT(INA,INC)
*     *************************
*
*     THIS NEGATES THE LOGICAL INA AND STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NLO) CALL FOXNTY(INA)
      NC(NBEG(INC)) = 1-NC(NBEG(INA))
      NTYP(INC) = NLO
*
      RETURN
      END
*
      SUBROUTINE LOPRI(INA,IUNIT)
*     ***************************
*
*     THIS SUBROUTINE PRINTS INA TO UNIT IUNIT
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NLO) CALL FOXNTY(INA)
*
      IF(NC(NBEG(INA)).EQ.0) THEN
         WRITE(IUNIT,'(A)') ' FALSE '
      ELSE
         WRITE(IUNIT,'(A)') ' TRUE  '
      ENDIF
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     VECTOR OPERATIONS                                                       *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE VECOP(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COPIES VECTOR INA INTO VECTOR INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VECOP')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE LVE(INA,INC)
*     ***********************
*
*     THIS SUBROUTINE RETURNS THE ALLOCATION LENGTH OF VECTOR IN INC.
*     THE REAL VARIABLE INA SUPPLIES THE NUMBER OF COMPONENTS OF VECTOR.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      A0 = CC(NBEG(INA))
      IF(A0.LE.0.D0) THEN
         PRINT*,'$$$ WARNING IN LVE, ARGUMENT MUST BE POSITIVE INTEGER.'
         CC(NBEG(INC)) = 2.D0
      ELSEIF(A0.LE.2147483647.D0) THEN
         NA0 = NINT(A0)
         IF(ABS(A0-NA0).GT.1.D-12)
     *      PRINT*,'$$$ WARNING IN LVE, ARGUMENT MUST BE INTEGER.'
         CC(NBEG(INC)) = DBLE(MAX(2,NA0))
      ELSE
         PRINT*,'$$$ WARNING IN LVE, ARGUMENT IS TOO LARGE.'
         CC(NBEG(INC)) = A0
      ENDIF
*
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE REURE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE UNITES THE TWO REALS INA AND INB INTO THE VECTOR INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(INB.EQ.INC) CALL DANFI(INA,INB,INC)
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IC = NBEG(INC)
      ICL = IC+1
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('REURE')
*
      CC(IC)   = CC(NBEG(INA))
      CC(IC+1) = CC(NBEG(INB))
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE REUVE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MERGES THE REAL INA WITH THE VECTOR INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(INB.EQ.INC) CALL DANFI(INA,INB,INC)
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IB = NBEG(INB)
      IC = NBEG(INC)
      ICL = IC+NEND(INB)-IB+1
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('REUVE')
*
      CC(IC) = CC(NBEG(INA))
      IC = IC+1
      DO 10 I=IB,NEND(INB)
      CC(IC) = CC(I)
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEURE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MERGES THE VECTOR INA WITH THE REAL INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(INB.EQ.INC) CALL DANFI(INA,INB,INC)
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA+1
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEURE')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)
 10   IC = IC+1
      CC(IC) = CC(NBEG(INB))
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEUVE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MERGES TWO VECTORS INA AND INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(INB.EQ.INC) CALL DANFI(INA,INB,INC)
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA+NEND(INB)-IB+1
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEUVE')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)
 10   IC = IC+1
*
      DO 20 I=IB,NEND(INB)
      CC(IC) = CC(I)
 20   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEAVE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS TWO VECTORS INA AND INB AND STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)
      LA = NEND(INA)-IA+1
      LB = NEND(INB)-IB+1
      ICL = IC+LA-1
*
      IF(ICL.GT.NMAX(INC).OR.LA.NE.LB) THEN
         IF(LA.NE.LB) THEN
            PRINT*,'$$$ ERROR IN VEAVE, LENGTH NOT MATCHED: ',INA,INB
            CALL FOXDEB
         ELSE
            CALL FOXERV('VEAVE')
         ENDIF
      ENDIF
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)+CC(IB)
      IB = IB+1
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VESVE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE SUBTRACTS TWO VECTORS INA AND INB AND STORES THE
*     RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)
      LA = NEND(INA)-IA+1
      LB = NEND(INB)-IB+1
      ICL = IC+LA-1
*
      IF(ICL.GT.NMAX(INC).OR.LA.NE.LB) THEN
         IF(LA.NE.LB) THEN
            PRINT*,'$$$ ERROR IN VESVE, LENGTH NOT MATCHED: ',INA,INB
            CALL FOXDEB
         ELSE
            CALL FOXERV('VESVE')
         ENDIF
      ENDIF
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)-CC(IB)
      IB = IB+1
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEMVE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MULTIPLIES TWO VECTORS INA AND INB AND STORES THE
*     RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)
      LA = NEND(INA)-IA+1
      LB = NEND(INB)-IB+1
      ICL = IC+LA-1
*
      IF(ICL.GT.NMAX(INC).OR.LA.NE.LB) THEN
         IF(LA.NE.LB) THEN
            PRINT*,'$$$ ERROR IN VEMVE, LENGTH NOT MATCHED: ',INA,INB
            CALL FOXDEB
         ELSE
            CALL FOXERV('VEMVE')
         ENDIF
      ENDIF
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)*CC(IB)
      IB = IB+1
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEDVE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE DIVIDES TWO VECTORS INA AND INB AND STORES THE
*     RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)
      LA = NEND(INA)-IA+1
      LB = NEND(INB)-IB+1
      ICL = IC+LA-1
*
      IF(ICL.GT.NMAX(INC).OR.LA.NE.LB) THEN
         IF(LA.NE.LB) THEN
            PRINT*,'$$$ ERROR IN VEDVE, LENGTH NOT MATCHED: ',INA,INB
            CALL FOXDEB
         ELSE
            CALL FOXERV('VEDVE')
         ENDIF
      ENDIF
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)/CC(IB)
      IB = IB+1
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEARE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS THE VECTOR INA AND THE SCALAR INB AND STORES
*     THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      CB = CC(NBEG(INB))
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEARE')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)+CB
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE REAVE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS THE REAL INA AND THE VECTOR INB AND STORES
*     THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IB = NBEG(INB)
      IC = NBEG(INC)
      ICL = IC+NEND(INB)-IB
      CA = CC(NBEG(INA))
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('REAVE')
*
      DO 10 I=IB,NEND(INB)
      CC(IC) = CA+CC(I)
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VESRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE SUBTRACTS THE VECTOR INA AND THE SCALAR INB AND STORES
*     THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      CB = CC(NBEG(INB))
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VESRE')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)-CB
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE RESVE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS THE REAL INA AND THE VECTOR INB AND STORES
*     THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IB = NBEG(INB)
      IC = NBEG(INC)
      ICL = IC+NEND(INB)-IB
      CA = CC(NBEG(INA))
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('RESVE')
*
      DO 10 I=IB,NEND(INB)
      CC(IC) = CA-CC(I)
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEMRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MULTIPLIES THE VECTOR INA AND THE SCALAR INB AND STORES
*     THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      CB = CC(NBEG(INB))
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEMRE')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)*CB
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE REMVE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MULTIPLIES THE REAL INA AND THE VECTOR INB AND STORES
*     THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IB = NBEG(INB)
      IC = NBEG(INC)
      ICL = IC+NEND(INB)-IB
      CA = CC(NBEG(INA))
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('REMVE')
*
      DO 10 I=IB,NEND(INB)
      CC(IC) = CA*CC(I)
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEDRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE DIVIDES THE VECTOR INA AND THE SCALAR INB AND STORES
*     THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      CB = CC(NBEG(INB))
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEDRE')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)/CB
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE REDVE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE DIVIDES THE REAL INA AND THE VECTOR INB AND STORES
*     THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IB = NBEG(INB)
      IC = NBEG(INC)
      ICL = IC+NEND(INB)-IB
      CA = CC(NBEG(INA))
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('REDVE')
*
      DO 10 I=IB,NEND(INB)
      CC(IC) = CA/CC(I)
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VESQR(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE SQUARE OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VESQR')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = CC(I)*CC(I)
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VESQRT(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE SQUARE ROOT OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VESQRT')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = SQRT(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEISRT(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE INVERSE SQUARE ROOT OF THE VECTOR INA AND
*     STORES THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEISRT')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = 1.D0/SQRT(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEISR3(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE RECIPROCAL TO THE POWER 3/2
*     OF THE VECTOR INA AND STORES THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEISR3')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = 1.D0/(CC(I)*SQRT(CC(I)))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEPOW(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE EXPONENTIATES THE VECTOR INA WITH RE INB
*     AND STORES THE RESULT TO INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      CB = CC(NBEG(INB))
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEPOW')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = (CC(I))**CB
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEEXP(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE EXPONENTIAL OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEEXP')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = EXP(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VELOG(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE LOGARITHM OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VELOG')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = LOG(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VESINE(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE SINE OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VESINE')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = SIN(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VECOSE(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE COSINE OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VECOSE')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = COS(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEDTAN(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE TANGENT OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEDTAN')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = TAN(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEASIN(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE ARC SINE OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEASIN')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = ASIN(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEACOS(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE ARC COSINE OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEACOS')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = ACOS(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEATAN(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE ARC TANGENT OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEATAN')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = ATAN(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VESINH(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE HYPERBOLIC SINE OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VESINH')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = SINH(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VECOSH(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE HYPERBOLIC COSINE OF THE VECTOR INA AND
*     STORES THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VECOSH')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = COSH(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VETANH(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE HYPERBOLIC TANGENT OF THE VECTOR INA AND
*     STORES THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VETANH')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = TANH(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEINT(INA,INC)
*     *************************
*
*     THIS SUBROUTINE TRUNCATES EACH ELEMENT OF THE VECTOR INA TO AN INTEGER
*     AND STORES THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEINT')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = DINT(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VENINT(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE NEAREST INTEGER FOR EACH ELEMENT OF THE
*     VECTOR INA AND STORES THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      ICL = IC+NEND(INA)-IA
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VENINT')
*
      DO 10 I=IA,NEND(INA)
      CC(IC) = DNINT(CC(I))
 10   IC = IC+1
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE VEPRI(INA,IUNIT)
*     ***************************
*
*     THIS SUBROUTINE PRINTS A VECTOR TO UNIT IUNIT
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IAE = NEND(INA)
*
      WRITE(IUNIT,'(A,250(5(G14.7E3,A):/2X))')
     * '  ',((CC(I),' ',I=5*J+IA,MIN(5*J+IA+4,IAE-1)),J=0,(IAE-IA)/5),
     *       (CC(I),I=MAX(IA+1,IAE),IAE),'  '
*
      RETURN
      END
*
      SUBROUTINE VELSET(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE SETS THE ELEMENT WITH NUMBER IN INB IN THE VECTOR
*     INA TO THE VALUE IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INC).NE.NRE) CALL FOXNTY(INC)
*
      IF(NTYP(INA).EQ.NVE) THEN
         LA = NEND(INA)-NBEG(INA)+1
      ELSEIF(NTYP(INA).EQ.NRE) THEN
         LA = 1
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      NB = NINT(CC(NBEG(INB)))
      IF(NB.LT.1.OR.NB.GT.LA) THEN
         PRINT*,'$$$ ERROR IN VELSET, INDEX IS OUT OF BOUND'
         PRINT*,'    INDEX, BOUND = ',NB,',',LA
         CALL FOXDEB
      ENDIF
*
      CC(NBEG(INA)+NB-1) = CC(NBEG(INC))
*
      RETURN
      END
*
      SUBROUTINE VELGET(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE PICKS THE ELEMENT WITH NUMBER IN INB OUT OF THE VECTOR
*     OR REAL IN INA AND STORES IT IN REAL INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).EQ.NVE) THEN
         LA = NEND(INA)-NBEG(INA)+1
      ELSEIF(NTYP(INA).EQ.NRE) THEN
         LA = 1
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      NB = NINT(CC(NBEG(INB)))
      IF(NB.LT.1.OR.NB.GT.LA) THEN
         PRINT*,'$$$ ERROR IN VELGET, INDEX IS OUT OF BOUND'
         PRINT*,'    INDEX, BOUND = ',NB,',',LA
         CALL FOXDEB
      ENDIF
*
      CC(NBEG(INC)) = CC(NBEG(INA)+NB-1)
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE VESUB(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PICKS THE SUB VECTOR OF INA SPECIFIED BY INB AND
*     STORES IT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IA = NBEG(INA)
      IF(NTYP(INA).EQ.NVE) THEN
         LA = NEND(INA)-IA+1
      ELSEIF(NTYP(INA).EQ.NRE) THEN
         LA = 1
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      NNB1 = NINT(CC(NBEG(INB)))
      NNB2 = NINT(CC(NBEG(INB)+1))
      N1 = MIN(NNB1,NNB2)
      N2 = MAX(NNB1,NNB2)
      IF(N1.LT.1.OR.N2.GT.LA) THEN
         PRINT*,'$$$ ERROR IN VESUB, INDEXES ARE OUT OF BOUND'
         PRINT*,'    INDEXES, BOUND = ',N1,',',N2,',',LA
         CALL FOXDEB
      ENDIF
*
      IC = NBEG(INC)-1
      DO 10 I=IA+N1-1,IA+N2-1
      IC = IC+1
      CC(IC) = CC(I)
 10   CONTINUE
*
      IF(IC.GT.NMAX(INC)) CALL FOXERV('VESUB')
      NEND(INC) = IC
      NTYP(INC) = NVE
      IF(N1.EQ.N2) NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE VERE(INA,INC)
*     ************************
*
*     THIS SUBROUTINE COMPUTES THE AVERAGE OF THE VECTOR INA AND
*     STORES IT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      LA = NEND(INA)-NBEG(INA)+1
*
*     WHEN LA IS TOO LARGE, IT WILL HELP TO COMPUTE AS   SUM = SUM+CC(I)/LA
*
      SUM = 0.D0
      DO 10 I=NBEG(INA),NEND(INA)
 10   SUM = SUM+CC(I)
*
      CC(NBEG(INC)) = SUM/LA
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE VECNST(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE MAXIMUM OF ABSOLUTE VALUE OF THE VECTOR
*     INA AND STORES THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      SUM = 0.D0
      DO 10 I=NBEG(INA),NEND(INA)
  10  SUM = MAX(SUM,ABS(CC(I)))
*
      CC(NBEG(INC)) = SUM
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE VEABS(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE ABSOLUTE VALUE OF THE VECTOR INA AND STORES
*     THE RESULT IN INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      SUM = 0.D0
      DO 10 I=NBEG(INA),NEND(INA)
  10  SUM = SUM+ABS(CC(I))
*
      CC(NBEG(INC)) = SUM
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE VEZERO(IMA,INNA,INB)
*     *******************************
*
*     THIS SUBROUTINE CHECKS IF ANY OF THE INNA COMPONENTS IN AN VE ARRAY IMA
*     EXCEEDS INB, IN WHICH CASE ALL RESPECTIVE COMPONENTS IN IMA ARE SET
*     TO ZERO.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INNA).NE.NRE) CALL FOXNTY(INNA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      IF(NTYP(IMA).GT.0) THEN
         PRINT*,'$$$ ERROR IN VEZERO, FIRST ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ENDIF
*
      NA = NINT(CC(NBEG(INNA)))
      CB = CC(NBEG(INB))
*
      DO 100 I=1,NA
      IF(NTYP(IMA+I).NE.NRE.AND.NTYP(IMA+I).NE.NVE) CALL FOXNTY(IMA+I)
      IA = NBEG(IMA+I)
      DO 100 J=0,NEND(IMA+I)-IA
      IF(ABS(CC(IA+J)).GT.CB) THEN
         DO 10 K=1,NA
         IP = NBEG(IMA+K)+J
         IF(IP.LE.NEND(IMA+K)) CC(IP) = 0.D0
 10      CONTINUE
      ENDIF
 100  CONTINUE
*
      RETURN
      END
*
      SUBROUTINE VELMIN(INA,INC)
*     **************************
*
*     STORES THE MINIMUM ELEMENT OF A VECTOR INA IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
*
      A = CC(IA)
      DO 10 I=IA+1,NEND(INA)
 10   A = MIN(A,CC(I))
*
      CC(NBEG(INC)) = A
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE VELMAX(INA,INC)
*     **************************
*
*     STORES THE MAXIMUM ELEMENT OF A VECTOR INA IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
*
      A = CC(IA)
      DO 10 I=IA+1,NEND(INA)
 10   A = MAX(A,CC(I))
*
      CC(NBEG(INC)) = A
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE VESET(INA,INB)
*     *************************
*
*     THIS SUBROUTINE FILLS THE VECTOR (INA) WITH THE REAL NUMBER (INB)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      CB = CC(NBEG(INB))
*
      DO 10 I=IA,NEND(INA)
 10   CC(I) = CB
*
      RETURN
      END
*
      SUBROUTINE VEDOT(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE COMPUTES THE DOT PRODUCT OF TWO VECTORS INA AND INB
*     AND STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NVE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      LA = NEND(INA)-IA+1
      LB = NEND(INB)-IB+1
*
      DOT = 0.D0
*
      IF(LA.LE.LB) THEN
         DO 10 I=IA,NEND(INA)
         DOT = DOT+CC(I)*CC(IB)
 10      IB = IB+1
      ELSE
         DO 20 I=IB,NEND(INB)
         DOT = DOT+CC(IA)*CC(I)
 20      IA = IA+1
      ENDIF
*
      CC(NBEG(INC)) = DOT
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE VEUNIT(INA,INC)
*     **************************
*
*     THIS SUBROUTINE NORMALIZES THE VECTOR INA AND STORES THE RESULTING
*     UNIT VECTOR IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      LA = NEND(INA)-IA+1
      IC = NBEG(INC)
      ICL = IC+LA-1
      IF(ICL.GT.NMAX(INC)) CALL FOXERV('VEUNIT')
*
      AMAX = 0.D0
      DO 10 I=IA,NEND(INA)
 10   AMAX = MAX(AMAX,ABS(CC(I)))
*
      IF(AMAX.GT.0.D0) THEN
         SUM = 0.D0
         DO 20 I=IC,ICL
         CC(I) = CC(IA)/AMAX
         SUM = SUM+CC(I)*CC(I)
 20      IA = IA+1
*
         SUM = SQRT(SUM)
         DO 30 I=IC,ICL
 30      CC(I) = CC(I)/SUM
      ELSE
         DO 40 I=IC,ICL
 40      CC(I) = 0.D0
      ENDIF
*
      NEND(INC) = ICL
      NTYP(INC) = NVE
*
      RETURN
      END
*%%
      SUBROUTINE VEFIL(IMA,A1,A2,NV,ND,NRN)
*     *************************************
*
*     THIS SUBROUTINE FILLS THE VECTOR ARRAY IMA FOR SCANNING.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION A1(NV),A2(NV)
      INTEGER NBMA(LNV)
*
      IF(NTYP(IMA).GT.0) THEN
         PRINT*,'$$$ ERROR IN VEFIL, FIRST ARGUMENT IS NOT ARRAY.'
         CALL FOXDEB
      ELSEIF(ND.LT.0.OR.NRN.LT.0) THEN
         PRINT*,'$$$ ERROR IN VEFIL, ND OR NRN IS NEGATIVE.'
         PRINT*,'    ND, NRN = ',ND,',',NRN
         CALL FOXDEB
      ELSEIF(ND+NRN.EQ.0) THEN
         RETURN
      ENDIF
*
      NTP = ND**NV+NRN
      LEN = NMAX(IMA+1)-NBEG(IMA+1)+1
      IF(LEN.LT.NTP) CALL FOXERV('VEFIL')
*
      DO 10 I=1,NV
 10   NBMA(I) = NBEG(IMA+I)
*
      LAST = ND**NV
      IF(ND.GT.1) THEN
         DO 30 I=1,NV
         KSTP = ND**I
         LEND = ND**(I-1)
            DO 30 J=1,ND
            XJ = A1(I)+DBLE(J-1)*(A2(I)-A1(I))/DBLE(ND-1)
               DO 30 K=1,LAST,KSTP
               IS = NBMA(I)+LEND*(J-1)+K-1
                  DO 20 L=1,LEND
                  CC(IS) = XJ
 20               IS = IS+1
 30      CONTINUE
      ELSEIF(ND.EQ.1) THEN
         DO 40 I=1,NV
 40      CC(NBMA(I)) = A1(I)
      ENDIF
*
      DO 50 J=1,NRN
      DO 50 I=1,NV
      IS = NBMA(I)+LAST+J-1
      CC(IS) = A1(I)+BRAN(0)*(A2(I)-A1(I))
 50   CONTINUE
*
      IF(NTP.EQ.1) THEN
         DO 60 I=1,NV
 60      NTYP(IMA+I) = NRE
      ELSE
         DO 70 I=1,NV
         NTYP(IMA+I) = NVE
 70      NEND(IMA+I) = NBMA(I)+NTP-1
      ENDIF
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     GRAPHICS OPERATIONS                                                     *
*                                                                             *
*******************************************************************************
*
      BLOCKDATA GRDT
*     **************
*
*     THIS BLOCKDATA SETS THE BASE FILE NAME FBASE AND OUTPUT FILE COUNTER ICNT
*     FOR GRAPHICS FILE OUTPUT
*
*-----GRAPH FILE OUTPUT -----------------------------------------------------
      INTEGER ICNT
      CHARACTER FBASE*20
      COMMON /GRFILE/ ICNT,FBASE
*----------------------------------------------------------------------------
*
      DATA FBASE / 'pic                 ' /
      DATA ICNT / 1 /
*
      END
*
      SUBROUTINE CLEAR(INA)
*     *********************
*
*     THIS SUBROUTINE CLEARS THE GRAPHICS OBJECT INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IA = NBEG(INA)
      CC(IA) = 0.D0
      NC(IA) = 0
*
      NTYP(INA) = NGR
      NEND(INA) = IA
      IF(NEND(INA).GT.NMAX(INA)) CALL FOXERV('CLEAR')
*
      RETURN
      END
*
      SUBROUTINE GRCOP(INA,INB)
*     *************************
*
*     THIS SUBROUTINE COPIES GRAPHICS OBJECT INA INTO INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NGR) CALL FOXNTY(INA)
*
      IB = NBEG(INB)-1
*
      DO 10 I=NBEG(INA),NEND(INA)
      IB = IB+1
      CC(IB) = CC(I)
 10   NC(IB) = NC(I)
*
      NTYP(INB) = NGR
      NEND(INB) = IB
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('GRCOP')
*
      RETURN
      END
*
      SUBROUTINE LGR(INA,INC)
*     ***********************
*
*     THIS SUBROUTINE RETURNS THE APPROXIMATE ALLOCATION LENGTH OF
*     GRAPHICS OBJECT IN INC. THE REAL VARIABLE INA SUPPLIES THE NUMBER
*     OF GR ELEMENTS OF GRAPHICS OBJECT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
*
      A0 = CC(NBEG(INA))
      IF(A0.LE.0.D0) THEN
         PRINT*,'$$$ WARNING IN LGR, ARGUMENT MUST BE POSITIVE INTEGER.'
         CC(NBEG(INC)) = 9.D0
      ELSEIF(A0.LE.2147483647.D0) THEN
         NA0 = NINT(A0)
         IF(ABS(A0-NA0).GT.1.D-12)
     *      PRINT*,'$$$ WARNING IN LGR, ARGUMENT MUST BE INTEGER.'
         CC(NBEG(INC)) = 200.D0+4.D0*DBLE(NA0)
      ELSE
         PRINT*,'$$$ WARNING IN LGR, ARGUMENT IS TOO LARGE.'
         CC(NBEG(INC)) = 200.D0+4.D0*A0
      ENDIF
*
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE GRUGR(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MERGES TWO GRAPHICS OBJECTS INA AND INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NGR) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NGR) CALL FOXNTY(INB)
      IF(INB.EQ.INC) CALL DANFI(INA,INB,INC)
*
      ICS = NBEG(INC)
      IC = ICS-1
*
      DO 10 I=NBEG(INA),NEND(INA)
      IC = IC+1
      CC(IC) = CC(I)
 10   NC(IC) = NC(I)
*
      DO 20 I=NBEG(INB)+1,NEND(INB)
      IC = IC+1
      CC(IC) = CC(I)
 20   NC(IC) = NC(I)
*
      IF(NC(NBEG(INB)).GT.0) NC(ICS) = NEND(INA)-NBEG(INA)+NC(NBEG(INB))
*
      NTYP(INC) = NGR
      NEND(INC) = IC
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('GRUGR')
*
      RETURN
      END
*
      SUBROUTINE GRMIMA(INA,INB,INC,IND,INE,INF,ING)
*     **********************************************
*
*     THIS SUBROUTINE RETURNS THE MINIMUM, MAXIMUM COORDINATES OF A GRAPHICS
*     OBJECT INA, IN (INB,INC) FOR X, (IND,INE) FOR Y, (INF,ING) FOR Z.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NGR) CALL FOXNTY(INA)
*
      CALL GRPRE(INA,NELE,NARI,1)
*
      NTYP(INB) = NRE
      NTYP(INC) = NRE
      NTYP(IND) = NRE
      NTYP(INE) = NRE
      NTYP(INF) = NRE
      NTYP(ING) = NRE
*
      CC(NBEG(INB)) = BDG(1,1)
      CC(NBEG(INC)) = BDG(1,2)
      CC(NBEG(IND)) = BDG(2,1)
      CC(NBEG(INE)) = BDG(2,2)
      CC(NBEG(INF)) = BDG(3,1)
      CC(NBEG(ING)) = BDG(3,2)
*
      NEND(INB) = NBEG(INB)
      NEND(INC) = NBEG(INC)
      NEND(IND) = NBEG(IND)
      NEND(INE) = NBEG(INE)
      NEND(INF) = NBEG(INF)
      NEND(ING) = NBEG(ING)
*
      RETURN
      END
*
      SUBROUTINE GROUTF(INA,INB)
*     **************************
*
*     THIS SUBROUTINE SETS THE OUTPUT FILE NAME AND FILE COUNTER
*     FOR GRAPHICS FILE OUTPUT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GRAPH FILE OUTPUT -----------------------------------------------------
      INTEGER ICNT
      CHARACTER FBASE*20
      COMMON /GRFILE/ ICNT,FBASE
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      LA = MIN(20,NEND(INA)-IA+1)
      IF(LA.GT.0) THEN
         FBASE = ''
         DO 10 I=1,LA
         FBASE(I:I) = CHAR(NC(IA+I-1))
 10      CONTINUE
      ENDIF
      ICNT = ABS(NINT(CC(NBEG(INB))))
*
      RETURN
      END
*
      SUBROUTINE GRMOVE(INA,INB,INC,IND)
*     **********************************
*
*     THIS SUBROUTINE APPENDS A MOVE TO THE GRAPHICS OBJECT IND
*     INA, INB, INC ARE COORDINATES
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      IF(NTYP(INC).NE.NRE) CALL FOXNTY(INC)
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      NC(NBEG(IND)) = NEND(IND)-NBEG(IND)+1
*
      ID = NEND(IND)
*
      CC(ID+1) = CC(NBEG(INA))
      CC(ID+2) = CC(NBEG(INB))
      CC(ID+3) = CC(NBEG(INC))
*
      NC(ID+1) = 100000003
      NC(ID+2) = 0
      NC(ID+3) = 0
*
      NEND(IND) = ID+3
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRMOVE')
*
      RETURN
      END
*
      SUBROUTINE GRDRAW(INA,INB,INC,IND)
*     **********************************
*
*     THIS SUBROUTINE APPENDS A DRAW TO THE GRAPHICS OBJECT IND
*     INA, INB, INC ARE COORDINATES
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      IF(NTYP(INC).NE.NRE) CALL FOXNTY(INC)
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      NC(NBEG(IND)) = NEND(IND)-NBEG(IND)+1
*
      ID = NEND(IND)
*
      CC(ID+1) = CC(NBEG(INA))
      CC(ID+2) = CC(NBEG(INB))
      CC(ID+3) = CC(NBEG(INC))
*
      NC(ID+1) = 200000003
      NC(ID+2) = 0
      NC(ID+3) = 0
*
      NEND(IND) = ID+3
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRDRAW')
*
      RETURN
      END
*
      SUBROUTINE GRDOT(INA,INB,INC,IND)
*     *********************************
*
*     THIS SUBROUTINE APPENDS A MOVE AND A DOT TO THE GRAPHICS OBJECT IND
*     INA, INB, INC ARE COORDINATES
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      IF(NTYP(INC).NE.NRE) CALL FOXNTY(INC)
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      NC(NBEG(IND)) = NEND(IND)-NBEG(IND)+1
*
      ID = NEND(IND)
*
      CC(ID+1) = CC(NBEG(INA))
      CC(ID+2) = CC(NBEG(INB))
      CC(ID+3) = CC(NBEG(INC))
*
      NC(ID+1) = 300000003
      NC(ID+2) = 0
      NC(ID+3) = 0
*
      NEND(IND) = ID+3
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRDOT')
*
      RETURN
      END
*
      SUBROUTINE GRTRI(INA,INB,INC,IND)
*     *********************************
*
*     THIS SUBROUTINE APPENDS A TRIANGLE TO THE GRAPHICS OBJECT IND
*     INA, INB, INC ARE COORDINATES OF THE THIRD/NEW POSITION TO FORM
*     THE TRIANGLE TOGETHER WITH THE LAST TWO GIVEN POSITIONS
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      IF(NTYP(INC).NE.NRE) CALL FOXNTY(INC)
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      NC(NBEG(IND)) = NEND(IND)-NBEG(IND)+1
*
      ID = NEND(IND)
*
      CC(ID+1) = CC(NBEG(INA))
      CC(ID+2) = CC(NBEG(INB))
      CC(ID+3) = CC(NBEG(INC))
*
      NC(ID+1) = 400000003
      NC(ID+2) = 0
      NC(ID+3) = 0
*
      NEND(IND) = ID+3
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRTRI')
*
      RETURN
      END
*
      SUBROUTINE GRPOLY(INA,INB,INC,IND)
*     **********************************
*
*     THIS SUBROUTINE APPENDS A DA CURVE OR SURFACE SUPPLIED BY A DA OR TM
*     ARRAY OF POLYNOMIALS TO THE GRAPHICS OBJECT IND.
*     THE CURVE OR SURFACE IS GIVEN BY A DA/TM ARRAY INA. SUPPLY AT LEAST
*     2 SPACIAL DA/TM VECTORS (FOR X,Y) OR 3 DA/TM VECTORS FOR X,Y,Z.
*     THE COLOR IS SPECIFIED BY INB.
*        RE OR VE: SOLID COLOR SET BY GRCOLR. IF -1, PREVIOUSLY SET COLOR.
*        DA ARRAY: AT LEAST 3 DA/TM VECTORS FOR RGB, AND 4TH FOR ALPHA OF RGBA.
*        TM        WHEN SUPPLYING ONLY 3 DA/TM VECTORS FOR RGB, THE ALPHA
*                  POLYNOMIAL IS SET TO BE 1 FOR THE FULL OPACITY.
*     THE CURVE OR SURFACE DA/TM VECTORS ARE PARAMETRIZED BY 2 DA VARIABLES.
*     THE DA VARIABLE NUMBER(S) OF THE PARAMETER(S) IS/ARE GIVEN BY INC.
*        RE: DA VARIABLE NUMBER OF THE PARAMETER (CURVE)
*        VE: DA VARIABLE NUMBERS OF THE 1ST AND THE 2ND PARAMETERS (SURFACE)
*     WHEN PRE-SPECIFYING THE DISCRETIZATION, HAVE INC TO BE AN ARRAY, AND
*        GIVE NP1,NP2 BY RE OR VE IN THE 2ND COMPONENT OF THE ARRAY.
*
*     THE POLYNOMIALS ARE STORED IN IND BY APPENDING THE FOLLOWING.
*     IDB IS THE MEMORY ADDRESS WHERE THE ADDITION STARTS FOR IND.
*     THE 2 HEADER SLOTS CONTAIN THE PROPERTIES OF THE POLYNOMIALS.
*     THERE ARE 7 POLYNOMIALS WITH NX, NY, NZ, NR, NG, NB, NA TERMS.
*     THE NON-ZERO COEFFICIENTS WITH THE EXPONENTS INFO FOLLOW THE HEADER,
*     IN THE ORGERING OF X, Y, Z (3D SPACIAL) AND R, G, B, A (RGBA COLOR).
*
*     ADDRESS I       NC(I)                               CC(I)
*
*     IDB         500000000+2+SUM_N'S            NX 00 NY 00 NZ
*     IDB+1         1000000 (CURVE)        NR 00 NG 00 NB 00 NA
*                OR 2000000 (SURFACE)
*
*     ...           P2 0 P1                COEFFICIENT OF V1**P1*V2**P2
*
*     "P2 0 P1" MEANS P2*100+P1, WHERE P1 AND P2 ARE EXPONENTS OF
*     THE 1ST AND THE 2ND PARAMETER VARIABLES.
*     IF THE POLYNOMIAL IS A TM, THE REMAINDER BOUND IS STORED IN CC(I)
*     WITH NC(I) = -1.
*     "NX 00 NY 00 NZ"        MEANS         NX*1E6+NY*1E3+NZ.
*     "NR 00 NG 00 NB 00 NA"  MEANS  NR*1E9+NG*1E6+NB*1E3+NA.
*     THE ARRAY IPLEN(J) CONTAINS NX, NY, NZ, NR, NG, NB, NA.
*     WHEN THE PRE-SPECIFIED DISCRETIZATION NP1,NP2 ARE GIVEN, THEY ARE STORED
*     IN NC(IDB+1) AS "100 NP1 000" (CURVE) OR "200 NP1 00 NP2" (SURFACE).
*
*     IF ANY TM POLYNOMIAL IS TM ARITHMETIC FAILURE, AND THE ARITHMETIC ERROR
*     HANDLING SWITCH LARI=1 TO CARRY ON THE COSY COMPUTATIONS, THE FOLLOWING
*     IS PROCESSED.
*     X,Y,Z: CC(I=IDB) = -1
*     RGBA:  THE COLOR POLYNOMIALS ARE IGNORED.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*-----NON-ABORTING ARITHMETIC -----------------------------------------------
      INTEGER LARI
      COMMON /ARIST/ LARI
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION RTM(7)
      INTEGER JJ(LNV),IPLEN(7)
      CHARACTER CODE*7
*
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      CODE = 'XYZRGBA'
      DO 10 J=1,7
      RTM(J) = 0.D0
      IPLEN(J) = 0
 10   CONTINUE
*
*     INITIAL PROCESSING OF INB (COLOR)
*
      NCARI = 0
      LB = 0
      IF(NTYP(INB).EQ.NRE) THEN
         IF(CC(NBEG(INB)).GE.0.D0) CALL GRCOLR(INB,IND)
      ELSEIF(NTYP(INB).EQ.NVE) THEN
         CALL GRCOLR(INB,IND)
      ELSEIF(NTYP(INB).NE.-1) THEN
         PRINT*,'$$$ ERROR IN GRPOLY, 2ND ARGUMENT MUST BE RE, VE '//
     *          'OR A DA/TM ARRAY.'
         CALL FOXDEB
      ELSE
         DO 20 J=1,MIN(NDIM(NMAX(INB)),4)
C        IF(NTYP(INB+J).EQ.NTM) THEN
C           IBE = NEND(INB+J)
C           NVTM = NC(IBE-1)
C           IF(NVTM.LT.0) THEN
C              NCARI = NCARI+1
C              IF(LARI.EQ.1) THEN
C                 PRINT*,'$$$ WARNING IN GRPOLY, TM ARITHMETIC FAILURE'
C    *                   //' -- '//CODE(J+3:J+3)//' OF RGBA'
C              ELSEIF(LARI.EQ.0) THEN
C                 PRINT*,'$$$ ERROR IN GRPOLY, TM ARITHMETIC FAILURE'
C    *                   //' -- '//CODE(J+3:J+3)//' OF RGBA'
C                 CALL FOXDEB
C              ENDIF
C           ELSE
C              RTM(J+3) = MAX(ABS(CC(IBE-2)),ABS(CC(IBE-1)))
C           ENDIF
C        ELSEIF(NTYP(INB+J).NE.NDA) THEN
         IF(NTYP(INB+J).NE.NDA) THEN
            GOTO 30
         ENDIF
         LB = LB+1
 20      CONTINUE
 30      IF(LB.LT.3) THEN
            PRINT*,'$$$ ERROR IN GRPOLY, NOT ENOUGH COLOR DA/TM '//
     *             'VARIABLES'
            CALL FOXDEB
         ENDIF
      ENDIF
*
*     INITIAL PROCESSING OF INA (SPACIAL POLYNOMIALS)
*
      IF(NTYP(INA).NE.-1) THEN
         PRINT*,'$$$ ERROR IN GRPOLY, 1ST ARGUMENT MUST BE AN ARRAY.'
         CALL FOXDEB
      ENDIF
*
      NSARI = 0
      LA = 0
      DO 40 J=1,MIN(NDIM(NMAX(INA)),3)
C     IF(NTYP(INA+J).EQ.NTM) THEN
C        IAE = NEND(INA+J)
C        NVTM = NC(IAE-1)
C        IF(NVTM.LT.0) THEN
C           NSARI = NSARI+1
C           IF(LARI.EQ.1) THEN
C              PRINT*,'$$$ WARNING IN GRPOLY, TM ARITHMETIC FAILURE'
C    *                //' -- '//CODE(J:J)//' OF XYZ'
C           ELSEIF(LARI.EQ.0) THEN
C              PRINT*,'$$$ ERROR IN GRPOLY, TM ARITHMETIC FAILURE'
C    *                //' -- '//CODE(J:J)//' OF XYZ'
C              CALL FOXDEB
C           ENDIF
C        ELSE
C           RTM(J) = MAX(ABS(CC(IAE-2)),ABS(CC(IAE-1)))
C        ENDIF
C     ELSEIF(NTYP(INA+J).NE.NDA) THEN
      IF(NTYP(INA+J).NE.NDA) THEN
         GOTO 50
      ENDIF
      LA = LA+1
 40   CONTINUE
 50   IF(LA.LT.2) THEN
         PRINT*,'$$$ ERROR IN GRPOLY, NOT ENOUGH SPACIAL DA/TM '//
     *          'VARIABLES'
         CALL FOXDEB
      ENDIF
*
      IF(NSARI.GT.0) THEN
*        PROCESS THE COLOR PART (INB) AHEAD OF THIS IF-BLOCK
         ID = NEND(IND)+1
         NC(ID) = 500000001
         CC(ID) = -1.D0
         NEND(IND) = ID
         IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRPOLY')
         RETURN
      ENDIF
*
*     PROCESSING OF INC
*
      INCV = INC
      INCN = 0
      IF(NTYP(INC).EQ.-1) THEN
         NDIMC = NDIM(NMAX(INC))
         IF(NDIMC.GE.1) INCV = INC+1
         IF(NDIMC.GE.2) INCN = INC+2
      ENDIF
*
*     DETERMINE IOBJ, IV1, IV2
*
      IC = NBEG(INCV)
      IF(NTYP(INCV).EQ.NRE) THEN
*      * CURVE *
         IOBJ = 1
         IV1 = NINT(CC(IC))
         IV2 = IV1
      ELSEIF(NTYP(INCV).EQ.NVE) THEN
*      * SURFACE *
         IOBJ = 2
         IV1 = NINT(CC(IC  ))
         IV2 = NINT(CC(IC+1))
         IF(IV2.LT.1.OR.IV2.GT.NVMAX) THEN
            PRINT*,'$$$ ERROR IN GRPOLY, 2ND PARAMETER IS OUT OF '//
     *             'DA RANGE'
            PRINT*,'    VARIABLE INDEX, BOUND = ',IV2,',',ABS(NVMAX)
            CALL FOXDEB
         ELSEIF(IV2.EQ.IV1) THEN
            PRINT*,'$$$ ERROR IN GRPOLY, 1ST AND 2ND PARAMETERS '//
     *             'MUST DIFFER'
            PRINT*,'    1ST, 2ND PARAMETER NUMBERS = ',IV1,IV2
            CALL FOXDEB
         ENDIF
      ELSE
         PRINT*,'$$$ ERROR IN GRPOLY, 3RD ARGUMENT MUST BE RE, VE '//
     *          'OR A RE/VE ARRAY.'
         CALL FOXNTY(INCV)
      ENDIF
*
      IF(IV1.LT.1.OR.IV1.GT.NVMAX) THEN
         PRINT*,'$$$ ERROR IN GRPOLY, 1ST PARAMETER IS OUT OF DA RANGE'
         PRINT*,'    VARIABLE INDEX, BOUND = ',IV1,',',ABS(NVMAX)
         CALL FOXDEB
      ENDIF
*
*     DETERMINE PRE-SPECIFIED NP1, NP2
*
      NP1 = 0
      NP2 = 0
*
      IF(INCN.NE.0) THEN
         IC = NBEG(INCN)
         IF(NTYP(INCN).EQ.NRE) THEN
            NP1 = MIN(NGP,MAX(0,NINT(CC(IC))))
         ELSEIF(NTYP(INCN).EQ.NVE) THEN
            NP1 = MIN(NGP,MAX(0,NINT(CC(IC  ))))
            NP2 = MIN(NGP,MAX(0,NINT(CC(IC+1))))
         ELSE
            PRINT*,'$$$ ERROR IN GRPOLY, 3RD ARGUMENT MUST BE RE, VE '//
     *             'OR A RE/VE ARRAY.'
            CALL FOXNTY(INCN)
         ENDIF
*
*        IF ANY POSITIVE INTEGER IS GIVEN, A NON-ZERO PAIR (NP1,NP2) IS STORED.
*        HOWEVER, FOR A CURVE (IOBJ=1), NP2 IS 0.
         IF(IOBJ.EQ.1) THEN
            NP2 = 0
         ELSEIF(IOBJ.EQ.2) THEN
            IF(NP1.GE.1.OR.NP2.GE.1) THEN
               NP1 = MAX(1,NP1)
               NP2 = MAX(1,NP2)
            ENDIF
         ENDIF
      ENDIF
*
*     APPEND THE CURVE OR SURFACE POLYNOMIAL INFORMATION TO IND
*
      NC(NBEG(IND)) = NEND(IND)-NBEG(IND)+1
*
      IS = NEND(IND)+2
      IV2F = IOBJ-1
*
*     SPACIAL POLYNOMIALS
*
      DO 130 J=1,LA
      IAB = NBEG(INA+J)
      IAE = NEND(INA+J)
C     IF(NTYP(INA+J).EQ.NTM) IAE = IAB+NC(IAE)-1
      IS0 = IS
*
      DO 120 I=IAB,IAE
      CALL DAENCW(IE1(NC(I)),IE2(NC(I)),JJ)
      DO 110 K=1,NVMAX
      IF(K.NE.IV1.AND.K.NE.IV2.AND.JJ(K).NE.0) GOTO 120
 110  CONTINUE
      IS = IS+1
      CC(IS) = CC(I)
      NC(IS) = JJ(IV1)+IV2F*100*JJ(IV2)
 120  CONTINUE
*
      IF(RTM(J).NE.0.D0) THEN
         IS = IS+1
         NC(IS) = -1
         CC(IS) = RTM(J)
      ENDIF
*
      IF(IS-IS0.GE.1000) THEN
         PRINT*,'$$$ ERROR IN GRPOLY, TOO MANY COEFFICIENTS'
         CALL FOXDEB
      ENDIF
      IPLEN(J) = IS-IS0
 130  CONTINUE
*
*     COLOR POLYNOMIALS
*
      IF(NCARI.GT.0) GOTO 170
*
      DO 160 J=1,LB
      IBB = NBEG(INB+J)
      IBE = NEND(INB+J)
C     IF(NTYP(INB+J).EQ.NTM) IBE = IBB+NC(IBE)-1
      IS0 = IS
*
      DO 150 I=IBB,IBE
      CALL DAENCW(IE1(NC(I)),IE2(NC(I)),JJ)
      DO 140 K=1,NVMAX
      IF(K.NE.IV1.AND.K.NE.IV2.AND.JJ(K).NE.0) GOTO 150
 140  CONTINUE
      IS = IS+1
      CC(IS) = CC(I)
      NC(IS) = JJ(IV1)+IV2F*100*JJ(IV2)
 150  CONTINUE
*
      IF(RTM(J+3).NE.0.D0) THEN
         IS = IS+1
         NC(IS) = -1
         CC(IS) = RTM(J+3)
      ENDIF
*
      IF(IS-IS0.GE.1000) THEN
         PRINT*,'$$$ ERROR IN GRPOLY, TOO MANY COEFFICIENTS'
         CALL FOXDEB
      ENDIF
      IPLEN(J+3) = IS-IS0
 160  CONTINUE
*
      IF(LB.EQ.3) THEN
         IPLEN(7) = IPLEN(7)+1
         IS = IS+1
         NC(IS) = 0
         CC(IS) = 1.D0
      ENDIF
*
 170  CONTINUE
      ID = NEND(IND)+1
      NC(ID) = 500000000+IS-ID+1
      CC(ID) = IPLEN(1)*1.D6 +IPLEN(2)*1.D3 +IPLEN(3)
      ID = ID+1
      NC(ID) = IOBJ*1000000+NP1*1000+NP2


      CC(ID) = IPLEN(4)*1.D9 +IPLEN(5)*1.D6 +IPLEN(6)*1.D3 +IPLEN(7)
*
      NEND(IND) = IS
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRPOLY')
*
      RETURN
      END
*
      SUBROUTINE GRCURV(IXF1,IXF2,IXF3,
     *                  ITI1,ITI2,ITI3,ITF1,ITF2,ITF3,IND)
*     ****************************************************
*
*     THIS SUBROUTINE APPENDS A CURVE TO GRAPHICS OBJECT IND
*     IXF1,IXF2,IXF3 ARE FINAL POINT COORDINATES
*     ITI1,ITI2,ITI3 ARE TANGENT VECTOR AT INITIAL POINT
*     ITF1,ITF2,ITF3 ARE TANGENT VECTOR AT FINAL POINT
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      PARAMETER(PI=3.141592653589793238462643383279502884197169399375D0)
      DOUBLE PRECISION XC(3),XF(3),TI(3),TF(3)
      DOUBLE PRECISION D(3,2),AP(3,2),AM(3,2),COEF(3,0:3)
      INTEGER IPLEN(3)
      LOGICAL LINE
*
      IF(NTYP(IXF1).NE.NRE) CALL FOXNTY(IXF1)
      IF(NTYP(IXF2).NE.NRE) CALL FOXNTY(IXF2)
      IF(NTYP(IXF3).NE.NRE) CALL FOXNTY(IXF3)
      IF(NTYP(ITI1).NE.NRE) CALL FOXNTY(ITI1)
      IF(NTYP(ITI2).NE.NRE) CALL FOXNTY(ITI2)
      IF(NTYP(ITI3).NE.NRE) CALL FOXNTY(ITI3)
      IF(NTYP(ITF1).NE.NRE) CALL FOXNTY(ITF1)
      IF(NTYP(ITF2).NE.NRE) CALL FOXNTY(ITF2)
      IF(NTYP(ITF3).NE.NRE) CALL FOXNTY(ITF3)
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      TI(1) = CC(NBEG(ITI1))
      TI(2) = CC(NBEG(ITI2))
      TI(3) = CC(NBEG(ITI3))
      TF(1) = CC(NBEG(ITF1))
      TF(2) = CC(NBEG(ITF2))
      TF(3) = CC(NBEG(ITF3))
*
      LINE = .FALSE.
      LINE = LINE.OR.(TI(1).EQ.0.D0.AND.TI(2).EQ.0.D0.AND.TI(3).EQ.0.D0)
     *           .OR.(TF(1).EQ.0.D0.AND.TF(2).EQ.0.D0.AND.TF(3).EQ.0.D0)
      IF(LINE) THEN
         CALL GRDRAW(IXF1,IXF2,IXF3,IND)
         RETURN
      ENDIF
*
      XF(1) = CC(NBEG(IXF1))
      XF(2) = CC(NBEG(IXF2))
      XF(3) = CC(NBEG(IXF3))
*
*     OBTAIN THE CURRENT POSITION XC()
      I = NBEG(IND)+NC(NBEG(IND))
      NCI = NC(I)/100000000
      IF(NCI.GE.1.AND.NCI.LE.4) THEN
         DO 10 J=1,3
         XC(J) = CC(I+J-1)
 10      CONTINUE
      ELSE
         DO 20 J=1,3
         XC(J) = 0.D0
 20      CONTINUE
         IF(NCI.EQ.5) CALL PLYPOS(I,XC)
      ENDIF
*
      LINE = LINE
     *   .OR.(XF(1).EQ.XC(1).AND.XF(2).EQ.XC(2).AND.XF(3).EQ.XC(3))
      IF(LINE) THEN
         CALL GRDRAW(IXF1,IXF2,IXF3,IND)
         RETURN
      ENDIF
*
      CALL FOXALL(IVI,1,3)
      CALL FOXALL(IVF,1,3)
      CALL FOXALL(IVD,1,3)
      CALL FOXALL(IRS,1,1)
      NTYP(IVI) = NVE
      NTYP(IVF) = NVE
      NTYP(IVD) = NVE
      ISI = NBEG(IVI)-1
      ISF = NBEG(IVF)-1
      ISD = NBEG(IVD)-1
      NEND(IVI) = ISI+3
      NEND(IVF) = ISF+3
      NEND(IVD) = ISD+3
      ISS = NBEG(IRS)
*
      DO 30 J=1,3
      CC(ISI+J) = TI(J)
      CC(ISF+J) = TF(J)
      CC(ISD+J) = XF(J)-XC(J)
 30   CONTINUE
*
      CALL VEDOT(IVD,IVD,IRS)
      RL = SQRT(CC(ISS))
      CALL VEUNIT(IVD,IVD)
      CALL VEUNIT(IVI,IVI)
      CALL VEDOT(IVD,IVI,IRS)
      RI = CC(ISS)
      CALL VEUNIT(IVF,IVF)
      CALL VEDOT(IVD,IVF,IRS)
      RT = CC(ISS)
      LINE = LINE
     *   .OR.(ABS(RI).GT.(1.D0-1.D-6).AND.ABS(RT).GT.(1.D0-1.D-6))
*
      CALL FOXDAL(IRS,1)
      CALL FOXDAL(IVD,1)
*
      IF(LINE) THEN
         CALL FOXDAL(IVF,1)
         CALL FOXDAL(IVI,1)
         CALL GRDRAW(IXF1,IXF2,IXF3,IND)
         RETURN
      ENDIF
*
      AI = ACOS(MAX(-1.D0,MIN(1.D0,RI)))
      AT = ACOS(MAX(-1.D0,MIN(1.D0,RT)))
*
*     CUBIC BEZIER POLYNOMIALS USING 4 SUPPORT POINTS XC(), D(,1), D(,2), XF().
*
      DO 40 J=1,3
      D(J,1) = XC(J)+CC(ISI+J)*0.1D0*(PI+AT)*RL
      D(J,2) = XF(J)-CC(ISF+J)*0.1D0*(PI+AI)*RL
      AP(J,1) = XF(J)+XC(J)
      AM(J,1) = XF(J)-XC(J)
      AP(J,2) = D(J,2)+D(J,1)
      AM(J,2) = D(J,2)-D(J,1)
      COEF(J,0) = 0.125D0*AP(J,1)+0.375D0*AP(J,2)
      COEF(J,3) = 0.125D0*AM(J,1)-0.375D0*AM(J,2)
      COEF(J,1) = 0.375D0*(AM(J,1)+AM(J,2))
      COEF(J,2) = 0.375D0*(AP(J,1)-AP(J,2))
 40   CONTINUE
*
      CALL FOXDAL(IVF,1)
      CALL FOXDAL(IVI,1)
*
      NC(NBEG(IND)) = NEND(IND)-NBEG(IND)+1
*
      IS = NEND(IND)+2
      DO 50 J=1,3
      IPLEN(J) = 0
      DO 50 K=0,3
      IF(COEF(J,K).NE.0.D0) THEN
         IPLEN(J) = IPLEN(J)+1
         IS = IS+1
         NC(IS) = K
         CC(IS) = COEF(J,K)
      ENDIF
 50   CONTINUE
*
      ID = NEND(IND)+1
      NC(ID) = 500000000+IS-ID+1
      CC(ID) = IPLEN(1)*1.D6 +IPLEN(2)*1.D3 +IPLEN(3)
      ID = ID+1
      NC(ID) = 1000000
      CC(ID) = 0.D0
*
      NEND(IND) = IS
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRCURV')
*
      RETURN
      END
*
      SUBROUTINE GRCHAR(INA,IND)
*     **************************
*
*     THIS SUBROUTINE APPENDS CHARACTERS TO THE GRAPHICS OBJECT IND.
*     STRING ELEMENT DENOTED BY INA. IND IS THE GRAPHICS OBJECT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NST) CALL FOXNTY(INA)
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      LA = NEND(INA)-NBEG(INA)+1
      IF(LA.LE.0) RETURN
*
      IF(LA.GT.1024) THEN
         PRINT*,'$$$ WARNING IN GRCHAR, LENGTH IS BEYOND 1024'
         LA = 1024
      ENDIF
*
      ID = NEND(IND)+1
      IS = ID-1
      DO 10 I=NBEG(INA),NBEG(INA)+LA-1
      IS = IS+1
      CC(IS) = DBLE(NC(I))
      NC(IS) = 0
 10   CONTINUE
*
      NC(ID) = 600000000+LA
*
      NEND(IND) = IS
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRCHAR')
*
      RETURN
      END
*
      SUBROUTINE GRCOLR(INA,IND)
*     **************************
*
*     THIS SUBROUTINE APPENDS A COLOR CHANGE TO GRAPHICS OBJECT IND.
*     INA IS THE COSY COLOR ID OR A VECTOR CONTAINING AN RGBA COLOR.
*     RGBA COLOR VALUES ARE BETWEEN 0 AND 1. A IS FOR ALPHA FOR OPACITY.
*     A=0: TRANSPARENT, INVISIBLE. A=1: OPAQUE, COLOR IS IN FULL STRENGTH.
*     THE DEFAULT SETTING IS R=G=B=0, A=1 (FULL BLACK).
*
*     THE CONVENTIONAL COSY INDEXED COLOR
*     1: BLACK, 2: BLUE, 3: RED, 4: YELLOW, 5: GREEN, 6: YELLOWISH GREEN,
*     7: CYAN, 8: MAGENTA, 9: NAVY, 10: BACKGROUND (WHITE)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION ACC(4)
*
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      DO 10 I=1,3
      ACC(I) = 0.D0
 10   CONTINUE
      ACC(4) = 1.D0
*
      IF(NTYP(INA).EQ.NRE) THEN
         I = NINT(CC(NBEG(INA)))
         IF(I.EQ.2) THEN
            ACC(2) = 0.2D0
            ACC(3) = 1.D0
         ELSEIF(I.EQ.3) THEN
            ACC(1) = 1.D0
         ELSEIF(I.EQ.4) THEN
            ACC(1) = 1.D0
            ACC(2) = 1.D0
         ELSEIF(I.EQ.5) THEN
            ACC(2) = 1.D0
         ELSEIF(I.EQ.6) THEN
            ACC(1) = 0.6D0
            ACC(2) = 0.9D0
            ACC(3) = 0.2D0
         ELSEIF(I.EQ.7) THEN
            ACC(2) = 1.D0
            ACC(3) = 1.D0
         ELSEIF(I.EQ.8) THEN
            ACC(1) = 1.D0
            ACC(3) = 1.D0
         ELSEIF(I.EQ.9) THEN
            ACC(2) = 0.2D0
            ACC(3) = 0.7D0
         ELSEIF(I.EQ.10) THEN
            ACC(1) = 1.D0
            ACC(2) = 1.D0
            ACC(3) = 1.D0
         ENDIF
      ELSEIF(NTYP(INA).EQ.NVE) THEN
         IA = NBEG(INA)
         LA = NEND(INA)-IA+1
         IF(LA.LT.3) THEN
            PRINT*,'$$$ ERROR IN GRCOLR, NOT ENOUGH RGB COMPONENTS'
            CALL FOXDEB
         ENDIF
*
         DO 20 I=1,MIN(4,LA)
         ACC(I) = MIN(MAX(0.D0,CC(IA+I-1)),1.D0)
 20      CONTINUE
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      ID = NEND(IND)+1
      CC(ID) = ACC(1)
      NC(ID) = 700000004
      ID = ID+1
      CC(ID) = ACC(2)
      NC(ID) = 0
      ID = ID+1
      CC(ID) = ACC(3)
      NC(ID) = 0
      ID = ID+1
      CC(ID) = ACC(4)
      NC(ID) = 0
*
      NEND(IND) = ID
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRCOLR')
*
      RETURN
      END
*
      SUBROUTINE GRWDTH(INA,IND)
*     **************************
*
*     THIS SUBROUTINE APPENDS A WIDHT CHANGE TO GRAPHICS OBJECT IND
*     INA IS THE WIDTH ID
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      ID = NEND(IND)+1
      CC(ID) = CC(NBEG(INA))
      NC(ID) = 800000001
*
      NEND(IND) = ID
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRWDTH')
*
      RETURN
      END
*
      SUBROUTINE GRSTYL(INA,IND)
*     **************************
*
*     THIS SUBROUTINE APPENDS A DRAWING STYLE CHANGE TO GRAPHICS OBJECT IND.
*     INA IS THE STYLE OPTION.
*
*     THE 2ND AND UP COMPONENTS OF VE INA CAN BE USED AS ADDITIONAL STYLE INFO.
*     THE 1ST COMPONENT: A FLAG FOR WIRE FRAME DRAWING FOR SURFACES.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
C     STYLE1 = 0.D0
C     STYLE2 = 0.D0
*
      IF(NTYP(INA).EQ.NRE) THEN
         IWIRE = NINT(CC(NBEG(INA)))
      ELSEIF(NTYP(INA).EQ.NVE) THEN
         LA = NEND(INA)-NBEG(INA)+1
         IWIRE = NINT(CC(NBEG(INA)))
C        STYLE1 = CC(NBEG(INA)+1)
C        IF(LA.GE.3) STYLE2 = CC(NBEG(INA)+2)
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      IDB = NEND(IND)+1
      ITEM = 1
      ID = IDB
      CC(ID) = IWIRE
*
*     ADD MORE STYLE OPTIONS
*
C     ITEM = 2
C     ID = ID+1
C     CC(ID) = STYLE1
C     NC(ID) = 0
C
C     ITEM = 3
C     ID = ID+1
C     CC(ID) = STYLE2
C     NC(ID) = 0
*
      NC(IDB) = 900000000+ITEM
*
      NEND(IND) = ID
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRSTYL')
*
      RETURN
      END
*
      SUBROUTINE GREPS(INA,IND)
*     *************************
*
*     THIS SUBROUTINE SETS DRAWING ERROR TOLERANCE TO GRAPHICS OBJECT IND.
*     THE SPACE TOLERANCE AND THE COLOR TOLERANCE CAN BE GIVEN IN INA
*     BY RE OR VE. IF 0 IS GIVEN, IT RESETS TO THE DEFAULT MODE.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      EPSCOL = 0.D0
*
      IF(NTYP(INA).EQ.NRE) THEN
         EPSXYZ = ABS(CC(NBEG(INA)))
      ELSEIF(NTYP(INA).EQ.NVE) THEN
         EPSXYZ = ABS(CC(NBEG(INA)))
         EPSCOL = ABS(CC(NBEG(INA)+1))
      ELSE
         CALL FOXNTY(INA)
      ENDIF
*
      IDB = NEND(IND)+1
      ID = IDB
      CC(ID) = EPSXYZ
      ID = ID+1
      CC(ID) = EPSCOL
      NC(ID) = 0
*
      NC(IDB) = 1000000002
*
      NEND(IND) = ID
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GREPS')
*
      RETURN
      END
*
      SUBROUTINE GRPROJ(INA,INB,IND)
*     ******************************
*
*     THIS SUBROUTINE SETS A ROTATION OF VIEW TO GRAPHICS OBJECT IND
*     INA AND INB ARE THE TWO ANGLES
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      ID = NEND(IND)+1
      CC(ID) = CC(NBEG(INA))
      NC(ID) = 1100000003
      ID = ID+1
      CC(ID) = CC(NBEG(INB))
      NC(ID) = 0
      ID = ID+1
      CC(ID) = 0.D0
      NC(ID) = 0
*
      NEND(IND) = ID
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRPROJ')
*
      RETURN
      END
*
      SUBROUTINE GRZOOM(IX1,IX2,IY1,IY2,IZ1,IZ2,IND)
*     **********************************************
*
*     THIS SUBROUTINE SETS A ZOOMING AREA SPECIFIED BY TWO POINTS
*     (X1,Y1,Z1) AND (X2,Y2,Z2) TO GRAPHICS OBJECT IND
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(IX1).NE.NRE) CALL FOXNTY(IX1)
      IF(NTYP(IX2).NE.NRE) CALL FOXNTY(IX2)
      IF(NTYP(IY1).NE.NRE) CALL FOXNTY(IY1)
      IF(NTYP(IY2).NE.NRE) CALL FOXNTY(IY2)
      IF(NTYP(IZ1).NE.NRE) CALL FOXNTY(IZ1)
      IF(NTYP(IZ2).NE.NRE) CALL FOXNTY(IZ2)
      IF(NTYP(IND).NE.NGR) CALL FOXNTY(IND)
*
      X1 = MIN(CC(NBEG(IX1)),CC(NBEG(IX2)))
      X2 = MAX(CC(NBEG(IX1)),CC(NBEG(IX2)))
      Y1 = MIN(CC(NBEG(IY1)),CC(NBEG(IY2)))
      Y2 = MAX(CC(NBEG(IY1)),CC(NBEG(IY2)))
      Z1 = MIN(CC(NBEG(IZ1)),CC(NBEG(IZ2)))
      Z2 = MAX(CC(NBEG(IZ1)),CC(NBEG(IZ2)))
*
      ID = NEND(IND)+1
      CC(ID) = X1
      NC(ID) = 1200000006
      ID = ID+1
      CC(ID) = Y1
      NC(ID) = 0
      ID = ID+1
      CC(ID) = Z1
      NC(ID) = 0
      ID = ID+1
      CC(ID) = X2
      NC(ID) = 0
      ID = ID+1
      CC(ID) = Y2
      NC(ID) = 0
      ID = ID+1
      CC(ID) = Z2
      NC(ID) = 0
*
      NEND(IND) = ID
      IF(NEND(IND).GT.NMAX(IND)) CALL FOXERV('GRZOOM')
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     COMPLEX OPERATIONS                                                      *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE CMCOP(INA,INB)
*     *************************
*
*     THIS SUBROUTINE COPIES COMPLEX NUMBERS INA INTO INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
*
      CC(IB  ) = CC(IA  )
      CC(IB+1) = CC(IA+1)
*
      NTYP(INB) = NCM
      NEND(INB) = IB+1
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CMCOP')
*
      RETURN
      END
*
      SUBROUTINE LCM(INA,INC)
*     ***********************
*
*     THIS SUBROUTINE RETURNS THE ALLOCATION LENGTH OF COMPLEX IN INC.
*     THE INCOMING ARGUMENT INA IS VOID.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      NTYP(INC) = NRE
      CC(NBEG(INC)) = 2.D0
*
      RETURN
      END
*
      SUBROUTINE CMACM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS TWO COMPLEX NUMBERS INA AND INB AND
*     STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)
*
      CC(IC  ) = CC(IA  )+CC(IB  )
      CC(IC+1) = CC(IA+1)+CC(IB+1)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMACM')
*
      RETURN
      END
*
      SUBROUTINE CMSCM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE SUBTRACTS TWO COMPLEX NUMBERS INA AND INB
*     AND STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)
*
      CC(IC  ) = CC(IA  )-CC(IB  )
      CC(IC+1) = CC(IA+1)-CC(IB+1)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMSCM')
*
      RETURN
      END
*
      SUBROUTINE CMMCM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MULTIPLIES TWO COMPLEX NUMBERS INA AND INB
*     AND STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
      B1 = CC(IB  )
      B2 = CC(IB+1)
*
      CC(IC  ) = A1*B1-A2*B2
      CC(IC+1) = A1*B2+A2*B1
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMMCM')
*
      RETURN
      END
*
      SUBROUTINE CMDCM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE DIVIDES TWO COMPLEX NUMBERS INA AND INB AND
*     STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
      B1 = CC(IB  )
      B2 = CC(IB+1)
*
      IF(B1.EQ.0.D0.AND.B2.EQ.0.D0) THEN
         PRINT*,'$$$ ERROR IN CMDCM, DIVISOR IS ZERO'
         CALL FOXDEB
      ENDIF
*
      DE = B1*B1+B2*B2
*
      CC(IC  ) = (A1*B1+A2*B2)/DE
      CC(IC+1) = (A2*B1-A1*B2)/DE
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMDCM')
*
      RETURN
      END
*
      SUBROUTINE CMARE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS A COMPLEX AND A REAL NUMBER INA AND INB AND
*     STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      CC(IC  ) = CC(IA  )+CC(NBEG(INB))
      CC(IC+1) = CC(IA+1)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMARE')
*
      RETURN
      END
*
      SUBROUTINE REACM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS A REAL AND A COMPLEX NUMBER INA AND INB AND
*     STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
*
      IB = NBEG(INB)
      IC = NBEG(INC)
*
      CC(IC  ) = CC(IB  )+CC(NBEG(INA))
      CC(IC+1) = CC(IB+1)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('REACM')
*
      RETURN
      END
*
      SUBROUTINE CMSRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE SUBTRACTS A COMPLEX AND A REAL NUMBER INA AND INB AND
*     STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      CC(IC  ) = CC(IA  )-CC(NBEG(INB))
      CC(IC+1) = CC(IA+1)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMSRE')
*
      RETURN
      END
*
      SUBROUTINE RESCM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE SUBTRACTS A REAL AND A COMPLEX NUMBER INA AND INB AND
*     STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
*
      IB = NBEG(INB)
      IC = NBEG(INC)
*
      CC(IC  ) = -CC(IB  )+CC(NBEG(INA))
      CC(IC+1) = -CC(IB+1)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('RESCM')
*
      RETURN
      END
*
      SUBROUTINE CMMRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MULTIPLIES A COMPLEX AND A REAL NUMBER INA AND INB AND
*     STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      B1 = CC(NBEG(INB))
*
      CC(IC  ) = CC(IA  )*B1
      CC(IC+1) = CC(IA+1)*B1
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMMRE')
*
      RETURN
      END
*
      SUBROUTINE REMCM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MULTIPLIES A REAL AND A COMPLEX NUMBER INA AND INB AND
*     STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
*
      IB = NBEG(INB)
      IC = NBEG(INC)
*
      A1 = CC(NBEG(INA))
*
      CC(IC  ) = CC(IB  )*A1
      CC(IC+1) = CC(IB+1)*A1
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('REMCM')
*
      RETURN
      END
*
      SUBROUTINE CMDRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE DIVIDES A COMPLEX AND A REAL NUMBER INA AND INB AND
*     STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      B1 = CC(NBEG(INB))
*
      IF(B1.EQ.0.D0) THEN
         PRINT*,'$$$ ERROR IN CMDRE, DIVISOR IS ZERO'
         CALL FOXDEB
      ENDIF
*
      CC(IC  ) = CC(IA  )/B1
      CC(IC+1) = CC(IA+1)/B1
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMDRE')
*
      RETURN
      END
*
      SUBROUTINE REDCM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE DIVIDES A REAL AND A COMPLEX NUMBER INA AND INB AND
*     STORES THE RESULT IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
*
      IB = NBEG(INB)
      IC = NBEG(INC)
*
      A1 = CC(NBEG(INA))
      B1 = CC(IB  )
      B2 = CC(IB+1)
*
      IF(B1.EQ.0.D0.AND.B2.EQ.0.D0) THEN
         PRINT*,'$$$ ERROR IN REDCM, DIVISOR IS ZERO'
         CALL FOXDEB
      ENDIF
*
      DE = B1*B1+B2*B2
*
      CC(IC  ) =  A1*B1/DE
      CC(IC+1) = -A1*B2/DE
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('REDCM')
*
      RETURN
      END
*
      SUBROUTINE IMUNIT(INA)
*     **********************
*
*     THIS SUBROUTINE TURNS INA INTO THE COMPLEX UNITY I
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IA = NBEG(INA)
*
      CC(IA  ) = 0.D0
      CC(IA+1) = 1.D0
*
      NTYP(INA) = NCM
      NEND(INA) = IA+1
      IF(NEND(INA).GT.NMAX(INA)) CALL FOXERV('IMUNIT')
*
      RETURN
      END
*
      SUBROUTINE CMCONJ(INA,INB)
*     **************************
*
*     THIS SUBROUTINE CONJUGATES THE COMPLEX NUMBER INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
*
      CC(IB  ) =  CC(IA  )
      CC(IB+1) = -CC(IA+1)
*
      NTYP(INB) = NCM
      NEND(INB) = IB+1
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CMCONJ')
*
      RETURN
      END
*
      SUBROUTINE CMCMPL(INA,INB)
*     **************************
*
*     THIS SUBROUTINE TURNS THE REAL INA OR THE VECTOR INA INTO COMPLEX INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NRE.AND.NTYP(INA).NE.NVE) CALL FOXNTY(INA)
*
      IB = NBEG(INB)
*
      CC(IB  ) = CC(NBEG(INA))
      CC(IB+1) = 0.D0
      IF(NTYP(INA).EQ.NVE) CC(IB+1) = CC(NBEG(INA)+1)
*
      NTYP(INB) = NCM
      NEND(INB) = IB+1
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CMCMPL')
*
      RETURN
      END
*
      SUBROUTINE CMRE(INA,INB)
*     ************************
*
*     THIS SUBROUTINE PICKS OUT THE REAL PART OF THE COMPLEX NUMBER INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      CC(NBEG(INB)) =  CC(NBEG(INA))
      NTYP(INB) = NRE
*
      RETURN
      END
*
      SUBROUTINE CMIM(INA,INB)
*     ************************
*
*     THIS SUBROUTINE PICKS OUT THE IMAGINARY PART OF THE COMPLEX NUMBER INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      CC(NBEG(INB)) =  CC(NBEG(INA)+1)
      NTYP(INB) = NRE
*
      RETURN
      END
*
      SUBROUTINE CMPRE(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PICKS OUT THE REAL PART OR THE IMAGINARY PART
*     OF THE COMPLEX NUMBER INA, AND STORES IT IN INC, ACCORDING TO
*     THE SPECIFICATION BY THE REAL NUMBER IN INB. 1: REAL, 2: IMAGINARY
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NRE) CALL FOXNTY(INB)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      NB = NINT(CC(NBEG(INB)))
*
      IF(NB.EQ.1) THEN
         CC(IC) = CC(IA)
      ELSEIF(NB.EQ.2) THEN
         CC(IC) = CC(IA+1)
      ELSE
         PRINT*,'$$$ WARNING IN CMPRE, WRONG INDEX'
         CC(IC) = 0.D0
      ENDIF
*
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE CMPVE(INA,INC)
*     *************************
*
*     THIS SUBROUTINE PICKS OUT THE REAL PART AND THE IMAGINARY PART
*     OF THE COMPLEX NUMBER INA, AND STORES IT IN THE VECTOR INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
      IF(NEND(INC).LE.IC) CALL FOXERV('CMPVE')
*
      CC(IC)   = CC(IA)
      CC(IC+1) = CC(IA+1)
*
      NEND(INC) = IC+1
      NTYP(INC) = NVE
*
      RETURN
      END
*
      SUBROUTINE CMSQR(INA,INC)
*     *************************
*
*     THIS SUBROUTINE SQUARES THE COMPLEX NUMBER INA TO INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
*
      CC(IC  ) = A1*A1-A2*A2
      CC(IC+1) = 2.D0*A1*A2
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMSQR')
*
      RETURN
      END
*
      SUBROUTINE CMSQRT(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE COMPLEX SQUARE ROOT
*
*     ACCORDING TO THE DEFINITION OF CSQRT IN FORTRAN 77,
*     THE REAL PART OF THE RESULT IS POSITIVE OR 0.
*     IF 0, THE IMAGINARY PART OF THE RESULT IS POSITIVE.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
      IF(A2.EQ.0.D0) THEN
         C1 = 0.D0
         C2 = 0.D0
         IF(A1.GT.0.D0) THEN
            C1 = SQRT(A1)
         ELSEIF(A1.LT.0.D0) THEN
            C2 = SQRT(-A1)
         ENDIF
      ELSE
         IF(A1.EQ.0.D0) THEN
            R = SQRT(ABS(A2))
         ELSE
            R = SQRT(SQRT(A1*A1+A2*A2))
         ENDIF
         PHI = 0.5D0*ATAN2(A2,A1)
         C1 = R*COS(PHI)
         C2 = R*SIN(PHI)
      ENDIF
*
      CC(IC  ) = C1
      CC(IC+1) = C2
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMSQRT')
*
      RETURN
      END
*
      SUBROUTINE CMEXP(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE COMPLEX EXPONENTIAL
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
*
      CC(IC  ) = EXP(A1)*COS(A2)
      CC(IC+1) = EXP(A1)*SIN(A2)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMEXP')
*
      RETURN
      END
*
      SUBROUTINE CMLOG(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE COMPLEX LOGARITHM
*
*     ACCORDING TO THE DEFINITION OF CLOG IN FORTRAN 77,
*     THE IMAGINARY PART OF THE RESULT IS IN (-PI,PI].
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
      IF(A1.EQ.0.D0.AND.A2.EQ.0.D0) THEN
         PRINT*,'$$$ ERROR IN CMLOG, LOGARITHM DOES NOT EXIST FOR ',INA
         CALL FOXDEB
      ENDIF
*
      CC(IC  ) = LOG(SQRT(A1*A1+A2*A2))
      CC(IC+1) = ATAN2(A2,A1)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMLOG')
*
      RETURN
      END
*
      SUBROUTINE CMSINE(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE COMPLEX SINE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
*
      CC(IC  ) = SIN(A1)*COSH(A2)
      CC(IC+1) = COS(A1)*SINH(A2)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMSINE')
*
      RETURN
      END
*
      SUBROUTINE CMCOSE(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE COMPLEX COSINE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
*
      CC(IC  ) =  COS(A1)*COSH(A2)
      CC(IC+1) = -SIN(A1)*SINH(A2)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMCOSE')
*
      RETURN
      END
*
      SUBROUTINE CMSINH(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE COMPLEX HYPERBOLIC SINE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
*
      CC(IC  ) = SINH(A1)*COS(A2)
      CC(IC+1) = COSH(A1)*SIN(A2)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMSINH')
*
      RETURN
      END
*
      SUBROUTINE CMCOSH(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE COMPLEX HYPERBOLIC COSINE
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
      IC = NBEG(INC)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
*
      CC(IC  ) = COSH(A1)*COS(A2)
      CC(IC+1) = SINH(A1)*SIN(A2)
*
      NTYP(INC) = NCM
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CMCOSH')
*
      RETURN
      END
*
      SUBROUTINE CMABS(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE ABSOLUTE VALUE OF THE COMPLEX NUMBER
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
*
      A1 = CC(IA  )
      A2 = CC(IA+1)
*
      CC(NBEG(INC)) = SQRT(A1*A1+A2*A2)
*
      NTYP(INC) = NRE
*
      RETURN
      END
*
      SUBROUTINE CMWERF(INA,INB)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE COMPLEX ERROR FUNCTION W OF COMPLEX INA
*
*     THE ALGORITHM IS BASED ON
*
*     W. GAUTSCHI, "ALGORITHM 363, COMPLEX ERROR FUNCTION",
*     COMM. ACM 12, P.635, 1969.
*
*     K.S. KOELBIG, "CERTIFICATION OF ALGORITHM 363",
*     COMM. ACM 15, PP.465-466, 1972.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      INTEGER CAPN,NU,N,NP1
      LOGICAL BLOG
      DOUBLE PRECISION LAMBDA
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      A0 = CC(NBEG(INA)  )
      B0 = CC(NBEG(INA)+1)
*
      A = ABS(A0)
      B = ABS(B0)
*
      IF(B.LT.4.29D0.AND.A.LT.5.33D0) THEN
         S = (1.D0-B/4.29D0)*SQRT(1.D0-A*A/28.41D0)
         H = 1.6D0*S
         H2 = 2.D0*H
         CAPN = 6+INT(23.D0*S)
         LAMBDA = H2**CAPN
         NU = 9+INT(21.D0*S)
      ELSE
         H = 0.D0
         CAPN = 0
         NU = 8
      ENDIF
*
      BLOG = H.EQ.0.D0.OR.LAMBDA.EQ.0.D0
      R1 = 0.D0
      R2 = 0.D0
      S1 = 0.D0
      S2 = 0.D0
*
      DO 10 N=NU,0,-1
      NP1 = N+1
      T1 = B+H+NP1*R1
      T2 = A-NP1*R2
      C = .5D0/(T1*T1+T2*T2)
      R1 = C*T1
      R2 = C*T2
      IF(H.GT.0.D0.AND.N.LE.CAPN) THEN
         T1 = LAMBDA+S1
         S1 = R1*T1-R2*S2
         S2 = R2*T1+R1*S2
         LAMBDA = LAMBDA/H2
      ENDIF
 10   CONTINUE
*
      IF(B.EQ.0.D0) THEN
         WR = DEXP(-A*A)
      ELSE
         WR = 1.12837916709551D0
         IF(BLOG) THEN
            WR = WR*R1
         ELSE
            WR = WR*S1
         ENDIF
      ENDIF
*
      WI = 1.12837916709551D0
      IF(BLOG) THEN
         WI = WI*R2
      ELSE
         WI = WI*S2
      ENDIF
*
      IF(B0.GE.0.D0.AND.A0.LT.0.D0) WI=-WI
      IF(B0.LT.0.D0) THEN
         E2 = 2.D0*DEXP(-A*A+B*B)
         WR = E2*DCOS(-2.D0*A*B)-WR
         WI = E2*DSIN(-2.D0*A*B)-WI
         IF(A0.GT.0.D0) WI=-WI
      ENDIF
*
      CC(NBEG(INB)  ) = WR
      CC(NBEG(INB)+1) = WI
*
      NTYP(INB) = NCM
      NEND(INB) = NBEG(INB)+1
      IF(NMAX(INB)-NBEG(INB).LT.1) CALL FOXERV('CMWERF')
*
      PRINT*,'$$$ WARNING IN CMWERF, THIS ROUTINE IS UNDER DEVELOPMENT.'
      PRINT*,'    CONTACT US AT BERZ@MSU.EDU'
*
      RETURN
      END
*
      SUBROUTINE CMPRI(INA,IUNIT)
*     ***************************
*
*     THIS SUBROUTINE PRINTS A COMPLEX NUMBER TO UNIT IUNIT
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
*
      IA = NBEG(INA)
*
      WRITE(IUNIT,'(A,2(G17.9E3,A))') ' (',CC(IA),',',CC(IA+1),')'
*
      RETURN
      END
*
      SUBROUTINE CMCSTO(A1,A2,FORM,C)
*     *******************************
*
*     THIS SUBROUTINE OUTPUTS THE COMPLEX NUMBER (A1,A2) WITH THE FORMAT FORM
*     TO THE STRING C
*
*     FORM IS ASSUMED :    (REFER TO RESCT)
*     - ENCLOSED BY ( ), AND THERE IS NO BLANK BETWEEN
*     - AFTER (, STARTS BY ONE OF IiFfEeDdGg FOLLOWED BY A NUMBER CHARACTER
*     - AFTERWARD UNTIL CLOSING ), HAS ONLY NUMBER CHARACTERS OR Ee OR .
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER FORM*20,C10*100,C20*100,C1*100,C2*100,C*200
*
      IL = ILAST(FORM,1,20)
      LPOI = IL
      IF(INDEX(FORM(1:IL),'.').NE.0) LPOI = INDEX(FORM(1:IL),'.')
      IF(LPOI.EQ.5) THEN
         READ(FORM(3:LPOI-1),'(I2)',ERR=70) IDGT
      ELSEIF(LPOI.EQ.4) THEN
         READ(FORM(3:LPOI-1),'(I1)',ERR=70) IDGT
      ELSE
         PRINT*,'$$$ ERROR IN CMCSTO, TOO MANY DIGITS WERE DEMANDED'
         GOTO 80
      ENDIF
      IF(IDGT.LE.0) GOTO 70
      IFORM = INDEX('Ii',FORM(2:2))
*
      C10 = '                                                  '//
     *      '                                                  '
      C20 = C10
      C = C10//C10
      CALL REWSTO(A1,FORM,IFORM,C10)
      CALL REWSTO(A2,FORM,IFORM,C20)
      C1 = C10
      C2 = C20
*
      WRITE(C,'(A)') '('//C1(1:IDGT)//','//C2(1:IDGT)//')'
*
      RETURN
*
 70   PRINT*,'$$$ ERROR IN CMCSTO, '//
     *       'SYSTEM CANNOT ACCEPT THE GIVEN FORMAT'
 80   PRINT*,'    FORMAT = '//FORM
      CALL FOXDEB
*
      END
*
*******************************************************************************
*                                                                             *
*     COMPLEX DA OPERATIONS                                                   *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE CDPRI(INA,IUNIT)
*     ***************************
*
*     THIS SUBROUTINE PRINTS THE COMPLEX DA VECTOR INA TO UNIT IUNIT.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      INTEGER JJ(LNV)
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
*
*     WRITE(IUNIT,'(A)') ' '
*     WRITE(IUNIT,'(2A)')  ' CD VECTOR'
*     WRITE(IUNIT,'(A)')   ' *********'
      IF(NEND(INA).LT.NBEG(INA)) THEN
         WRITE(IUNIT,'(A)') '   ALL COMPONENTS ZERO '
      ELSE
         WRITE(IUNIT,'(A)')
     *    '     I  COEFFICIENTS'//
     *    '                           ORDER EXPONENTS'
      ENDIF
      IOUT = 0
      DO 100 IOA = 0,NOMAX
      DO 100 II=NBEG(INA),NEND(INA),2
      NCIA = NC(II)
      IF(IEO(NCIA).NE.IOA) GOTO 100
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
      IOUT = IOUT+1
      WRITE(IUNIT,
     *      '(I6,2(1X,G22.16),I4,2X,5(2I2,1X),20(/54X,5(2I2,1X)))')
     *      IOUT,CC(II),CC(II+1),IOA,(JJ(III),III=1,NVMAX)
*
 100  CONTINUE
*
      WRITE(IUNIT,'(A)') '                                      '
*
      RETURN
      END
*
      SUBROUTINE CDCOP(INA,INB)
*     *************************
*
*     THIS SUBROUTINE COPIES THE CD VECTOR A TO THE CD VECTOR B
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IB = NBEG(INB)-2
*
      DO 100 IA = NBEG(INA),NEND(INA),2
      IF(IEO(NC(IA)).GT.NOCUT) GOTO 100
      IB = IB+2
      CC(IB  ) = CC(IA  )
      NC(IB  ) = NC(IA  )
      CC(IB+1) = CC(IA+1)
      NC(IB+1) = NC(IA+1)
 100  CONTINUE
*
      NTYP(INB) = NCD
      NEND(INB) = IB+1
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CDCOP')
*
      RETURN
      END
*
      SUBROUTINE LCD(INA,INC)
*     ***********************
*
*     THIS SUBROUTINE RETURNS THE ALLOCATION LENGTH OF CD VECTOR IN INC.
*     THE VECTOR INA SUPPLIES THE ORDER AND THE DIMENSION OF CD.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      CALL LDA(INA,INC)
      CC(NBEG(INC)) = 2.D0*CC(NBEG(INC))
*
      RETURN
      END
*
      SUBROUTINE CDCNST(INA,INC)
*     **************************
*
*     THIS SUBROUTINE STORES THE CONSTANT PART OF CD VARIABLE INA IN INC
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
*
      NTYP(INC) = NCM
      NEND(INC) = NBEG(INC)+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDCNST')
      IC = NBEG(INC)
*
      IPOA = NBEG(INA)
*
      IF(NEND(INA).EQ.IPOA-1) THEN
         CC(IC  ) = 0.D0
         CC(IC+1) = 0.D0
         RETURN
      ELSE
         IF(NC(IPOA).EQ.1) THEN
            CC(IC  ) = CC(IPOA  )
            CC(IC+1) = CC(IPOA+1)
            RETURN
         ELSE
            CC(IC  ) = 0.D0
            CC(IC+1) = 0.D0
            RETURN
         ENDIF
      ENDIF
      END
*
      SUBROUTINE CDACD(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS TWO CD VECTORS INA AND INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCD) CALL FOXNTY(INB)
      IF((INA.EQ.INC).OR.(INB.EQ.INC)) CALL DANFI(INA,INB,INC)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)-2
      IAMAX = NEND(INA)
      IBMAX = NEND(INB)
      NA = NC(IA)
      NB = NC(IB)
*
      IF(IA.GT.IAMAX) THEN
         IC = IC+1
         DO 17 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = CC(IS)
  17     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDACD')
         RETURN
      ELSEIF(IB.GT.IBMAX) THEN
         IC = IC+1
         DO 18 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  18     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDACD')
         RETURN
      ENDIF
*
      IF(NA-NB) 30,20,40
*
*     ADDING TWO TERMS
*     ****************
*
  20  CONTINUE
      CC1 = CC(IA  )+CC(IB  )
      CC2 = CC(IA+1)+CC(IB+1)
      IF(ABS(CC1)+ABS(CC2).LT.EPS) GOTO 25
      IC = IC+2
      CC(IC  ) = CC1
      CC(IC+1) = CC2
      NC(IC  ) = NA
      NC(IC+1) = 0
  25  CONTINUE
      IA = IA+2
      IB = IB+2
      IF(IA.GT.IAMAX) THEN
         IC = IC+1
         DO 27 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = CC(IS)
  27     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDACD')
         RETURN
      ELSEIF(IB.GT.IBMAX) THEN
         IC = IC+1
         DO 28 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  28     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDACD')
         RETURN
      ENDIF
      NA = NC(IA)
      NB = NC(IB)
      IF(NA-NB) 30,20,40
*
*     STORING TERM A
*     **************
*
  30  CONTINUE
      IC = IC+2
      CC(IC  ) = CC(IA  )
      CC(IC+1) = CC(IA+1)
      NC(IC  ) = NA
      NC(IC+1) = 0
      IA = IA+2
      IF(IA.GT.IAMAX) THEN
         IC = IC+1
         DO 35 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = CC(IS)
  35     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDACD')
         RETURN
      ENDIF
      NA = NC(IA)
      IF(NA-NB) 30,20,40
*
*     STORING TERM B
*     **************
*
  40  CONTINUE
      IC = IC+2
      CC(IC  ) = CC(IB  )
      CC(IC+1) = CC(IB+1)
      NC(IC  ) = NB
      NC(IC+1) = 0
      IB = IB+2
      IF(IB.GT.IBMAX) THEN
         IC = IC+1
         DO 45 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  45     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDACD')
         RETURN
      ENDIF
      NB = NC(IB)
      IF(NA-NB) 30,20,40
*
      END
*%%
      SUBROUTINE CDSCD(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE SUBTRACTS TWO CD VECTORS INA AND INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCD) CALL FOXNTY(INB)
      IF((INA.EQ.INC).OR.(INB.EQ.INC)) CALL DANFI(INA,INB,INC)
*
      IA = NBEG(INA)
      IB = NBEG(INB)
      IC = NBEG(INC)-2
      IAMAX = NEND(INA)
      IBMAX = NEND(INB)
      NA = NC(IA)
      NB = NC(IB)
*
      IF(IA.GT.IAMAX) THEN
         IC = IC+1
         DO 17 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = -CC(IS)
  17     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDSCD')
         RETURN
      ELSEIF(IB.GT.IBMAX) THEN
         IC = IC+1
         DO 18 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  18     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDSCD')
         RETURN
      ENDIF
*
      IF(NA-NB) 30,20,40
*
*     SUBTRACTING TWO TERMS
*     *********************
*
  20  CONTINUE
      CC1 = CC(IA  )-CC(IB  )
      CC2 = CC(IA+1)-CC(IB+1)
      IF(ABS(CC1)+ABS(CC2).LT.EPS) GOTO 25
      IC = IC+2
      CC(IC  ) = CC1
      CC(IC+1) = CC2
      NC(IC  ) = NA
      NC(IC+1) = 0
  25  CONTINUE
      IA = IA+2
      IB = IB+2
      IF(IA.GT.IAMAX) THEN
         IC = IC+1
         DO 27 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = -CC(IS)
  27     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDSCD')
         RETURN
      ELSEIF(IB.GT.IBMAX) THEN
         IC = IC+1
         DO 28 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  28     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDSCD')
         RETURN
      ENDIF
      NA = NC(IA)
      NB = NC(IB)
      IF(NA-NB) 30,20,40
*
*     STORING TERM A
*     **************
*
  30  CONTINUE
      IC = IC+2
      CC(IC  ) = CC(IA  )
      CC(IC+1) = CC(IA+1)
      NC(IC  ) = NA
      NC(IC+1) = 0
      IA = IA+2
      IF(IA.GT.IAMAX) THEN
         IC = IC+1
         DO 35 IS=IB,IBMAX
         IC = IC+1
         CC(IC) = -CC(IS)
  35     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDSCD')
         RETURN
      ENDIF
      NA = NC(IA)
      IF(NA-NB) 30,20,40
*
*     STORING TERM B
*     **************
*
  40  CONTINUE
      IC = IC+2
      CC(IC  ) = -CC(IB  )
      CC(IC+1) = -CC(IB+1)
      NC(IC  ) = NB
      NC(IC+1) = 0
      IB = IB+2
      IF(IB.GT.IBMAX) THEN
         IC = IC+1
         DO 45 IS=IA,IAMAX
         IC = IC+1
         CC(IC) = CC(IS)
  45     NC(IC) = NC(IS)
         NTYP(INC) = NCD
         NEND(INC) = IC
         IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDSCD')
         RETURN
      ENDIF
      NB = NC(IB)
      IF(NA-NB) 30,20,40
*
      END
*%%
      SUBROUTINE CDMCD(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE PERFORMS A CD MULTIPLICATION OF THE CD VECTORS A AND B.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION NBNO(0:LNO),NENO(0:LNO)
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCD) CALL FOXNTY(INB)
*
      CALL CDCLR
*
*     RESORTING B BY ORDER
*     ********************
*
      JMEM = IMEM
      DO 20 I=0,NOCUT
      NBNO(I) = JMEM+1
      NENO(I) = JMEM-1
      JMEM = JMEM+2*NMMAX
  20  CONTINUE
*
      CALL MEMCHK(JMEM)
*
      DO 40 IB=NBEG(INB),NEND(INB),2
*
      NCIB = NC(IB)
      NOIB = IEO(NCIB)
      IPOS = NENO(NOIB)+2
      NENO(NOIB) = IPOS
*
      CC(IPOS  ) = CC(IB  )
      CC(IPOS+1) = CC(IB+1)
      NC(IPOS  ) = NC(IB  )
      NC(IPOS+1) = 0
*
  40  CONTINUE
*
*     PERFORMING ACTUAL MULTIPLICATION
*     ********************************
*
      DO 100 IA=NBEG(INA),NEND(INA),2
*
      NCIA = NC(IA)
      I1IA = IE1(NCIA)
      I2IA = IE2(NCIA)
      CCIR = CC(IA  )
      CCII = CC(IA+1)
*
      DO 100 NOIB = 0,NOCUT-IEO(NCIA)
      DO 100 IB = NBNO(NOIB),NENO(NOIB),2
*
      NCIB = NC(IB)
      ICC = 2*(IA2(I2IA+IE2(NCIB))+IA1(I1IA+IE1(NCIB)))-1
      CDA(ICC  ) = CDA(ICC  )+CCIR*CC(IB  )-CCII*CC(IB+1)
      CDA(ICC+1) = CDA(ICC+1)+CCIR*CC(IB+1)+CCII*CC(IB  )
*
 100  CONTINUE
*
      CALL CDPAC(INC)
*
      RETURN
      END
*
      SUBROUTINE CDDCD(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE DIVIDES THE CD VECTORS INA AND INB.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(IDIV,1,2*NMMAX)
      CALL CDMUI(INB,IDIV)
      CALL CDMCD(INA,IDIV,INC)
      CALL FOXDAL(IDIV,1)
*
      RETURN
      END
*
      SUBROUTINE CMACD(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS THE CM INA AND THE CD VECTOR INB.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL CDACM(INB,INA,INC)
      RETURN
      END
*
      SUBROUTINE CDACM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE ADDS THE CD VECTOR INA AND THE CM INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(INA.EQ.INC)       CALL DANFI(INA,INB,INC)
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
*
      IBEA = NBEG(INA)
      IENA = NEND(INA)
      IBEC = NBEG(INC)
*
      IF(IENA.LT.IBEA) THEN
         IF(ABS(CC(NBEG(INB)))+ABS(CC(NBEG(INB)+1)).GT.EPS) THEN
            CC(IBEC  ) = CC(NBEG(INB)  )
            CC(IBEC+1) = CC(NBEG(INB)+1)
            NC(IBEC  ) = 1
            NC(IBEC+1) = 0
            NEND(INC) = IBEC+1
            NTYP(INC) = NCD
            RETURN
         ELSE
            NEND(INC) = IBEC-1
            NTYP(INC) = NCD
            RETURN
         ENDIF
      ELSEIF(NC(IBEA).EQ.1) THEN
         CCR = CC(IBEA  )+CC(NBEG(INB)  )
         CCI = CC(IBEA+1)+CC(NBEG(INB)+1)
         IBBA = IBEA+2
      ELSE
         CCR = CC(NBEG(INB)  )
         CCI = CC(NBEG(INB)+1)
         IBBA = IBEA
      ENDIF
*
      IF(ABS(CCR)+ABS(CCI).GT.EPS) THEN
         CC(IBEC  ) = CCR
         CC(IBEC+1) = CCI
         NC(IBEC  ) = 1
         NC(IBEC+1) = 0
         IC = IBEC+1
      ELSE
         IC = IBEC-1
      ENDIF
*
      DO 10 I=IBBA,IENA
      IC = IC+1
      CC(IC) = CC(I)
      NC(IC) = NC(I)
  10  CONTINUE
*
      IF(IC.GT.NMAX(INC)) CALL FOXERV('CDACM')
*
      NEND(INC) = IC
      NTYP(INC) = NCD
*
      RETURN
      END
*
      SUBROUTINE CDSCM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE SUBTRACTS THE CM INB FROM THE CD VECTOR INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
      ICC = NBEG(INB)
      CC(ICC  ) = -CC(ICC  )
      CC(ICC+1) = -CC(ICC+1)
      CALL CDACM(INA,INB,INC)
      CC(ICC  ) = -CC(ICC  )
      CC(ICC+1) = -CC(ICC+1)
*
      RETURN
      END
*
      SUBROUTINE CMSCD(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE SUBTRACTS THE CD VECTOR INB FROM THE CM INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCM) CALL FOXNTY(INA)
      ICC = NBEG(INA)
      CC(ICC  ) = -CC(ICC  )
      CC(ICC+1) = -CC(ICC+1)
      CALL CDACM(INB,INA,INC)
      CALL CDADI(INC,INC)
      CC(ICC  ) = -CC(ICC  )
      CC(ICC+1) = -CC(ICC+1)
*
      RETURN
      END
*%%
      SUBROUTINE CDADI(INA,INB)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE ADDITIVE INVERSE OF INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
*
      IB = NBEG(INB)-1
*
      DO 100 IA = NBEG(INA),NEND(INA)
      IB = IB+1
      CC(IB) = -CC(IA)
      NC(IB) = NC(IA)
 100  CONTINUE
*
      NTYP(INB) = NCD
      NEND(INB) = IB
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CDADI')
*
      RETURN
      END
*
      SUBROUTINE CMMCD(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL CDMCM(INB,INA,INC)
*
      RETURN
      END
*
      SUBROUTINE CDMCM(INA,INB,INC)
*     *****************************
*
*     THIS SUBROUTINE MULTIPLIES CD VECTOR INA WITH CM INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
      CKOR = CC(NBEG(INB)  )
      CKOI = CC(NBEG(INB)+1)
*
      IC = NBEG(INC)-2
*
      DO 100 IA = NBEG(INA),NEND(INA),2
      CCR = CC(IA  )*CKOR-CC(IA+1)*CKOI
      CCI = CC(IA+1)*CKOR+CC(IA  )*CKOI
      IF(ABS(CCR)+ABS(CCI).LT.EPS) GOTO 100
      IC = IC+2
      CC(IC  ) = CCR
      CC(IC+1) = CCI
      NC(IC  ) = NC(IA  )
      NC(IC+1) = NC(IA+1)
 100  CONTINUE
*
      NTYP(INC) = NCD
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDMCM')
*
      RETURN
      END
*%%
      SUBROUTINE CDDCM(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE DIVIDES THE CD VECTOR INA BY THE CM INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INB).NE.NCM) CALL FOXNTY(INB)
      ICC = NBEG(INB)
      CCR = CC(ICC  )
      CCI = CC(ICC+1)
      IF(ABS(CCR)+ABS(CCI).EQ.0.D0) THEN
         PRINT*, '$$$ ERROR IN CDDCM, DIVISOR IS ZERO'
         CALL FOXDEB
      ENDIF
*
      CC(ICC  ) =  CCR/(CCR*CCR+CCI*CCI)
      CC(ICC+1) = -CCI/(CCR*CCR+CCI*CCI)
      CALL CDMCM(INA,INB,INC)
      CC(ICC  ) = CCR
      CC(ICC+1) = CCI
*
      RETURN
      END
*
      SUBROUTINE CMDCD(INA,INB,INC)
*     ******************************
*
*     THIS SUBROUTINE DIVIDES THE CM INA BY THE CD VECTOR INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL CDMUI(INB,INC)
      CALL CDMCM(INC,INA,INC)
*
      RETURN
      END
*%%
      SUBROUTINE REACD(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL FOXALL(ICM,1,2)
      CALL CMCMPL(INA,ICM)
      CALL CMACD (ICM,INB,INC)
      CALL FOXDAL(ICM,1)
      RETURN
      END
*
      SUBROUTINE RESCD(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL FOXALL(ICM,1,2)
      CALL CMCMPL(INA,ICM)
      CALL CMSCD (ICM,INB,INC)
      CALL FOXDAL(ICM,1)
      RETURN
      END
*
      SUBROUTINE REMCD(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL FOXALL(ICM,1,2)
      CALL CMCMPL(INA,ICM)
      CALL CMMCD (ICM,INB,INC)
      CALL FOXDAL(ICM,1)
      RETURN
      END
*
      SUBROUTINE REDCD(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL FOXALL(ICM,1,2)
      CALL CMCMPL(INA,ICM)
      CALL CMDCD (ICM,INB,INC)
      CALL FOXDAL(ICM,1)
      RETURN
      END
*
      SUBROUTINE CDARE(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL FOXALL(ICM,1,2)
      CALL CMCMPL(INB,ICM)
      CALL CDACM (INA,ICM,INC)
      CALL FOXDAL(ICM,1)
      RETURN
      END
*
      SUBROUTINE CDSRE(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL FOXALL(ICM,1,2)
      CALL CMCMPL(INB,ICM)
      CALL CDSCM (INA,ICM,INC)
      CALL FOXDAL(ICM,1)
      RETURN
      END
*
      SUBROUTINE CDMRE(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL FOXALL(ICM,1,2)
      CALL CMCMPL(INB,ICM)
      CALL CDMCM (INA,ICM,INC)
      CALL FOXDAL(ICM,1)
      RETURN
      END
*
      SUBROUTINE CDDRE(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL FOXALL(ICM,1,2)
      CALL CMCMPL(INB,ICM)
      CALL CDDCM (INA,ICM,INC)
      CALL FOXDAL(ICM,1)
      RETURN
      END
*
      SUBROUTINE DAACM(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INA,ICD)
      CALL CDACM (ICD,INB,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE DASCM(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INA,ICD)
      CALL CDSCM (ICD,INB,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE DAMCM(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INA,ICD)
      CALL CDMCM (ICD,INB,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE DADCM(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INA,ICD)
      CALL CDDCM (ICD,INB,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE CMADA(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INB,ICD)
      CALL CMACD (INA,ICD,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE CMSDA(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INB,ICD)
      CALL CMSCD (INA,ICD,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE CMMDA(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INB,ICD)
      CALL CMMCD (INA,ICD,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE CMDDA(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INB,ICD)
      CALL CMDCD (INA,ICD,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE DAACD(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INA,ICD)
      CALL CDACD (ICD,INB,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE DASCD(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INA,ICD)
      CALL CDSCD (ICD,INB,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE DAMCD(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INA,ICD)
      CALL CDMCD (ICD,INB,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE DADCD(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INA,ICD)
      CALL CDDCD (ICD,INB,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE CDADA(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INB,ICD)
      CALL CDACD (INA,ICD,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE CDSDA(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INB,ICD)
      CALL CDSCD (INA,ICD,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE CDMDA(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INB,ICD)
      CALL CDMCD (INA,ICD,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE CDDDA(INA,INB,INC)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL CDCMPL(INB,ICD)
      CALL CDDCD (INA,ICD,INC)
      CALL FOXDAL(ICD,1)
      RETURN
      END
*
      SUBROUTINE CDMUI(INA,INC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE MULTIPLICATIVE INVERSE OF INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
*
      CALL FOXALL(ICON,1,2*NMMAX)
      CALL FOXALL(IDEN,1,2*NMMAX)
      CALL FOXALL(ISCR,1,2*NMMAX)
      CALL FOXALL(INCR,1,2*NMMAX)
      CALL FOXALL(INCI,1,2*NMMAX)
      CALL FOXALL(IIUN,1,2)
      CALL IMUNIT(IIUN)
      CALL CDCONJ(INA ,ICON)
      CALL CDMCD (INA ,ICON,ISCR)
      CALL CDRE  (ISCR,IDEN)
      CALL CDRE  (INA ,ISCR)
      CALL DADDA (ISCR,IDEN,INCR)
      CALL CDIM  (INA ,ISCR)
      CALL DADDA (ISCR,IDEN,INCI)
      CALL CDCMPL(INCI,ISCR)
      CALL CDMCM (ISCR,IIUN,INCI)
      CALL CDCMPL(INCR,ISCR)
      CALL CDSCD (ISCR,INCI,INC )
      CALL FOXDAL(IIUN,1)
      CALL FOXDAL(INCI,1)
      CALL FOXDAL(INCR,1)
      CALL FOXDAL(ISCR,1)
      CALL FOXDAL(IDEN,1)
      CALL FOXDAL(ICON,1)
*
      RETURN
      END
*
      SUBROUTINE CDSQR(INA,INC)
*     *************************
*
*     THIS SUBROUTINE SQUARES INA
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL CDMCD(INA,INA,INC)
*
      RETURN
      END
*
      SUBROUTINE CDCLR
*     ****************
*
*     THIS SUBROUTINE SETS ALL THE SCRATCH SPACE IN CDA TO ZERO
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IF(LFLT.EQ.0) THEN
         DO 10 I=1,2*NMMAX
  10     CDA(I) = 0.D0
      ELSE
         DO 20 I=1,NFLT
         NCF2 = 2*NCFLT(I)
         CDA(NCF2-1) = 0.D0
  20     CDA(NCF2  ) = 0.D0
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE CDPAC(INC)
*     *********************
*
*     THIS SUBROUTINE PACKS THE INFORMATION IN THE SCRATCH VECTOR CDA
*     INTO THE CD VECTOR INC.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      IC = NBEG(INC)-2
*
      IF(LFLT.EQ.0) THEN
         DO 100 I=1,2*NMMAX,2
         CCR = CDA(I  )
         CCI = CDA(I+1)
         IF(ABS(CCR)+ABS(CCI).LT.EPS) GOTO 100
         IC = IC+2
         CC(IC  ) = CCR
         CC(IC+1) = CCI
         NC(IC  ) = (I+1)/2
         NC(IC+1) = 0
 100     CONTINUE
      ELSE
         DO 200 I=1,NFLT
         NCF  = NCFLT(I)
         NCF2 = 2*NCF
         CCR = CDA(NCF2-1)
         CCI = CDA(NCF2)
         IF(ABS(CCR)+ABS(CCI).LT.EPS) GOTO 200
         IC = IC+2
         CC(IC  ) = CCR
         CC(IC+1) = CCI
         NC(IC  ) = NCF
         NC(IC+1) = 0
 200     CONTINUE
      ENDIF
*
      NTYP(INC) = NCD
      NEND(INC) = IC+1
      IF(NEND(INC).GT.NMAX(INC)) CALL FOXERV('CDPAC')
*
      RETURN
      END
*%%
      SUBROUTINE CDRE(INA,INB)
*     ************************
*
*     THIS SUBROUTINE EXTRACTS THE REAL PART OF THE CD VECTOR A TO THE DA
*     VECTOR B
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IB = NBEG(INB)-1
*
      DO 100 IA = NBEG(INA),NEND(INA),2
      IF(CC(IA).EQ.0.D0) GOTO 100
      IB = IB+1
      CC(IB) = CC(IA)
      NC(IB) = NC(IA)
 100  CONTINUE
*
      NTYP(INB) = NDA
      NEND(INB) = IB
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CDRE')
*
      RETURN
      END
*
      SUBROUTINE CDIM(INA,INB)
*     ************************
*
*     THIS SUBROUTINE EXTRACTS THE IMAG. DA PART OF THE CD VECTOR A TO THE
*     DA VECTOR B
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IB = NBEG(INB)-1
*
      DO 100 IA = NBEG(INA),NEND(INA),2
      IF(CC(IA+1).EQ.0.D0) GOTO 100
      IB = IB+1
      CC(IB) = CC(IA+1)
      NC(IB) = NC(IA)
 100  CONTINUE
*
      NTYP(INB) = NDA
      NEND(INB) = IB
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CDIM')
*
      RETURN
      END
*
      SUBROUTINE CDCONJ(INA,INB)
*     **************************
*
*     THIS SUBROUTINE CONJUGATES THE CD VECTOR A TO THE CD VECTOR B
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      IB = NBEG(INB)-2
*
      DO 100 IA = NBEG(INA),NEND(INA),2
      IB = IB+2
      CC(IB  ) =  CC(IA  )
      CC(IB+1) = -CC(IA+1)
      NC(IB  ) =  NC(IA  )
      NC(IB+1) =  0
 100  CONTINUE
*
      NTYP(INB) = NCD
      NEND(INB) = IB+1
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CDCONJ')
*
      RETURN
      END
*
      SUBROUTINE CDCMPL(INA,INB)
*     **************************
*
*     THIS SUBROUTINE TURNS DA INA INTO CD VECTOR INB
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(INA.EQ.INB) CALL DANFI(INA,INB,0)
      IF(NTYP(INA).NE.NDA) CALL FOXNTY(INA)
      IB = NBEG(INB)-2
*
      DO 100 IA = NBEG(INA),NEND(INA)
      IB = IB+2
      CC(IB  ) =  CC(IA  )
      CC(IB+1) =  0.D0
      NC(IB  ) =  NC(IA  )
      NC(IB+1) =  0
 100  CONTINUE
*
      NTYP(INB) = NCD
      NEND(INB) = IB+1
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CDCMPL')
*
      RETURN
      END
*
      SUBROUTINE RECD(IIV,INC)
*     ************************
*
*     THIS SUBROUTINE CREATES AN IDENTITY CD VECTOR INC OF IV-TH
*     INDEPENDENT VARIABLE. THE IMAGINARY PART OF INC IS 0.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL FOXALL(ISC,1,1)
      CALL DA(IIV,ISC)
      CALL CDCMPL(ISC,INC)
      CALL FOXDAL(ISC,1)
*
      RETURN
      END
*
      SUBROUTINE CDNF(INA,IMU1,IMU2,IMU3,IIRES,IJRES,INRES,INB)
*     *********************************************************
*
*     THIS SUBROUTINE MULTIPLIES EACH MONOMIAL IN INA WITH
*     1 / (1 - EXP(I*( (N1-N2)*MU1 + (N3-N4)*MU2 + (N5-N6)*MU3) ) ) ),
*     EXCEPT WHEN THE TRIPLET (N1-N2), (N3-N4), (N5-N6) OCCURS IN IRES.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION IRES(3,100)
      INTEGER JJ(LNV)
*
      IF(LEW.EQ.1) THEN
         PRINT*,'$$$ ERROR IN CDNF, '//
     *          'WEIGHTED DA (DANOTW) IS NOT SUPPORTED'
         CALL FOXDEB
      ENDIF
*
      IF(NTYP(INA).NE.NCD)   CALL FOXNTY(INA)
      IF(NTYP(IMU1).NE.NRE) CALL FOXNTY(IMU1)
      TMU1 = CC(NBEG(IMU1))
      IF(NTYP(IMU2).NE.NRE) CALL FOXNTY(IMU2)
      TMU2 = CC(NBEG(IMU2))
      IF(NTYP(IMU3).NE.NRE) CALL FOXNTY(IMU3)
      TMU3 = CC(NBEG(IMU3))
      IF(NTYP(IJRES).NE.NRE) CALL FOXNTY(IJRES)
      IF(NINT(CC(NBEG(IJRES))).NE.3) THEN
         PRINT*,'$$$ ERROR IN CDNF, WRONG FIRST DIM OF RESONANCE ARRAY'
         CALL FOXDEB
      ENDIF
      IF(NTYP(INRES).NE.NRE) CALL FOXNTY(INRES)
      NRES = NINT(CC(NBEG(INRES)))
*
      DO 10 I=1,NRES
      DO 10 J=1,3
      IAD = NBEG(IIRES+1)-4+J+3*I
      IF(NTYP(IAD).NE.NRE) CALL FOXNTY(IAD)
      IRES(J,I) = NINT(CC(NBEG(IAD)))
      PRINT*,'J, I, IRES: ',J,I,IRES(J,I)
  10  CONTINUE
*
      IB = NBEG(INB)-2
*
      DO 50 IA=NBEG(INA),NEND(INA),2
*
      NCIA = NC(IA)
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
      N1 = JJ(1)-JJ(2)
      N2 = JJ(3)-JJ(4)
      N3 = JJ(5)-JJ(6)
*
      DO 20 I=1,NRES
      IF((N1.EQ.IRES(1,I)).AND.(N2.EQ.IRES(2,I)).AND.(N3.EQ.IRES(3,I)))
     *    GOTO 50
  20  CONTINUE
*
      PHI = N1*TMU1+N2*TMU2+N3*TMU3
*
      CCRI = 1.D0-COS(PHI)
      CCII =-SIN(PHI)
      ABS2 = CCRI*CCRI+CCII*CCII
      IF(ABS2.LT.EPS*EPS) GOTO 50
*
      CCR  =  CCRI/ABS2
      CCI  = -CCII/ABS2
      IB = IB+2
      CC(IB  ) = CC(IA  )*CCR-CC(IA+1)*CCI
      CC(IB+1) = CC(IA+1)*CCR+CC(IA  )*CCI
      NC(IB  ) = NCIA
      NC(IB+1) = 0
*
  50  CONTINUE
*
      NTYP(INB) = NCD
      NEND(INB) = IB+1
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CDNF')
*
      RETURN
      END
*
      SUBROUTINE CDF2(INA,IMU1,IMU2,IMU3,INB)
*     ***************************************
*
*     THIS SUBROUTINE MULTIPLIES EACH MONOMIAL IN INA WITH
*     EXP(I*( (N1-N2)*MU1 + (N3-N4)*MU2 + (N5-N6)*MU3) ) ) )
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      INTEGER JJ(LNV)
*
      IF(LEW.EQ.1) THEN
         PRINT*,'$$$ ERROR IN CDF2, '//
     *          'WEIGHTED DA (DANOTW) IS NOT SUPPORTED'
         CALL FOXDEB
      ENDIF
*
      IF(NTYP(INA).NE.NCD)   CALL FOXNTY(INA)
      IF(NTYP(IMU1).NE.NRE) CALL FOXNTY(IMU1)
      TMU1 = CC(NBEG(IMU1))
      IF(NTYP(IMU2).NE.NRE) CALL FOXNTY(IMU2)
      TMU2 = CC(NBEG(IMU2))
      IF(NTYP(IMU3).NE.NRE) CALL FOXNTY(IMU3)
      TMU3 = CC(NBEG(IMU3))
*
      IB = NBEG(INB)-2
*
      DO 50 IA=NBEG(INA),NEND(INA),2
*
      NCIA = NC(IA)
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
      N1 = JJ(1)-JJ(2)
      N2 = JJ(3)-JJ(4)
      N3 = JJ(5)-JJ(6)
*
      PHI = N1*TMU1+N2*TMU2+N3*TMU3
*
      CCR = COS(PHI)
      CCI = SIN(PHI)
*
      IB = IB+2
      CC(IB  ) = CC(IA  )*CCR-CC(IA+1)*CCI
      CC(IB+1) = CC(IA+1)*CCR+CC(IA  )*CCI
      NC(IB  ) = NCIA
      NC(IB+1) = 0
*
  50  CONTINUE
*
      NTYP(INB) = NCD
      NEND(INB) = IB+1
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CDF2')
*
      RETURN
      END
*
      SUBROUTINE CDNFDA(INA,IR,IT,IIC,IIM,IEPS,INB)
*     *********************************************
*
*     THIS SUBROUTINE MULTIPLIES EACH MONOMIAL IN INA WITH
*     1 / (R(IC)*EXP(I*MU(IC)) - PROD_I R(J)**N(J) EXP(I*VEC(MU)*VEC(N)) ) ),
*     EXCEPT WHEN THE NORM OF THE DENOMINATOR IS BELOW TEPS.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      LOGICAL FL1,FL2
      DIMENSION R(LNV),T(LNV)
      INTEGER JJ(LNV)
*
      IF(LEW.EQ.1) THEN
         PRINT*,'$$$ ERROR IN CDNFDA, '//
     *          'WEIGHTED DA (DANOTW) IS NOT SUPPORTED'
         CALL FOXDEB
      ENDIF
*
      IF(NTYP(INA).NE.NCD)   CALL FOXNTY(INA)
      IF(NTYP(IEPS).NE.NRE)  CALL FOXNTY(IEPS)
      TEPS = MAX(EPS,CC(NBEG(IEPS)))
      IF(NTYP(IIM).NE.NRE) CALL FOXNTY(IIM)
      IM = NINT(CC(NBEG(IIM)))
      IF(NTYP(IIC).NE.NRE) CALL FOXNTY(IIC)
      IC = NINT(CC(NBEG(IIC)))
*
      INR = NBEG(IR+1)-1
      INT = NBEG(IT+1)-1
      DO 10 J=1,IM
      IF(NTYP(IR+J).NE.NRE) CALL FOXNTY(IR+J)
      IF(NTYP(IT+J).NE.NRE) CALL FOXNTY(IT+J)
      R(J) = CC(NBEG(IR+J))
      T(J) = CC(NBEG(IT+J))
  10  CONTINUE
*
      IB = NBEG(INB)-2
      NEX = 0
      FL1 = .TRUE.
      DO 20 J=1,IM
  20  IF(R(J).NE.1.D0) FL1 = .FALSE.
      DO 21 J=2,IM,2
  21  IF(T(J-1).NE.-T(J)) FL1 = .FALSE.
*
      DO 50 IA=NBEG(INA),NEND(INA),2
      NCIA = NC(IA)
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
*
*      IF(FL1) THEN
*        FL2 = .TRUE.
*        DO 30 J=2,IM,2
*          N = JJ(J-1)-JJ(J)
*          IF(J-1.EQ.IC) THEN
*            IF(N.NE.1) FL2 = .FALSE.
*          ELSEIF(J.EQ.IC) THEN
*            IF(N.NE.-1) FL2 = .FALSE.
*          ELSE
*            IF(N.NE.0) FL2 = .FALSE.
*          ENDIF
*  30    CONTINUE
*        IF(FL2) GOTO 50
*      ENDIF
*
      RRR = 1.D0
      PHI = 0.D0
      DO 40 J=1,IM
      N = JJ(J)
      PHI = PHI+T(J)*N
      RRR = RRR*(R(J)**N)
  40  CONTINUE
*
      CCRI = R(IC)*COS(T(IC))-RRR*COS(PHI)
      CCII = R(IC)*SIN(T(IC))-RRR*SIN(PHI)
      ABS2 = CCRI*CCRI+CCII*CCII
*
      IF((ABS2.LT.TEPS*TEPS).AND.
     1((ABS2.GT.1D-16).OR.((ABS(CC(IA))+ABS(CC(IA+1))).GT.1E-8)))
     2 NEX = NEX+1
      IF(ABS2.LT.TEPS*TEPS) GOTO 50
*
      CCR  =  CCRI/ABS2
      CCI  = -CCII/ABS2
      IB = IB+2
      CC(IB  ) = CC(IA  )*CCR-CC(IA+1)*CCI
      CC(IB+1) = CC(IA+1)*CCR+CC(IA  )*CCI
      NC(IB  ) = NCIA
      NC(IB+1) = 0
*
  50  CONTINUE
*
*      IF(NEX.NE.0)
*     1WRITE(6,'(1X,A19,I7,A21)')
*     2'--> INFO IN CDNFDA,',NEX,' TERMS NOT ELIMINATED'
*
      NTYP(INB) = NCD
      NEND(INB) = IB+1
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CDNFDA')
*
      RETURN
      END
*%%
      SUBROUTINE CDNFDS(INA,IR,IT,ITS,IIM,IEPS,INB)
*     *********************************************
*
*     THIS SUBROUTINE MULTIPLIES EACH MONOMIAL IN INA WITH
*     1 / (1 - EXP(I*MUS) * PROD_I R(J)**N(J) EXP(I*VEC(MU)*VEC(N)) ) ),
*     EXCEPT WHEN THE NORM OF THE DENOMINATOR IS BELOW TEPS.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION R(LNV),T(LNV)
      INTEGER JJ(LNV)
*
      IF(LEW.EQ.1) THEN
         PRINT*,'$$$ ERROR IN CDNFDS, '//
     *          'WEIGHTED DA (DANOTW) IS NOT SUPPORTED'
         CALL FOXDEB
      ENDIF
*
      IF(NTYP(INA).NE.NCD)   CALL FOXNTY(INA)
      IF(NTYP(IEPS).NE.NRE)  CALL FOXNTY(IEPS)
      TEPS = MAX(EPS,CC(NBEG(IEPS)))
      IF(NTYP(IIM).NE.NRE) CALL FOXNTY(IIM)
      IM = NINT(CC(NBEG(IIM)))
      IF(NTYP(ITS).NE.NRE) CALL FOXNTY(ITS)
      TS = CC(NBEG(ITS))
*
      INR = NBEG(IR+1)-1
      INT = NBEG(IT+1)-1
      DO 10 J=1,IM
      IF(NTYP(IR+J).NE.NRE) CALL FOXNTY(IR+J)
      IF(NTYP(IT+J).NE.NRE) CALL FOXNTY(IT+J)
      R(J) = CC(NBEG(IR+J))
      T(J) = CC(NBEG(IT+J))
  10  CONTINUE
*
      IB = NBEG(INB)-2
*
      DO 50 IA=NBEG(INA),NEND(INA),2
*
      NCIA = NC(IA)
      CALL DAENC(IE1(NCIA),IE2(NCIA),JJ)
      RRR = 1.D0
      PHI = TS
      DO 40 J=1,IM
      N = JJ(J)
      PHI = PHI+T(J)*N
      RRR = RRR*(R(J)**N)
  40  CONTINUE
*
      CCRI = 1.D0-RRR*COS(PHI)
      CCII = 0.D0-RRR*SIN(PHI)
      ABS2 = CCRI*CCRI+CCII*CCII
      IF(ABS2.LT.TEPS*TEPS) GOTO 50
*
      CCR  =  CCRI/ABS2
      CCI  = -CCII/ABS2
      IB = IB+2
      CC(IB  ) = CC(IA  )*CCR-CC(IA+1)*CCI
      CC(IB+1) = CC(IA+1)*CCR+CC(IA  )*CCI
      NC(IB  ) = NCIA
      NC(IB+1) = 0
*
  50  CONTINUE
*
      NTYP(INB) = NCD
      NEND(INB) = IB+1
      IF(NEND(INB).GT.NMAX(INB)) CALL FOXERV('CDNFDS')
*
      RETURN
      END
*
      SUBROUTINE CDNOR(INA,INORM)
*     ***************************
*
*     THIS SUBROUTINE COMPUTES THE NORM OF THE CD VECTOR A
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
      ANORM = 0.D0
      DO 100 I=NBEG(INA),NEND(INA)
      ANORM = MAX(ANORM ,ABS(CC(I)))
 100  CONTINUE
*
      NTYP(INORM) = NRE
      CC(NBEG(INORM)) = ANORM
      RETURN
      END
*
      SUBROUTINE CDWERF(INA,INC)
*     **************************
*
*     THIS SUBROUTINE COMPUTES THE COMPLEX ERROR FUNCTION W OF COMPLEX DA INA
*
*     USE 7.1.20 PAGE 298 HANDBOOK MATH. FUNCTIONS.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XFR(0:LNO),XFI(0:LNO)
*
      IF(NTYP(INA).NE.NCD) CALL FOXNTY(INA)
*
      CALL FOXALL(ICR,1,2*NMMAX)
      CALL FOXALL(ICI,1,2*NMMAX)
      CALL FOXALL(ICD,1,2*NMMAX)
      CALL FOXALL(IUI,1,2)
      CALL FOXALL(IA0,1,2)
      CALL FOXALL(IW0,1,2)
      CALL FOXALL(IW1,1,2)
      CALL FOXALL(IW2,1,2)
      CALL FOXALL(ISR,1,1)
*
      NTYP(ISR) = NRE
      NBSR = NBEG(ISR)
*
      CALL IMUNIT(IUI)
      CALL CDCNST(INA,IA0)
*
      CALL CMWERF(IA0,IW0)
      XFR(0) = CC(NBEG(IW0)  )
      XFI(0) = CC(NBEG(IW0)+1)
*
      CALL CMMCM(IA0,IW0,IW1)
      CC(NBSR) = -2.D0
      CALL REMCM(ISR,IW1,IW1)
      CC(NBSR) = 1.D0/SQRT(ATAN(1.D0))
      CALL REMCM(ISR,IUI,IW2)
      CALL CMACM(IW1,IW2,IW1)
      XFR(1) = CC(NBEG(IW1)  )
      XFI(1) = CC(NBEG(IW1)+1)
*
      FACTO = 1.D0
      DO 10 I=2,NOCUT
      FACTO = FACTO*DBFLOAT(I)
      CALL CMMCM(IA0,IW1,IW2)
      CC(NBSR) = -2.D0
      CALL REMCM(ISR,IW2,IW2)
      CC(NBSR) = -2.D0*DBFLOAT(I-1)
      CALL REMCM(ISR,IW0,IW0)
      CALL CMACM(IW2,IW0,IW2)
      XFR(I) = CC(NBEG(IW2)  )/FACTO
      XFI(I) = CC(NBEG(IW2)+1)/FACTO
      CALL CMCOP(IW1,IW0)
      CALL CMCOP(IW2,IW1)
 10   CONTINUE
*
      CALL POLFUN(INA,1.D0,XFR,NOCUT,ICR)
      CALL POLFUN(INA,1.D0,XFI,NOCUT,ICI)
*
      CALL CDMCM(ICI,IUI,ICD)
      CALL CDACD(ICR,ICD,INC)
*
      CALL FOXDAL(ISR,1)
      CALL FOXDAL(IW2,1)
      CALL FOXDAL(IW1,1)
      CALL FOXDAL(IW0,1)
      CALL FOXDAL(IA0,1)
      CALL FOXDAL(IUI,1)
      CALL FOXDAL(ICD,1)
      CALL FOXDAL(ICI,1)
      CALL FOXDAL(ICR,1)
*
      PRINT*,'$$$ WARNING IN CDWERF, THIS ROUTINE IS UNDER DEVELOPMENT.'
      PRINT*,'    CONTACT US AT BERZ@MSU.EDU'
*
      RETURN
      END
*
      SUBROUTINE POLFUN(INA,CNORM,XF,NP,INC)
*     **************************************
*
*     THIS SUBROUTINE EVALUATES
*     C = XF(0) + XF(1)*AA + ...+ XF(NP)*AA^NP
*     WHERE AA IS THE NON-CONSTANT PART OF (A/CNORM)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----DIFFERENTIAL ALGEBRA --------------------------------------------------
      PARAMETER(LEA=100000,LIA=1400000,LNO=99,LNV=40)
      INTEGER IE1(LEA),IE2(LEA),IEO(LEA),IA1(0:LIA),IA2(0:LIA),
     *        NCFLT(LEA),IEW(LNV),IED(LNV),LEW,LEWI,IESP,
     *        NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT,ITM,LENTM,ITMPR
      DOUBLE PRECISION CDA(2*LEA),EPS,EPSMAC,TMT,TMS,EPSM,TOLTMR
      COMMON /DACOM/ CDA,EPS,EPSMAC,IE1,IE2,IEO,IA1,IA2,NCFLT,
     *       IEW,IED,LEW,LEWI,IESP,NOMAX,NVMAX,NMMAX,NOCUT,LFLT,NFLT
      COMMON /TMCOM/ TMT,TMS,EPSM,TOLTMR,ITM,LENTM,ITMPR
*----------------------------------------------------------------------------
*
      DIMENSION XF(0:LNO)
*
      IF(NTYP(INA).NE.NDA.AND.NTYP(INA).NE.NCD) CALL FOXNTY(INA)
*
      CALL FOXAAL(ISCO,1,2*MAX(NMMAX,2))
      CALL FOXAAL(ISINV,1,2*NMMAX)
      CALL FOXALL(ICA,1,2*NMMAX)
      CALL FOXALL(ISC,1,2*NMMAX)
      CALL FOXALL(IDA,1,2)
      CALL FOXALL(ISR,1,1)
*
      NTYP(ISR) = NRE
      NBSR = NBEG(ISR)
*
      CC(NBSR) = 1.D0
      CALL DA(ISR,IDA)
*
      IF(NP.EQ.0) THEN
         NTYP(ISCO+1) = NDA
         NBSC = NBEG(ISCO+1)
         IF(ABS(XF(NP)).GT.EPS) THEN
            CC(NBSC) = XF(NP)
            NC(NBSC) = 1
            NEND(ISCO+1) = NBSC
         ELSE
            NEND(ISCO+1) = NBSC-1
         ENDIF
      ELSE
         CC(NBSR) = XF(NP)
         CALL DAMRE(IDA,ISR,ISC)
*
         DO 10 I=NP-1,1,-1
         CC(NBSR) = XF(I)
         CALL DAARE(ISC,ISR,ISCO+1)
         CALL DAMDA(ISCO+1,IDA,ISC)
 10      CONTINUE
*
         CC(NBSR) = XF(0)
         CALL DAARE(ISC,ISR,ISCO+1)
      ENDIF
*
      CC(NBSR) = 1.D0/CNORM
      IF(NTYP(INA).EQ.NDA) THEN
         CALL DAMRE(INA,ISR,ICA)
         CALL DACNST(ICA,IDA)
         CALL DASRE(ICA,IDA,ISINV+1)
      ELSEIF(NTYP(INA).EQ.NCD) THEN
         CALL CDMRE(INA,ISR,ICA)
         CALL CDCNST(ICA,IDA)
         CALL CDSCM(ICA,IDA,ISINV+1)
      ENDIF
*
      CC(NBSR) = 1.D0
      CALL POLVAL(ISR,ISCO,ISR,ISINV,ISR,ISCO,ISR)
*
      IF(NTYP(INA).EQ.NDA) THEN
         CALL DACOP(ISCO+1,INC)
      ELSEIF(NTYP(INA).EQ.NCD) THEN
         CALL CDCOP(ISCO+1,INC)
      ENDIF
*
      CALL FOXDAL(ISR,1)
      CALL FOXDAL(IDA,1)
      CALL FOXDAL(ISC,1)
      CALL FOXDAL(ICA,1)
      CALL FOXADA(ISINV,1)
      CALL FOXADA(ISCO,1)
*
      RETURN
      END
*
      SUBROUTINE RKCO(IH,IA,IB,IC,ID)
*     *******************************
*
*     THIS SUBROUTINE FILLS THE RUNGA KUTTA COEFFICIENTS NEEDED BY FOXY'S
*     RUNGE KUTTA INTEGRATOR
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*
      INH = NBEG(IH)
      INA = NBEG(IA+1)-1
      INB = NBEG(IB+1)-14
      INC = NBEG(IC+1)-1
      IND = NBEG(ID+1)-1
*
      CC(INH)          =            1.D0 /           9.D0
*
      CC(INA+ 1)       =            0.D0
      CC(INA+ 2)       =            1.D0 /          18.D0
      CC(INA+ 3)       =            1.D0 /          12.D0
      CC(INA+ 4)       =            1.D0 /           8.D0
      CC(INA+ 5)       =            5.D0 /          16.D0
      CC(INA+ 6)       =            3.D0 /           8.D0
      CC(INA+ 7)       =           59.D0 /         400.D0
      CC(INA+ 8)       =           93.D0 /         200.D0
      CC(INA+ 9)       =   5490023248.D0 /  9719169821.D0
      CC(INA+10)       =           13.D0 /          20.D0
      CC(INA+11)       =   1201146811.D0 /  1299019798.D0
      CC(INA+12)       =            1.D0
      CC(INA+13)       =            1.D0
*
      CC(INB+ 2+13* 1) =            1.D0 /          18.D0
      CC(INB+ 3+13* 1) =            1.D0 /          48.D0
      CC(INB+ 3+13* 2) =            1.D0 /          16.D0
      CC(INB+ 4+13* 1) =            1.D0 /          32.D0
      CC(INB+ 4+13* 2) =            0.D0
      CC(INB+ 4+13* 3) =            3.D0 /          32.D0
      CC(INB+ 5+13* 1) =            5.D0 /          16.D0
      CC(INB+ 5+13* 2) =            0.D0
      CC(INB+ 5+13* 3) = -         75.D0 /          64.D0
      CC(INB+ 5+13* 4) =           75.D0 /          64.D0
      CC(INB+ 6+13* 1) =            3.D0 /          80.D0
      CC(INB+ 6+13* 2) =            0.D0
      CC(INB+ 6+13* 3) =            0.D0
      CC(INB+ 6+13* 4) =            3.D0 /          16.D0
      CC(INB+ 6+13* 5) =            3.D0 /          20.D0
      CC(INB+ 7+13* 1) =     29443841.D0 /   614563906.D0
      CC(INB+ 7+13* 2) =            0.D0
      CC(INB+ 7+13* 3) =            0.D0
      CC(INB+ 7+13* 4) =     77736538.D0 /   692538347.D0
      CC(INB+ 7+13* 5) = -   28693883.D0 /  1125000000.D0
      CC(INB+ 7+13* 6) =     23124283.D0 /  1800000000.D0
      CC(INB+ 8+13* 1) =     16016141.D0 /   946692911.D0
      CC(INB+ 8+13* 2) =            0.D0
      CC(INB+ 8+13* 3) =            0.D0
      CC(INB+ 8+13* 4) =     61564180.D0 /   158732637.D0
      CC(INB+ 8+13* 5) =     22789713.D0 /   633445777.D0
      CC(INB+ 8+13* 6) =    545815736.D0 /  2771057229.D0
      CC(INB+ 8+13* 7) = -  180193667.D0 /  1043307555.D0
      CC(INB+ 9+13* 1) =     39632708.D0 /   573591083.D0
      CC(INB+ 9+13* 2) =            0.D0
      CC(INB+ 9+13* 3) =            0.D0
      CC(INB+ 9+13* 4) = -  433636366.D0 /   683701615.D0
      CC(INB+ 9+13* 5) = -  421739975.D0 /  2616292301.D0
      CC(INB+ 9+13* 6) =    100302831.D0 /   723423059.D0
      CC(INB+ 9+13* 7) =    790204164.D0 /   839813087.D0
      CC(INB+ 9+13* 8) =    800635310.D0 /  3783071287.D0
      CC(INB+10+13* 1) =    246121993.D0 /  1340847787.D0
      CC(INB+10+13* 2) =     0.D0
      CC(INB+10+13* 3) =     0.D0
      CC(INB+10+13* 4) = -37695042795.D0 / 15268766246.D0
      CC(INB+10+13* 5) = -  309121744.D0 /  1061227803.D0
      CC(INB+10+13* 6) = -   12992083.D0 /   490766935.D0
      CC(INB+10+13* 7) =   6005943493.D0 /  2108947869.D0
      CC(INB+10+13* 8) =    393006217.D0 /  1396673457.D0
      CC(INB+10+13* 9) =    123872331.D0 /  1001029789.D0
      CC(INB+11+13* 1) = - 1028468189.D0 /   846180014.D0
      CC(INB+11+13* 2) =     0.D0
      CC(INB+11+13* 3) =     0.D0
      CC(INB+11+13* 4) =   8478235783.D0 /   508512852.D0
      CC(INB+11+13* 5) =   1311729495.D0 /  1432422823.D0
      CC(INB+11+13* 6) = -10304129995.D0 /  1701304382.D0
      CC(INB+11+13* 7) = -48777925059.D0 /  3047939560.D0
      CC(INB+11+13* 8) =  15336726248.D0 /  1032824649.D0
      CC(INB+11+13* 9) = -45442868181.D0 /  3398467696.D0
      CC(INB+11+13*10) =   3065993473.D0 /   597172653.D0
      CC(INB+12+13* 1) =    185892177.D0 /   718116043.D0
      CC(INB+12+13* 2) =     0.D0
      CC(INB+12+13* 3) =     0.D0
      CC(INB+12+13* 4) = - 3185094517.D0 /   667107341.D0
      CC(INB+12+13* 5) = -  477755414.D0 /  1098053517.D0
      CC(INB+12+13* 6) = -  703635378.D0 /   230739211.D0
      CC(INB+12+13* 7) =   5731566787.D0 /  1027545527.D0
      CC(INB+12+13* 8) =   5232866602.D0 /   850066563.D0
      CC(INB+12+13* 9) = - 4093664535.D0 /   808688257.D0
      CC(INB+12+13*10) =   3962137247.D0 /  1805957418.D0
      CC(INB+12+13*11) =     65686358.D0 /   487910083.D0
      CC(INB+13+13* 1) =    403863854.D0 /   491063109.D0
      CC(INB+13+13* 2) =            0.D0
      CC(INB+13+13* 3) =            0.D0
      CC(INB+13+13* 4) = - 5068492393.D0 /   434740067.D0
      CC(INB+13+13* 5) = -  411421997.D0 /   543043805.D0
      CC(INB+13+13* 6) =    652783627.D0 /   914296604.D0
      CC(INB+13+13* 7) =  11173962825.D0 /   925320556.D0
      CC(INB+13+13* 8) = -13158990841.D0 /  6184727034.D0
      CC(INB+13+13* 9) =   3936647629.D0 /  1978049680.D0
      CC(INB+13+13*10) = -  160528059.D0 /   685178525.D0
      CC(INB+13+13*11) =    248638103.D0 /  1413531060.D0
      CC(INB+13+13*12) =            0.D0
*
      CC(INC+ 1)       =     14005451.D0 /   335480064.D0
      CC(INC+ 2)       =            0.D0
      CC(INC+ 3)       =            0.D0
      CC(INC+ 4)       =            0.D0
      CC(INC+ 5)       =            0.D0
      CC(INC+ 6)       = -   59238493.D0 /  1068277825.D0
      CC(INC+ 7)       =    181606767.D0 /   758867731.D0
      CC(INC+ 8)       =    561292985.D0 /   797845732.D0
      CC(INC+ 9)       = - 1041891430.D0 /  1371343529.D0
      CC(INC+10)       =    760417239.D0 /  1151165299.D0
      CC(INC+11)       =    118820643.D0 /   751138087.D0
      CC(INC+12)       = -  528747749.D0 /  2220607170.D0
      CC(INC+13)       =            1.D0 /           4.D0
*
      CC(IND+ 1)       =     13451932.D0 /   455176623.D0
      CC(IND+ 2)       =            0.D0
      CC(IND+ 3)       =            0.D0
      CC(IND+ 4)       =            0.D0
      CC(IND+ 5)       =            0.D0
      CC(IND+ 6)       = -  808719846.D0 /   976000145.D0
      CC(IND+ 7)       =   1757004468.D0 /  5645159321.D0
      CC(IND+ 8)       =    656045339.D0 /   265891186.D0
      CC(IND+ 9)       = - 3867574721.D0 /  1518517206.D0
      CC(IND+10)       =    465885868.D0 /   322736535.D0
      CC(IND+11)       =     53011238.D0 /   667516719.D0
      CC(IND+12)       =            2.D0 /          45.D0
      CC(IND+13)       =            0.D0 /           1.D0
*
      RETURN
      END
*
