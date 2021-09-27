*******************************************************************************
*                                                                             *
*                                                                             *
*                                PROGRAM FOXY                                 *
*                        COMPILER AND EXECUTER PACKAGE                        *
*                                                                             *
*                           PART OF THE COSY SYSTEM                           *
*                                                                             *
*                                 VERSION 10.0                                *
*                                                                             *
*                           UPDATED IN MARCH, 2017                            *
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
*FACE                                                                    *FACE*
*                                                                             *
*******************************************************************************
*
*FACE SUBROUTINE ENTRY                                                   *FACE
*     ****************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----CODE ------------------------------------------------------------------
      PARAMETER(LCOD=2000000,LSYM=700000,LCOM=353)
      PARAMETER(LCON=20000,LLEV=100,LPRO=500,LSCR=50,NSCR=50000)
      INTEGER NCOD(LCOD),ICOD,ICOO, NSYM(LSYM,4),ISYM, ICON
      INTEGER NLEV(LLEV,2),ILEV, NPRO(LPRO,3),IPRO,JPRO
      INTEGER LERR,IBEG,MSCR,ILIS,IREF,ICOM,ICOP,IEXE
      DOUBLE PRECISION DCON(LCON)
      CHARACTER CSYM(LSYM)*32,COMM(LCOM)*8,CCON(LCON)*512,CREF*800
      COMMON /CODE/ NCOD,ICOD,ICOO,  DCON,ICON,  NSYM,ISYM,
     *              NLEV,ILEV,  NPRO,IPRO,JPRO,
     *              LERR,IBEG,MSCR,ILIS,IREF,ICOM,ICOP,IEXE
      COMMON /CODC/ CSYM, COMM, CCON, CREF
*----------------------------------------------------------------------------
*
      PARAMETER(LTYP=11,LSYNR=30,LOPS=20,LFUN=100,LSUB=200)
      CHARACTER CTYID(LTYP)*2
      COMMON /TYCID/ CTYID
*
*GEN  HEADER INSERT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
*     FOXY COMPILER V.10.0 GENERATED 31-Mar-2017 14:39:45
*
*     ALLOWED TYPES
*     *************
*     RE
*     ST
*     LO
*     CM
*     VE
*     DA
*     CD
*     GR
*
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
      NRE =   1
      NST =   2
      NLO =   3
      NCM =   4
      NVE =   5
      NDA =   6
      NCD =   7
      NGR =   8
*
      CTYID(  1) = 'RE'
      CTYID(  2) = 'ST'
      CTYID(  3) = 'LO'
      CTYID(  4) = 'CM'
      CTYID(  5) = 'VE'
      CTYID(  6) = 'DA'
      CTYID(  7) = 'CD'
      CTYID(  8) = 'GR'
*
*GEN  HEADER INSERT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
*     COMPILE: C, PRECOMPILE: P, EDIT INPUT: I, EDIT OUTPUT: O
*
*      CALL HEADER(6,' C O S Y    I N F I N I T Y  ',                     *NORM
*     *              '    FOXY LANGUAGE SYSTEM     ',                     *NORM
*     *              '        VERSION  10.0        ',                     *NORM
*     *              '        (C) MSU 2017         ')                     *NORM
*MPI  CALL MPI_INIT(MPERR)                                               *MPI
*
      CALL INPUT                                                         *NORM
*MPI  CALL INPUT                                                         *MPI
      CALL CODCOM                                                        *NORM
*MPI  CALL CODCOM                                                        *MPI
      CALL CODEXE(NCOD, DCON, CCON, LERR)                                *NORM
*MPI  CALL CODEXE(NCOD, DCON, CCON, LERR)                                *MPI
      CALL FOXSTP(0)                                                     *NORM
*MPI  CALL FOXSTP(0)                                                     *MPI
*
      END
*
      SUBROUTINE CODCOM
*     *****************
*
*     THIS SUBROUTINE COMPILES THE FOXY LANGUAGE TO METACODE. THE SYNTAX
*     DESCRIBED IN THE ARRAY CSYN IS USED.
*
*-----CODE ------------------------------------------------------------------
      PARAMETER(LCOD=2000000,LSYM=700000,LCOM=353)
      PARAMETER(LCON=20000,LLEV=100,LPRO=500,LSCR=50,NSCR=50000)
      INTEGER NCOD(LCOD),ICOD,ICOO, NSYM(LSYM,4),ISYM, ICON
      INTEGER NLEV(LLEV,2),ILEV, NPRO(LPRO,3),IPRO,JPRO
      INTEGER LERR,IBEG,MSCR,ILIS,IREF,ICOM,ICOP,IEXE
      DOUBLE PRECISION DCON(LCON)
      CHARACTER CSYM(LSYM)*32,COMM(LCOM)*8,CCON(LCON)*512,CREF*800
      COMMON /CODE/ NCOD,ICOD,ICOO,  DCON,ICON,  NSYM,ISYM,
     *              NLEV,ILEV,  NPRO,IPRO,JPRO,
     *              LERR,IBEG,MSCR,ILIS,IREF,ICOM,ICOP,IEXE
      COMMON /CODC/ CSYM, COMM, CCON, CREF
*----------------------------------------------------------------------------
*
*     MAKE SURE LAIN IS LARGE ENOUGH LIKE ABOUT 2000.
*     USE LAIN=4*512 WHEN CCON LENGTH IS 512. USE LAIN=5*256 IF 256.
      PARAMETER(LSYN=27,LNA=200,LLIS=1,LAIN=4*512)
      PARAMETER(LTYP=11,LSYNR=30,LOPS=20,LFUN=100,LSUB=200)
*
      INTEGER NA(LNA),NC(LNA)
      CHARACTER CTYP*1,CINS*3,A*(LAIN),CBLA*40,CID*3,CSEA*32
      CHARACTER CSYN(LSYN)*36,PRO*1200
      DATA LCHECK, LMAP / 0, 0 /
*
*     NCOD : CONTAINS THE META CODE CONSISTING OF VARIABLE-LENGTH RECORDS.
*            LAST ENTRY IS ICOD, THE ONE BEFORE LAST IS ICOO.
*            CODING IS AS FOLLOWS:
*     IC+1        : POINTER TO CODE IDENTIFIER (CF CODEXE)
*     IC+2        : POINTER TO NEXT COMMAND TO EXECUTE
*     IC+3        : POINTER TO ALTERNATE NEXT COMMAND
*     IC+4        : POINTER TO NEXT COMMAND IN NCOD (IN)
*     IC+5        : POINTER TO BEGINNING OF VARIABLE USE SECTION (IV)
*     IC+6 ... IV : ADDRESSES OF ARGUMENTS TO COMMAND
*     IV+1 ... IN : VARIABLE USE SECTION
*
*     CSYM: NAME OF THE SYMBOL
*
*     NSYM: ISYM PARAMETERS DESCRIBING THE SYMBOLS
*     1:    FOR VAR: 0, FOR FUN: START ADDRESS, FOR SUB: -START ADDRESS
*     2:    IPRO LEVEL AT WHICH IT OCCURED
*     3:    FOR VAR AND FUN: LOCAL VARIABLE ADDRESS, FOR SUB: 0
*     4:    ADDRESS OF CODE LIST LINE WHERE DEFINED
*
*     DCON: ICON DOUBLE PRECISION CONSTANTS
*     CCON: ICON CHARACTER CONSTANTS
*
*     NLEV: ILEV CURRENT RECURSION LEVELS.
*     1:    ID NUMBER OF LEVEL (LAST ENTRY IN CSYN)
*     2:    POINTER TO NCOD WHERE LEVEL BEGINS
*
*     NPRO:
*     1:    POINTER TO LAST VARIABLE SYMBOL BEFORE ROUTINE
*     2:    POINTER TO LAST D V STATEMENT OF ROUTINE
*     3:    MSCR BEFORE ROUTINE IS CALLED
*
*     LERR: ERROR FLAG. 0 IF CODE IS ERROR FREE, 1 ELSE
*     IBEG: FIRST STATEMENT IN CODE TO BE EXECUTED
*
*============================================================
      DATA CSYN( 1) / '$$$$$$$$$$      : =                 '/
      DATA CSYN( 2) / 'BEGIN           B E ;               '/
      DATA CSYN( 3) / 'END             E E ;               '/
      DATA CSYN( 4) / 'VARIABLE        D V D E - ;         '/
      DATA CSYN( 5) / 'FUNCTION        B F D -             '/
      DATA CSYN( 6) / 'ENDFUNCTION     E F ;               '/
      DATA CSYN( 7) / 'PROCEDURE       B P D - ;           '/
      DATA CSYN( 8) / 'ENDPROCEDURE    E P ;               '/
      DATA CSYN( 9) / '$$$$$$$$$$      C P E - ;           '/
      DATA CSYN(10) / '$$$$$$$$$$      C S E - ;           '/
*     C P: CALL USER PROCEDURE, C S: CALL INTRINSIC PROCEDURE
*
      DATA CSYN(11) / 'IF              B I E ;             '/
      DATA CSYN(12) / 'ELSEIF          O I E ;             '/
      DATA CSYN(13) / 'ENDIF           E I ;               '/
      DATA CSYN(14) / 'WHILE           B W E ;             '/
      DATA CSYN(15) / 'ENDWHILE        E W ;               '/
      DATA CSYN(16) / 'LOOP            B L N E E - ;       '/
      DATA CSYN(17) / 'ENDLOOP         E L ;               '/
      DATA CSYN(18) / 'FIT             B O N - ;           '/
      DATA CSYN(19) / 'ENDFIT          E O E E E E - ;     '/
      DATA CSYN(20) / 'DEBUG           D C ;               '/
*
      DATA CSYN(21) / 'WRITE           W F E - ;           '/
      DATA CSYN(22) / 'READ            R F E - ;           '/
      DATA CSYN(23) / 'SAVE            W C E ;             '/
      DATA CSYN(24) / 'INCLUDE         I C E ;             '/
      DATA CSYN(25) / 'PLOOP           B M N E E - ;       '/
      DATA CSYN(26) / 'ENDPLOOP        E M E N - ;         '/
      DATA CSYN(27) / 'GUIIO           G I E - ;           '/
*============================================================
*GEN  PROCEDURE INSERT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DATA PRO(   1:  30)  / 'MEMALLMEMFREMEMDPVMEMWRTSCRLEN' /
      DATA PRO(  31:  60)  / 'CPUSECPWTIMEPNPRO PROOT QUIT  ' /
      DATA PRO(  61:  90)  / 'SLEEPMOS    ARGGETOPENF OPENFB' /
      DATA PRO(  91: 120)  / 'CLOSEFREWF  BACKF READS READB ' /
      DATA PRO( 121: 150)  / 'WRITEBREADM WRITEMDAINI DANOT ' /
      DATA PRO( 151: 180)  / 'DANOTWDAEPS DAEPSMEPSMINDAFSET' /
      DATA PRO( 181: 210)  / 'DAFILTDAPEW DAREA DAPRV DAREV ' /
      DATA PRO( 211: 240)  / 'DAFLO CDFLO DAGMD RERAN DARAN ' /
      DATA PRO( 241: 270)  / 'DADIU DADMU DADER DAINT DAPLU ' /
      DATA PRO( 271: 300)  / 'DASCL DATRN DASGN DAPEE DAPEA ' /
      DATA PRO( 301: 330)  / 'DACODEDANORODANORSDACLIWDACQLC' /
      DATA PRO( 331: 360)  / 'DAPEP DANOW DAEST MTREE CDF2  ' /
      DATA PRO( 361: 390)  / 'CDNF  CDNFDACDNFDSLINV  LDET  ' /
      DATA PRO( 391: 420)  / 'LEV   MBLOCKLSLINESUBSTRSTCRE ' /
      DATA PRO( 421: 450)  / 'RECST VELSETVELGETVEDOT VEUNIT' /
      DATA PRO( 451: 480)  / 'VEZEROIMUNITLTRUE LFALSEINTPOL' /
      DATA PRO( 481: 510)  / 'CLEAR GRMOVEGRDRAWGRDOT GRTRI ' /
      DATA PRO( 511: 540)  / 'GRPOLYGRCURVGRCHARGRCOLRGRWDTH' /
      DATA PRO( 541: 570)  / 'GRPROJGRZOOMGRMIMAGREPS GRSTYL' /
      DATA PRO( 571: 600)  / 'GROUTFGUISETRKCO  POLSETPOLVAL' /
      DATA PRO( 601: 630)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 631: 660)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 661: 690)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 691: 720)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 721: 750)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 751: 780)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 781: 810)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 811: 840)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 841: 870)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 871: 900)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 901: 930)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 931: 960)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 961: 990)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO( 991:1020)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO(1021:1050)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO(1051:1080)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO(1081:1110)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO(1111:1140)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO(1141:1170)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA PRO(1171:1200)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
*GEN  PROCEDURE INSERT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DATA CINS / 'NDE' /
      DATA CBLA / '                                        ' /
*
*      WRITE(6,'(1X,A)') '--- BEGINNING COMPILATION '                     *NORM
*
      DO 10 J=1,LSYN
  10  COMM(J) = CSYN(J)(17:19)//'     '
      DO 15 J=LSYN+1,LSYNR
  15  COMM(J) = '$$$$$$$$'
      DO 20 J=1,LSUB
      IF(PRO(6*(J-1)+1:6*J).NE.'      ') THEN
         COMM(LSYNR+LOPS+LFUN+J) = PRO(6*(J-1)+1:6*J)//'  '
      ELSE
         COMM(LSYNR+LOPS+LFUN+J) = '$$$$$$$$'
      ENDIF
  20  CONTINUE
*
      ICOD = 0
      ISYM = 0
      ICON = 0
      ILEV = 0
      IPRO = 0
      MSCR = 0
      LA   = 0
      ILIS = 0
*
  100 CONTINUE
*
      CALL GETCOM(A,LAIN,LA, NA,LNA,IA, IASS,IEND,LERR)
*
      ILIS  = ILIS + 1
      IREF  = 0
      ILEVO = ILEV
      JO    = 0
  105 JO = JO + 1
      IF(A(JO:JO).EQ.' ') GOTO 105
      JOO   = JO
      JO    = MAX(JO,LA-78)
      JO    = MAX(1,JO)
*
      IF(IA.LT.2) THEN
         CALL MSG('### ERROR, NO COMMAND FOUND '//A(JO:LA))
         LERR = 1
         GOTO 500
*      ELSEIF(A(NA(IA):NA(IA)).NE.';') THEN
*         CALL MSG( '### ERROR, LAST ENTRY NOT ";" '//A(JO:LA))
*         LERR = 1
*         GOTO 500
      ENDIF
      IF(LCHECK.EQ.1) WRITE(6,*) ' COMMAND: ',
*    *                A(1:LA),(NA(J),J=1,MIN(3,IA))
     *                A(1:LA),(NA(J),J=1,IA)
      IF(IEND.EQ.1) GOTO 1000
*
*     IDENTIFYING COMMAND
*     *******************
*
      IL = ILAST(A,NA(1),NA(2)-1)
*
      IF(IASS.EQ.1) THEN
         IID = 1
         CID = 'ASS'
         NC1 = 1
         NC3 = 0
         GOTO 120
      ENDIF
*
      IF(IL-NA(1)+1.LE.15) THEN
         DO 110 IID=1,LSYN
         IF(A(NA(1):IL).NE.CSYN(IID)(1:IL-NA(1)+1)) GOTO 110
         IF(CSYN(IID)(IL-NA(1)+2:IL-NA(1)+2).NE.' ') GOTO 110
         CID = CSYN(IID)(17:19)
         NC1 = IID
         NC3 = 0
         GOTO 120
 110     CONTINUE
      ENDIF
*
      IF(IL-NA(1)+1.LE.32) THEN
         DO 115 IIID=ISYM,LCON+1,-1
         IF(NSYM(IIID,1).GE.0) GOTO 115
         IF(A(NA(1):IL).NE.CSYM(IIID)(1:IL-NA(1)+1)) GOTO 115
         IF(CSYM(IIID)(IL-NA(1)+2:IL-NA(1)+2).NE.' ') GOTO 115
         DO 113 IID=1,LSYN
         IF(CSYN(IID)(17:19).EQ.'C P') GOTO 114
 113     CONTINUE
 114     CONTINUE
         NC1 = IID
         NC3 = -NSYM(IIID,1)
         CID = 'C P'
         CREF(IREF+1:IREF+IL-NA(1)+4) = ' %'//CSYM(IIID)(1:IL-NA(1)+1)
     *                                  //':'
         IREF = IREF + IL - NA(1) + 4
         IF(NSYM(IIID,4).LT.100000) THEN
            WRITE(CREF(IREF+1:IREF+5),'(I5)') NSYM(IIID,4)
            IREF = IREF + 5
         ELSE
            WRITE(CREF(IREF+1:IREF+8),'(I8)') NSYM(IIID,4)
            IREF = IREF + 8
         ENDIF
         GOTO 120
 115     CONTINUE
      ENDIF
*
      IF(IL-NA(1)+1.LE.6) THEN
         CSEA(1:6) = '      '
         CSEA(1:IL-NA(1)+1) = A(NA(1):IL)
         II = (INDEX(PRO,CSEA(1:6))+5)/6
         IF(II.NE.0) THEN
            DO 117 IID=1,LSYN
            IF(CSYN(IID)(17:19).EQ.'C S') GOTO 118
 117        CONTINUE
 118        CONTINUE
            NC1 = LSYNR+LOPS+LFUN + II
            NC3 = 0
            CID = 'C S'
            GOTO 120
         ENDIF
      ENDIF
*
      CALL MSG( '### ERROR, UNKNOWN COMMAND: '//A(JO:LA))
      LERR = 1
      GOTO 500
*
 120  CONTINUE
*
*     CODE READ AND WRITE
*     *******************
*
      IF(ICOD.EQ.0) THEN
         IF(CID.NE.'I C'.AND.CID.NE.'B E') THEN
            CALL MSG( '### ERROR, ILLEGAL COMMAND AT BEGINNING')
            CALL FOXSTP(1)
         ENDIF
      ENDIF
*
      IF(CID.EQ.'W C') THEN
         IF(A(NA(2):NA(2)).NE.'''') THEN
            CALL MSG( '### ERROR IN SAVE, FILENAME NOT STRING')
            CALL FOXSTP(1)
         ENDIF
         INA3 = NA(3)
  121    INA3 = INA3 - 1
         IF(A(INA3:INA3).NE.'''') GOTO 121
         OPEN(11,FILE=A(NA(2)+1:INA3-1)//'.bin',STATUS='UNKNOWN',
     *           FORM='UNFORMATTED',ERR=101)
         WRITE (11) LCOD,LSYM,LCOM,LCON,LLEV,LPRO,LSCR,NSCR
         WRITE (11) NCOD,ICOD,ICOO
         WRITE (11) DCON,ICON
         WRITE (11) NSYM,ISYM
         WRITE (11) NLEV,ILEV
         WRITE (11) NPRO,IPRO,JPRO
         WRITE (11) LERR,IBEG,MSCR,ILIS,IREF,ICOM,ICOP,IEXE
         WRITE (11) CSYM, COMM, CCON, CREF
         CLOSE (11)
         CALL MSG( '--- BIN FILE WRITTEN: '//A(NA(2)+1:INA3-1))
         GOTO 500
  101    CONTINUE
         CALL MSG( '### ERROR WHILE OPENING FILE: '//A(JO:LA))
         LERR = 1
         GOTO 500
      ELSEIF(CID.EQ.'I C') THEN
         IF(A(NA(2):NA(2)).NE.'''') THEN
            CALL MSG( '### ERROR IN INCLUDE, FILENAME NOT STRING')
            CALL FOXSTP(1)
         ENDIF
         INA3 = NA(3)
  122    INA3 = INA3 - 1
         IF(A(INA3:INA3).NE.'''') GOTO 122
         OPEN(11,FILE=A(NA(2)+1:INA3-1)//'.bin',STATUS='UNKNOWN',
     *           FORM='UNFORMATTED',ERR=102)
         READ (11) LLCOD,LLSYM,LLCOM,LLCON,LLLEV,LLPRO,LLSCR,LNSCR
         IF((LLCOD.NE.LCOD).OR.(LLSYM.NE.LSYM).OR.(LLCOM.NE.LCOM).OR.
     *      (LLCON.NE.LCON).OR.(LLLEV.NE.LLEV).OR.(LLPRO.NE.LPRO).OR.
     *      (LLSCR.NE.LSCR).OR.(LNSCR.NE.NSCR)) THEN
            CALL MSG( '### ERROR, FILE INCOMPATIBLE, RECOMPILE SOURCE:'
     *              //A(JO:LA))
            LERR = 1
            GOTO 500
         ENDIF
         READ (11) NCOD,ICOD,ICOO
         READ (11) DCON,ICON
         READ (11) NSYM,ISYM
         READ (11) NLEV,ILEV
         READ (11) NPRO,IPRO,JPRO
         READ (11) LERR,IBEG,MSCR,ILIS,IREF,ICOM,ICOP,IEXE
         READ (11) CSYM, COMM, CCON, CREF
*         CALL MSG( '--- BIN FILE READ: '//A(NA(2)+1:INA3-1))
         CLOSE(11)
         GOTO 500
  102    CONTINUE
         CALL MSG( '### ERROR WHILE OPENING FILE: '//A(JO:LA))
         LERR = 1
         GOTO 500
      ENDIF
*
      IF(ICOD+IA*10.GT.LCOD) THEN
         CALL MSG( '!!! MEMORY EXHAUSTION, INCREASE LCOD')
         CALL FOXSTP(1)
      ENDIF
*
      ICOM = ICOD
      ICOP = ICOO
      JTYP = 11
*
*     ASSIGNMENT
*     **********
*
      IF(CID.EQ.'ASS') THEN
         IF(IEXE.NE.3) NCOD(NPRO(IPRO,2)+2) = ICOD
         IEXE = 3
         CALL ARICOM(A,1,LA,0)
         GOTO 500
      ENDIF
*
*     UPDATE NPRO, IPRO IF NECESSARY
*     ******************************
*
      IF(CID.EQ.'B E') THEN
         JPRO = 1
         ISYM = LCON + 1
         CSYM(ISYM)(1:1) = '$'
         IPRO = IPRO + 1
         NPRO(IPRO,1) = ISYM + 1
         NPRO(IPRO,2) = ICOD
         NPRO(IPRO,3) = 0
         DO 129 IIID = ISYM+1,ISYM+LSCR+1
  129    CSYM(IIID)(1:1) = '$'
         ISYM = ISYM + LSCR + 1
      ELSEIF(CID.EQ.'B F'.OR.CID.EQ.'B P') THEN
         IPRO = IPRO + 1
         NPRO(IPRO,1) = ISYM + 1
         NPRO(IPRO,2) = ICOD
         NPRO(IPRO,3) = MSCR
         DO 128 IIID = ISYM+1,ISYM+LSCR+1
  128    CSYM(IIID)(1:1) = '$'
         ISYM = ISYM + LSCR + 1
      ENDIF
*
*     CHECKING DECLARATION - EXECUTION PLACEMENT
*     ******************************************
*
      IF(CID.EQ.'D V') THEN
         IF(IEXE.NE.1) THEN
            CALL MSG( '### COMMAND PLACEMENT ERROR: '//A(JO:LA))
            LERR = 1
            GOTO 500
         ENDIF
      ELSEIF(CID.EQ.'B E'.OR.CID.EQ.'B P'.OR.CID.EQ.'B F') THEN
         IF(IEXE.EQ.3) THEN
            CALL MSG( '### COMMAND PLACEMENT ERROR: '//A(JO:LA))
            LERR = 1
            GOTO 500
         ENDIF
         IEXE = 1
      ELSEIF(CID.EQ.'E E'.OR.CID.EQ.'E P'.OR.CID.EQ.'E F') THEN
         IF(IEXE.NE.3) THEN
            CALL MSG( '### NO EXECUTABLE IN SEGMENT: '//A(JO:LA))
            LERR = 1
         ENDIF
         IEXE = 2
      ELSE
         IF(IEXE.NE.3) NCOD(NPRO(IPRO,2)+2) = ICOD
         IEXE = 3
      ENDIF
*
*     LOOP OVER ALL SUPPLIED PARAMETERS
*     *********************************
*
      IF(IA-2.GT.LSCR) THEN
         CALL MSG('### ERROR, TOO MANY ARGUMENTS: '//A(1:LA))
         WRITE(6,'(1X,A)') '### ERRORS IN CODE, NO EXECUTION'
         CALL FOXSTP(1)
      ENDIF
*
      IA1 = 2
*
      DO 150 JA=IA1,IA
*
      IC = JA - IA1 + 1
      CTYP = CSYN(IID)(2*JTYP-1:2*JTYP-1)
*
      IF(JA.EQ.IA) THEN
         IF(CTYP.NE.';'.AND.CTYP.NE.'-'.AND.CID.NE.'C P'
     *      .AND.CID.NE.'E M') THEN
            CALL MSG('### ERROR, TOO FEW PARAMETERS: '//A(JO:LA))
            LERR = 1
            GOTO 500
         ENDIF
         GOTO 150
      ELSEIF(CTYP.EQ.';') THEN
         CALL MSG('### ERROR, TOO MANY PARAMETERS: '//A(JO:LA))
         LERR = 1
         GOTO 500
      ELSEIF(CTYP.EQ.'-') THEN
         CONTINUE
      ELSE
         ITYP = INDEX(CINS,CTYP)
         JTYP = JTYP + 1
      ENDIF
*
      GOTO (130,135,140), ITYP
*
      CALL MSG('@@@ ERROR IN LANGUAGE SYNTAX WITH : '//CID)
      CALL FOXSTP(1)
*
*     NAME
*     ----
*
  130 CONTINUE
      IF(LCHECK.EQ.1) WRITE(6,*) '   N: ',A(NA(JA):NA(JA+1)-1)
      CSEA = CBLA(1:32)
      CSEA = A(NA(JA):NA(JA+1)-1)
      DO 132 J=ISYM,LCON+1,-1
      IF(CSYM(J).EQ.CSEA) GOTO 133
  132 CONTINUE
      CALL MSG( '### ERROR, UNKNOWN NAME '//A(JO:LA))
      LERR = 1
      GOTO 500
  133 CONTINUE
      NC(IC) = J
      GOTO 150
*
*     DECLARATION
*     -----------
*
  135 CONTINUE
      IF(LCHECK.EQ.1) WRITE(6,*) '   D: ',A(NA(JA):NA(JA+1)-1)
      IF(INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ',A(NA(JA):NA(JA))).EQ.0) THEN
         CALL MSG( '### ERROR IN SYNTAX OF NAME '//A(JO:LA))
         LERR = 1
         GOTO 500
      ENDIF
      DO 136 J=NA(JA),NA(JA+1)-1
      IF(INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_ ',
     *   A(J:J)).EQ.0) THEN
         CALL MSG( '### ERROR IN SYNTAX OF NAME '//A(JO:LA))
         LERR = 1
         GOTO 500
      ENDIF
 136  CONTINUE
      IF(CID.EQ.'D V'.OR.
     *  (CID.EQ.'B F'.AND.JA.NE.2).OR.(CID.EQ.'B P'.AND.JA.NE.2)) THEN
         ISYM = ISYM + 1
         IF(ISYM.GT.LSYM) THEN
            CALL MSG( '!!! MEMORY EXHAUSTION, INCREASE LSYM')
            CALL FOXSTP(1)
         ENDIF
         CSYM(ISYM) = CBLA(1:32)
         CSYM(ISYM)(1:NA(JA+1)-NA(JA)) = A(NA(JA):NA(JA+1)-1)
         NSYM(ISYM,1) = 0
         NSYM(ISYM,2) = IPRO
         NSYM(ISYM,3) = ISYM
         NSYM(ISYM,4) = ILIS
      ELSEIF(CID.EQ.'B F') THEN
         CSYM(NPRO(IPRO,1)) = CBLA(1:32)
         CSYM(NPRO(IPRO,1))(1:NA(JA+1)-NA(JA)) = A(NA(JA):NA(JA+1)-1)
         NSYM(NPRO(IPRO,1),1) = 0
         NSYM(NPRO(IPRO,1),2) = IPRO
         NSYM(NPRO(IPRO,1),3) = NPRO(IPRO,1)
         NSYM(NPRO(IPRO,1),4) = ILIS
      ELSEIF(CID.EQ.'B P') THEN
         CSYM(NPRO(IPRO,1)) = CBLA(1:32)
         CSYM(NPRO(IPRO,1))(1:NA(JA+1)-NA(JA)) = A(NA(JA):NA(JA+1)-1)
         NSYM(NPRO(IPRO,1),1) = -ICOD
         NSYM(NPRO(IPRO,1),2) = IPRO
         NSYM(NPRO(IPRO,1),3) = NPRO(IPRO,1)
         NSYM(NPRO(IPRO,1),4) = ILIS
      ELSE
         CALL MSG( '@@@ ERROR IN LANGUAGE SYNTAX WITH D')
         CALL FOXSTP(1)
      ENDIF
      NC(IC) = ISYM
      GOTO 150
*
*     EXPRESSION
*     ----------
*
  140 CONTINUE
      IF(LCHECK.EQ.1) WRITE(6,*) '   E: ',A(NA(JA):NA(JA+1)-1)
      IL = ILAST(A,NA(JA),NA(JA+1)-1)
      IF(IL-NA(JA).GE.32) GOTO 142
      DO 141 IIID=ISYM,LCON+1,-1
      IF(CSYM(IIID)(1:1).EQ.'$') GOTO 141
      IF(A(NA(JA):IL).NE.CSYM(IIID)(1:IL-NA(JA)+1)) GOTO 141
      IF(CSYM(IIID)(IL-NA(JA)+2:IL-NA(JA)+2).NE.' ') GOTO 141
      NC(IC) = IIID
      GOTO 150
  141 CONTINUE
  142 CONTINUE
*
      IF(A(NA(JA+1)-1:NA(JA+1)-1).EQ.' ')
     *   A(NA(JA+1)-1:NA(JA+1)-1) = ';'
      MSCR = MAX(MSCR,JA-IA1+1)
      CALL ARICOM(A,NA(JA),NA(JA+1)-1,JA-IA1+1+NPRO(IPRO,1))
      IF(A(NA(JA+1)-1:NA(JA+1)-1).EQ.';')
     *   A(NA(JA+1)-1:NA(JA+1)-1) = ' '
      NC(IC) = NPRO(IPRO,1) + JA - IA1 + 1
      GOTO 150
*
  150 CONTINUE
*
      IF(CID.EQ.'D V') NPRO(IPRO,2) = ICOD
*
*     MAIN PART OF COMMAND
*     ********************
*
      IP = ICOD + 5
      IV = ICOD + 3 + IA
*
      IF(CID.EQ.'B E') THEN
         IN = IV + 6
      ELSEIF(CID.EQ.'B F') THEN
         IN = IV + 4
      ELSEIF(CID.EQ.'B P') THEN
         IN = IV + 4
      ELSEIF(CID.EQ.'O I') THEN
         IN = IV + 2
      ELSEIF(CID.EQ.'B W') THEN
         IN = IV + 1
      ELSE
         IN = IV
      ENDIF
*
      NCOD(ICOD+1) = NC1
      NCOD(ICOD+2) = IN
      NCOD(ICOD+3) = NC3
      NCOD(ICOD+4) = IN
      IF(IV.EQ.IN) THEN
         NCOD(ICOD+5) = -ILIS
      ELSE
         NCOD(ICOD+5) = IV
      ENDIF
*
*     ANALYZE CORRECTNESS OF NESTING STRUCTURE
*     ****************************************
*
      JLEV = ICHAR(CID(3:3))
*
      IF(CID(1:1).EQ.'B') THEN
         ILEV = ILEV + 1
         IF(ILEV.GT.LLEV) THEN
            CALL MSG( '!!! MEMORY EXHAUSTION, INCREASE LLEV')
            CALL FOXSTP(1)
         ENDIF
         NLEV(ILEV,1) = JLEV
         NLEV(ILEV,2) = ICOD
         IF(CID.EQ.'B W') NCOD(IV+1) = ICOM
         IF(CID.EQ.'B E'.AND.ICOD.NE.0) THEN
            CALL MSG('### COMMAND PLACEMENT ERROR: '//A(JO:LA))
            LERR = 1
            GOTO 500
         ENDIF
      ELSEIF(CID(1:1).EQ.'E') THEN
         IF(ILEV.EQ.0) THEN
            CALL MSG('### ERROR, NO NESTING OPEN: '//A(JO:LA))
            LERR = 1
            GOTO 500
         ELSEIF(NLEV(ILEV,1).NE.JLEV) THEN
            DO 153 JSYN=1,LSYN
               IF('E '//CHAR(NLEV(ILEV,1)).EQ.CSYN(JSYN)(17:19))
     *         CALL MSG('### ERROR, ILLEGAL NESTING, EXPECTED '
     *            //CSYN(JSYN)(1:16))
  153       CONTINUE
            LERR = 1
            GOTO 500
         ELSE
            NCOD(ICOD+3) = NLEV(ILEV,2)
            NCOD(NLEV(ILEV,2)+3) = ICOD
            IF(CID.EQ.'E E') THEN
               IVEE = NCOD(NLEV(ILEV,2)+5)
               NCOD(IVEE+1) = MSCR
               NCOD(IVEE+2) = LSCR
               NCOD(IVEE+3) = NSCR
               NCOD(IVEE+4) = ICON
               NCOD(IVEE+5) = LCON
               NCOD(IVEE+6) = JPRO
               ILEV = ILEV - 1
               IF(ILEV.NE.0) THEN
                  CALL MSG( '### ERROR, NESTING B '//
     *                     CHAR(NLEV(ILEV,1))//' STILL OPEN')
                  LERR = 1
                  GOTO 500
               ENDIF
               DO 155 JSYM=NPRO(IPRO,1)+1,ISYM
 155           CSYM(JSYM)(1:1) = '$'
*              ISYM = NPRO(IPRO,1)
               MSCR = NPRO(IPRO,3)
               IPRO = IPRO - 1
               IF(IPRO.NE.0) THEN
                  CALL MSG( '### ERROR, PROCEDURES STILL OPEN')
                  CALL FOXSTP(1)
               ENDIF
               GOTO 1000
            ELSEIF(CID.EQ.'E F'.OR.CID.EQ.'E P') THEN
               JPRO = JPRO + 1
               IVEE = NCOD(NLEV(ILEV,2)+5)
               NCOD(IVEE+1) = MSCR
               NCOD(IVEE+2) = JPRO
               NCOD(IVEE+3) = NPRO(IPRO,1)
               NCOD(IVEE+4) = ISYM - NPRO(IPRO,1)
               MSCR = 0
               DO 156 JSYM=NPRO(IPRO,1)+1,ISYM
 156           CSYM(JSYM)(1:1) = '$'
               MSCR = NPRO(IPRO,3)
               IF(CID.EQ.'E F') NSYM(NPRO(IPRO,1),1) = NLEV(ILEV,2)
               IPRO = IPRO - 1
            ELSEIF(CID.EQ.'E I') THEN
               NCOD(NLEV(ILEV,2)+3) = ICOD
               ICEI = NCOD(NLEV(ILEV,2)+1)
               IVEI = NCOD(NLEV(ILEV,2)+5)
  177          CONTINUE
               IF(CSYN(ICEI)(17:19).EQ.'B I') GOTO 178
               NCOD(NCOD(IVEI+2)+2) = ICOD
               ICEI = NCOD(NCOD(IVEI+1)+1)
               IVEI = NCOD(NCOD(IVEI+1)+5)
               GOTO 177
  178          CONTINUE
               NCOD(ICOD+3) = IVEI
            ELSEIF(CID.EQ.'E W') THEN
               NCOD(ICOD+3) = NCOD(NLEV(ILEV,2)+7)
            ENDIF
            ILEV = ILEV - 1
         ENDIF
      ELSEIF(CID(1:1).EQ.'O') THEN
         IF(JLEV.NE.NLEV(ILEV,1)) THEN
         CALL MSG('### ERROR, ILLEGAL NESTING, NO B '//
     *                    CHAR(JLEV)//' FOUND FOR '//A(JO:LA))
            LERR = 1
            GOTO 500
         ENDIF
         NCOD(NLEV(ILEV,2)+3) = ICOM
         NCOD(IV+1) = NLEV(ILEV,2)
         NCOD(IV+2) = ICOP
         NLEV(ILEV,2) = ICOD
         IN = IV + 2
      ENDIF
*
      DO 190 J=1,IA-IA1
  190 NCOD(IP+J) = NC(J)
*
      ICOO = ICOD
      ICOD = NCOD(ICOD+4)
*
      GOTO 500
*
*     FINISH INPUT OF LINE
*     ********************
*
  500 CONTINUE
*
      IF(LLIS.EQ.1) THEN
         WRITE(2,'(1X,I6,1X,70A/18X,100(60A/18X))')
     *         ILIS,(' ',' ',J=1,ILEVO),(A(J:J),J=JOO,LA),
     *         (CREF(J:J),J=1,IREF)
      ENDIF
*
      GOTO 100
*
*     END OF INPUT
*     ************
*
 1000 CONTINUE
*
      IF(CID.NE.'E E'.AND.CID.NE.'W C') THEN
         CALL MSG( '### ERROR, ILLEGAL END: '//A(JO:LA))
         LERR = 1
      ENDIF
*
      IF(LMAP.EQ.1) THEN
         OPEN(3,FILE='COSYdebug.cod',STATUS='UNKNOWN')
         IC   = 0
         IMAP = 0
 1010    CONTINUE
         IF(IC.GE.ICOD) GOTO 1020
         IMAP = IMAP + 1
         J1   = IC + 2
         J2   = NCOD(IC+4)
         WRITE(3,'(1X,2I6, 2X,A,4I6,5X,4I6:10(/1X,22X,9I6))')
     *      IMAP,IC,COMM(NCOD(IC+1)),(NCOD(J),J=J1,J2)
         IC   = J2
         GOTO 1010
 1020    CONTINUE
         CLOSE(3)
      ENDIF
*
      CLOSE(1)
      CLOSE(2)
*
      RETURN
      END
*
      SUBROUTINE ARICOM(A, IA1,IA2, IAST)
*     ***********************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER(LANA=1000,LDEC=200,LNAM=32)
      PARAMETER(LTYP=11,LSYNR=30,LOPS=20,LFUN=100,LSUB=200)
*
*-----CODE ------------------------------------------------------------------
      PARAMETER(LCOD=2000000,LSYM=700000,LCOM=353)
      PARAMETER(LCON=20000,LLEV=100,LPRO=500,LSCR=50,NSCR=50000)
      INTEGER NCOD(LCOD),ICOD,ICOO, NSYM(LSYM,4),ISYM, ICON
      INTEGER NLEV(LLEV,2),ILEV, NPRO(LPRO,3),IPRO,JPRO
      INTEGER LERR,IBEG,MSCR,ILIS,IREF,ICOM,ICOP,IEXE
      DOUBLE PRECISION DCON(LCON)
      CHARACTER CSYM(LSYM)*32,COMM(LCOM)*8,CCON(LCON)*512,CREF*800
      COMMON /CODE/ NCOD,ICOD,ICOO,  DCON,ICON,  NSYM,ISYM,
     *              NLEV,ILEV,  NPRO,IPRO,JPRO,
     *              LERR,IBEG,MSCR,ILIS,IREF,ICOM,ICOP,IEXE
      COMMON /CODC/ CSYM, COMM, CCON, CREF
*----------------------------------------------------------------------------
*
*     THIS SUBROUTINE DECODES THE CODE STORED ON CHARACTER A FROM IA1 TO IA2
*     INTO ELEMENTARY OPERATIONS AND STORES THE RESULT IN NCOD.
*
*     IF IAST = 0, THE CODE IS A FULL ASSIGNMENT INCLUDING EQUAL SIGN;
*     IF IAST > 0, THE CODE IS AN ARITHMETIC EXPRESSION WITHOUT EQUAL SIGN, AND
*                  THE RESULT IS STORED IN ADDRESS IAST
*
*     THE CODING IS AS FOLLOWS:
*                  1    L (CODE FOR OPERATION):
*                       L BETWEEN  31 AND  50: INTRINSIC OPERATOR
*                       L BETWEEN  51 AND 150: INTRINSIC FUNCTIONS, WHERE
*                       L =  51:               COPY
*                       L = 351:               ARRAY LOOKUP
*                       L = 352:               COPY INTO ARRAY
*                       L = 353:               EXTERNAL FUNCTION
*                  2    ADDRESS OF NEXT OPERATION
*                  3,4  ADDRESS, PROCEDURE LEVEL OF RESULT
*                  5,6  ADDRESS, PROCEDURE LEVEL OF FIRST OPERAND
*                  7,8, ADDRESS, PROCEDURE LEVEL OF FURTHER OPERANDS
*
*                  CONSTANTS ARE STORED AT PROCEDURE LEVEL 0.
*
      INTEGER NANA(LANA,6),NPRI(20)
      CHARACTER A*(*)
      CHARACTER OPS*20,FUN*600,CANA*10
      CHARACTER BLANK*32,NUM*10,LET*26,CSEA*32,AER*50,CVAL*(LDEC)
      SAVE ICALL, MPRI
*
      DATA LCHECK / 0 /
*GEN  OPERATOR INSERT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DATA OPS  / '+-*/^<>=#&|%        ' /
      DATA NPRI / 3,3,4,4,5,2,2,2,2,2,6,7,0,0,0,0,0,0,0,0 /
*GEN  OPERATOR INSERT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*GEN  FUNCTION INSERT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DATA FUN(   1:  30)  / 'COPY  RE    ST    LO    CM    ' /
      DATA FUN(  31:  60)  / 'VE    DA    CD    LRE   LST   ' /
      DATA FUN(  61:  90)  / 'LLO   LCM   LVE   LDA   LCD   ' /
      DATA FUN(  91: 120)  / 'LGR   TYPE  LENGTHVARMEMVARPOI' /
      DATA FUN( 121: 150)  / 'EXP   LOG   SIN   COS   TAN   ' /
      DATA FUN( 151: 180)  / 'ASIN  ACOS  ATAN  SINH  COSH  ' /
      DATA FUN( 181: 210)  / 'TANH  SQRT  ISRT  ISRT3 SQR   ' /
      DATA FUN( 211: 240)  / 'ERF   WERF  VMIN  VMAX  ABS   ' /
      DATA FUN( 241: 270)  / 'NORM  CONS  REAL  IMAG  CMPLX ' /
      DATA FUN( 271: 300)  / 'CONJ  INT   NINT  NOT   TRIM  ' /
      DATA FUN( 301: 330)  / 'LTRIM GRIU  $$$$$$$$$$$$$$$$$$' /
      DATA FUN( 331: 360)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA FUN( 361: 390)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA FUN( 391: 420)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA FUN( 421: 450)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA FUN( 451: 480)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA FUN( 481: 510)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA FUN( 511: 540)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA FUN( 541: 570)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
      DATA FUN( 571: 600)  / '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$' /
*GEN  FUNCTION INSERT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
      DATA CANA / '=OFAVSC#,;' /
      DATA NUM  / '1234567890' /
      DATA LET  / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
*
      DATA AER  / '                                                  ' /
      DATA ICALL  / 0 /
      DATA BLANK / '                                ' /
*
      IF(ICALL.EQ.0) THEN
         ICALL = 1
         MPRI = 1
         DO 11 J=1,LOPS
  11     MPRI = MAX(MPRI,NPRI(J))
         MPRI = MPRI + 2
*
         DO 16 J=1,LOPS
         IF(OPS(J:J).NE.' ') THEN
            COMM(LSYNR+J) = OPS(J:J)//'       '
         ELSE
            COMM(LSYNR+J) = '$$$$$$$$'
         ENDIF
  16     CONTINUE
         DO 18 J=1,LFUN
         IF(FUN(6*(J-1)+1:6*J).NE.'      ') THEN
            COMM(LSYNR+LOPS+J) = FUN(6*(J-1)+1:6*J)//'  '
         ELSE
            COMM(LSYNR+LOPS+J) = '$$$$$$$$'
         ENDIF
  18     CONTINUE
         COMM( 51) = 'COPY    '
*
         COMM(351) = 'ARRAY   '
         COMM(352) = 'COPINARR'
         COMM(353) = 'FUNCTION'
      ENDIF
*
*-----------------------------------------------------------------------
*
      IF(LCHECK.EQ.1) WRITE(6,'(1X,A)') A(IA1:IA2)
*
*     SETTING START VALUES
*     ********************
*
      INCR   = 0
      INCRO  = 0
      ILEFT  = IAST
      ILEFTO = IPRO
      ISCR   = 0
      IF(IAST.NE.0) ISCR = IAST - NPRO(IPRO,1)
      IOP    = 1
      IPAR   = 0
      IEQU   = 0
      IFUVA  = 0
      I      = IA1 - 1
*
      IANA  = 1
      NANA(IANA,1) = INDEX(CANA,';')
      NANA(IANA,2) = 0
      NANA(IANA,3) = 0
      NANA(IANA,4) = 0
      NANA(IANA,5) = 0
      NANA(IANA,6) = IANA + 1
*
*     PROCESSING LEFT SIDE OF EQUAL SIGN
*     **********************************
*
      IF(ILEFT.NE.0) THEN
         I = IA1
         IEQU = 1
         IANA  = IANA + 1
         NANA(IANA,1) = INDEX(CANA,'=')
         NANA(IANA,2) = ILEFT
         NANA(IANA,3) = IPRO
         NANA(IANA,4) = 1
         NANA(IANA,5) = IANA - 1
         NANA(IANA,6) = IANA + 1
         GOTO 100
      ENDIF
*
  20  I = I + 1
      IF(A(I:I).EQ.' ') GOTO 20
      IF(INDEX(LET,A(I:I)).EQ.0) THEN
         AER = 'ERROR, ILLEGAL VARIABLE NAME'
         GOTO 1000
      ENDIF
      CSEA = BLANK
      IS = 1
      CSEA(IS:IS) = A(I:I)
  21  I = I + 1
      IF(A(I:I).EQ.' ') GOTO 21
      IF(INDEX(LET//NUM//'_',A(I:I)).EQ.0.OR.IS.GT.LNAM) GOTO 22
      IS = IS + 1
      CSEA(IS:IS) = A(I:I)
      GOTO 21
  22  CONTINUE
      DO 23 ILEFT=ISYM,LCON+1,-1
      IF(CSEA.EQ.CSYM(ILEFT)) GOTO 24
  23  CONTINUE
      AER = 'ERROR, UNKNOWN VARIABLE'
      GOTO 1000
  24  CONTINUE
      IF(NSYM(ILEFT,1).NE.0) THEN
         AER = 'ERROR, EXPECTED VARIABLE '
         GOTO 1000
      ENDIF
      ILEFTO = NSYM(ILEFT,2)
      ILEFT  = NSYM(ILEFT,3)
  25  CONTINUE
      IF(A(I:I).EQ.' ') THEN
         I = I + 1
         GOTO 25
      ENDIF
      IF(A(I:I).EQ.':') THEN
         CONTINUE
      ELSEIF(A(I:I).EQ.'(') THEN
         IFUVA = 1
      ELSE
         AER = 'ERROR, EXPECTED OPENING PARENTHESIS OR ASSIGNMENT'
         GOTO 1000
      ENDIF
      GOTO 100
*
*     EQUAL SIGN
*     **********
*
  50  CONTINUE
      IF(A(I:I+1).NE.':=') THEN
         AER = 'ERROR, EXPECTED ASSIGNMENT 1'
         GOTO 1000
      ELSEIF(IEQU.NE.0) THEN
         AER = 'ERROR, ONLY ONE ASSIGNMENT'
         GOTO 1000
      ELSEIF(IPAR.NE.0) THEN
         AER = 'ERROR, UNBALANCED PARENTHESIS AT ASSIGNMENT'
         GOTO 1000
      ENDIF
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,'=')
      NANA(IANA,2) = ILEFT
      NANA(IANA,3) = ILEFTO
      NANA(IANA,4) = 1
      NANA(IANA,5) = IANA - 1
      NANA(IANA,6) = IANA + 1
      IEQU = 1
      I = I + 1
  52  I = I + 1
      IF(A(I:I).EQ.' ') GOTO 52
      GOTO 100
*
*     SWITCH TO EQUAL, OPEN PARENTHESIS, SIGN, CONSTANT, FUNCTION OR VARIABLE
*     ***********************************************************************
*
 100  CONTINUE
      IF(A(I:I).EQ.' ') GOTO 100
      IF(A(I:I).EQ.':') GOTO 50
      IF(A(I:I).EQ.'(') GOTO 110
      IF(A(I:I).EQ.'+'.OR.A(I:I).EQ.'-') GOTO 170
      IF(INDEX(NUM,A(I:I)).NE.0.OR.A(I:I).EQ.'.'.OR.A(I:I).EQ.'''')
     *                    GOTO 120
      IF(A(I:I).EQ.'[') GOTO 130
      GOTO 101
*
*     DECODING NAME OF FUNCTION OR VARIABLE
*     *************************************
*
 101  CONTINUE
      IF(INDEX(LET,A(I:I)).EQ.0) THEN
         AER = 'ERROR, ILLEGAL NAME OF VARIABLE OR FUNCTION'
         GOTO 1000
      ENDIF
      CSEA = BLANK
      CSEA(1:1) = A(I:I)
      IS = 1
 102  I = I + 1
      IF(A(I:I).EQ.' ') GOTO 102
      IS = IS + 1
      IF(INDEX(LET//NUM//'_',A(I:I)).EQ.0.OR.IS.GT.LNAM) GOTO 103
      CSEA(IS:IS) = A(I:I)
      GOTO 102
 103  CONTINUE
      DO 104 III=ISYM,LCON+1,-1
      IF(CSEA.EQ.CSYM(III)) THEN
         IF(NSYM(III,1).EQ.0) THEN
            II = NSYM(III,3)
            IP = NSYM(III,2)
            GOTO 160
         ELSEIF(NSYM(III,1).GT.0) THEN
            II = -NSYM(III,1)
            IP = NSYM(III,2)
            CREF(IREF+1:IREF+ILAST(CSYM(III),1,32)+3) = ' %'//
     *         CSYM(III)(1:ILAST(CSYM(III),1,32))//':'
               IREF = IREF + ILAST(CSYM(III),1,32) + 3
            IF(NSYM(III,4).LT.100000) THEN
               WRITE(CREF(IREF+1:IREF+5),'(I5)') NSYM(III,4)
               IREF = IREF + 5
            ELSE
               WRITE(CREF(IREF+1:IREF+8),'(I8)') NSYM(III,4)
               IREF = IREF + 8
            ENDIF
            GOTO 150
         ELSE
            AER = 'ERROR, EXPECTED VARIABLE OR FUNCTION'
            GOTO 1000
         ENDIF
      ENDIF
 104  CONTINUE
      II = (INDEX(FUN,CSEA(1:6))+5)/6
      IP = 0
      IF(II.NE.0) GOTO 150
      AER = 'ERROR, UNKNOWN VARIABLE OR FUNCTION'
      GOTO 1000
*
*     OPEN PARENTHESIS
*     ****************
*
 110  CONTINUE
      IF(A(I:I).NE.'(') THEN
         AER = 'ERROR, EXPECTED OPENING PARENTHESIS'
         GOTO 1000
      ENDIF
 112  I = I + 1
      IF(A(I:I).EQ.' ') GOTO 112
      IPAR = IPAR + 1
      IF(IPAR.GT.20) THEN
         AER = 'ERROR, PARENTHESES NESTED TOO DEEPLY'
         GOTO 1000
      ENDIF
      GOTO 100
*
*     CONSTANT
*     ********
*
 120  CONTINUE
      ICON = ICON + 1
      IF(ICON.GT.LCON) THEN
         CALL MSG( '!!! MEMORY EXHAUSTION, INCREASE LCON')
         CALL FOXSTP(1)
      ENDIF
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,'C')
      NANA(IANA,2) = ICON
      NANA(IANA,3) = 0
      NANA(IANA,4) = 0
      NANA(IANA,5) = IANA - 1
      NANA(IANA,6) = IANA + 1
      IF(A(I:I).EQ.'''') THEN
         ICHA = 1
         CCON(ICON)(1:1) = '$'
         JCON = 1
 125     I = I + 1
         IF(JCON.GE.512) THEN
            DCON(ICON) = JCON
            I = I + 1
            GOTO 200
         ELSEIF(A(I:I).EQ.'''') THEN
            IF(ICHA.EQ.1) THEN
               IF(A(I+1:I+1).EQ.'''') THEN
                  ICHA = 2
               ELSE
                  DCON(ICON) = JCON
                  I = I + 1
                  GOTO 200
               ENDIF
            ELSEIF(ICHA.EQ.2) THEN
               ICHA = 1
            ENDIF
         ENDIF
         IF(ICHA.EQ.1) THEN
            JCON = JCON + 1
            CCON(ICON)(JCON:JCON) = A(I:I)
         ENDIF
         GOTO 125
      ELSE
         CCON(ICON)(1:1) = ' '
      ENDIF
      DO 121 IJ=1,LDEC
      CVAL(IJ:IJ) = ' '
 121  CONTINUE
      IJ = 0
      I = I - 1
 122  I = I + 1
      IF(A(I:I).EQ.' ') GOTO 122
      IF(INDEX(NUM//'.DE',A(I:I)).NE.0) THEN
         IJ = IJ+1
         IF(IJ.GT.LDEC) GOTO 124
         CVAL(IJ:IJ) = A(I:I)
         GOTO 122
      ELSEIF(INDEX('+-',A(I:I)).NE.0) THEN
         IF(INDEX('DE',A(I-1:I-1)).EQ.0) GOTO 123
         IJ = IJ+1
         IF(IJ.GT.LDEC) GOTO 124
         CVAL(IJ:IJ) = A(I:I)
         GOTO 122
      ELSE
         GOTO 123
      ENDIF
      GOTO 122
 123  CONTINUE
      CALL VALCH(CVAL,1,IJ,VAL,IER)
      IF(IER.EQ.1) THEN
         I = I-1
         AER = 'ERROR, CONSTANT EXPRESSION IS WRONG'
         GOTO 1000
      ENDIF
      DCON(ICON) = VAL
      GOTO 200
 124  CONTINUE
      AER = 'ERROR, TOO MANY DIGITS IN CONSTANT'
      GOTO 1000
*
*     ADDITIONAL OPTION
*     *****************
*
 130  CONTINUE
      AER = 'ERROR, OPERATION NOT SUPPORTED YET'
      GOTO 1000
*
*     CLOSING PARENTHESIS
*     *******************
*
 140  CONTINUE
      IF(A(I:I).NE.')') THEN
         AER = 'ERROR, EXPECTED CLOSING PARENTHESIS'
         GOTO 1000
      ENDIF
      IF(IPAR.EQ.IFUVA) IFUVA = 0
      IPAR = IPAR - 1
      IF(IPAR.LT.0) THEN
         AER = 'ERROR, PARENTHESIS CLOSED TOO EARLY'
         GOTO 1000
      ENDIF
 141  I = I + 1
      IF(A(I:I).EQ.' ') GOTO 141
      IF(IPAR.EQ.0.AND.IEQU.EQ.0.AND.A(I:I).NE.':') THEN
         AER = 'ERROR, EXPECTED ASSIGNMENT 2'
         GOTO 1000
      ENDIF
      GOTO 200
*
*     FUNCTION
*     ********
*
 150  CONTINUE
      IF(A(I:I).NE.'(') THEN
         AER = 'ERROR, EXPECTED OPENING PARENTHESIS'
         GOTO 1000
      ENDIF
      IF(IFUVA.EQ.0) IFUVA = IPAR + 1
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,'F')
      NANA(IANA,2) = II
      NANA(IANA,3) = IP
      NANA(IANA,4) = MPRI*(IPAR+1)+MPRI-1
      NANA(IANA,5) = IANA - 1
      NANA(IANA,6) = IANA + 1
      IOP = MAX(IOP,MPRI*(IPAR+1)+MPRI-1)
      GOTO 110
*
*     VARIABLE
*     ********
*
 160  CONTINUE
      IF(A(I:I).NE.'(') THEN
         IANA = IANA + 1
         NANA(IANA,1) = INDEX(CANA,'V')
         NANA(IANA,2) = II
         NANA(IANA,3) = IP
         NANA(IANA,4) = 0
         NANA(IANA,5) = IANA - 1
         NANA(IANA,6) = IANA + 1
         GOTO 200
      ELSE
         IF(IFUVA.EQ.0) IFUVA = IPAR + 1
         IANA = IANA + 1
         NANA(IANA,1) = INDEX(CANA,'A')
         NANA(IANA,2) = II
         NANA(IANA,3) = IP
         NANA(IANA,4) = MPRI*(IPAR+1)+MPRI
         NANA(IANA,5) = IANA - 1
         NANA(IANA,6) = IANA + 1
         IOP = MAX(IOP,MPRI*(IPAR+1)+MPRI)
         GOTO 110
      ENDIF
*
*     PLUS OR MINUS SIGN IN FRONT OF VARIABLE,FUNCTION OR CONSTANT
*     ************************************************************
*
 170  CONTINUE
      IF(A(I:I).EQ.'+') THEN
 172     I = I + 1
         IF(A(I:I).EQ.' ') GOTO 172
         GOTO 100
      ELSEIF(A(I:I).EQ.'-') THEN
 174     I = I + 1
         IF(A(I:I).EQ.' ') GOTO 174
         ICON = ICON + 1
         DCON(ICON)      = -1.D0
         CCON(ICON)(1:1) = ' '
         IANA = IANA + 1
         NANA(IANA,1) = INDEX(CANA,'C')
         NANA(IANA,2) = ICON
         NANA(IANA,3) = 0
         NANA(IANA,4) = 0
         NANA(IANA,5) = IANA - 1
         NANA(IANA,6) = IANA + 1
         IANA = IANA + 1
         NANA(IANA,1) = INDEX(CANA,'O')
         NANA(IANA,2) = 3
         NANA(IANA,3) = MPRI*(IPAR+1)+NPRI(INDEX(OPS,'-'))
         NANA(IANA,4) = MPRI*(IPAR+1)+NPRI(INDEX(OPS,'-'))
         NANA(IANA,5) = IANA - 1
         NANA(IANA,6) = IANA + 1
         IB = IB + 10
         IOP = MAX(IOP,MPRI*(IPAR+1)+NPRI(INDEX(OPS,'-')))
         GOTO 100
      ELSE
         AER = 'ERROR WITH PLUS OR MINUS SIGN'
         GOTO 1000
      ENDIF
*
*     COMMA (IN FUNCTION OR ARRAY)
*     ****************************
*
 180  CONTINUE
      IF(A(I:I).NE.',') THEN
         AER = 'ERROR, COMMA EXPECTED'
         GOTO 1000
      ENDIF
      IF(IFUVA.EQ.0) THEN
         AER = 'ERROR, COMMA MISPLACED'
         GOTO 1000
      ENDIF
 182  I = I + 1
      IF(A(I:I).EQ.' ') GOTO 182
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,',')
      NANA(IANA,2) = 0
      NANA(IANA,3) = 0
      NANA(IANA,4) = MPRI*(IPAR+1)
      NANA(IANA,5) = IANA - 1
      NANA(IANA,6) = IANA + 1
      GOTO 100
*
*     SWITCH TO OPERATOR, CLOSING PARENTHESIS, COMMA OR END
*     *****************************************************
*
 200  CONTINUE
 201  CONTINUE
      IF(A(I:I).EQ.' ') THEN
         I = I + 1
         GOTO 201
      ENDIF
      CSEA = BLANK
      IF(INDEX(OPS,A(I:I)).NE.0) GOTO 210
      IF(A(I:I).EQ.':') GOTO 50
      IF(A(I:I).EQ.')') GOTO 140
      IF(A(I:I).EQ.',') GOTO 180
      IF(A(I:I).EQ.';') GOTO 220
      AER = 'ERROR, OPERATOR, COMMA OR CLOSING PARENS EXPECTED'
      GOTO 1000
*
*     OPERATORS
*     *********
*
 210  CONTINUE
      II = INDEX(OPS,A(I:I))
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,'O')
      NANA(IANA,2) = II
      NANA(IANA,3) = MPRI*(IPAR+1)+NPRI(II)
      NANA(IANA,4) = MPRI*(IPAR+1)+NPRI(II)
      NANA(IANA,5) = IANA - 1
      NANA(IANA,6) = IANA + 1
      IOP = MAX(IOP,MPRI*(IPAR+1)+NPRI(II))
 212  I = I + 1
      IF(A(I:I).EQ.' ') GOTO 212
      GOTO 100
*
*     END OF SOURCE CODE
*     ******************
*
 220  CONTINUE
      IF(IPAR.NE.0) THEN
         AER = 'ERROR, UNBALANCED PARENTHESES'
         GOTO 1000
      ENDIF
      IF(IEQU.EQ.0) THEN
         AER = 'ERROR, ASSIGNMENT INCOMPLETE'
         GOTO 1000
      ENDIF
      IANA = IANA + 1
      NANA(IANA,1) = INDEX(CANA,';')
      NANA(IANA,2) = 0
      NANA(IANA,3) = 0
      NANA(IANA,4) = 0
      NANA(IANA,5) = IANA - 1
      NANA(IANA,6) = 0
*
      GOTO 300
*
*     EXTRACTION OF OPERATIONS FROM INTERIM CODE
*     ******************************************
*
 300  CONTINUE
*
      IF(LCHECK.EQ.1) WRITE(6,'(1X,A)') 'BEGINNING EXTRACTION'
*
      DO 350 IO = IOP,1,-1
*     ********************
*
*     EXTRACTING NEXT ELEMENTARY OPERATION
*     ************************************
*
      IA = 1
 310  CONTINUE
      IF(NANA(IA,4).NE.IO) THEN
         IA = NANA(IA,6)
         IF(IA.EQ.IANA) GOTO 350
         GOTO 310
      ELSEIF(CANA(NANA(IA,1):NANA(IA,1)).EQ.',') THEN
         IA = NANA(IA,6)
         IF(IA.EQ.IANA) GOTO 350
         GOTO 310
      ENDIF
*
      IF(LCHECK.EQ.1) THEN
         WRITE(6,'(A,I3,71A)') ' IO = ',IO , ' ',('-',J=1,70)
         WRITE(6,'(100(6X,6(A1,I5,I5,A1)/))')
     *   (CANA(NANA(J,1):NANA(J,1)),NANA(J,2),NANA(J,4),'|',J=1,IANA)
      ENDIF
*
*     FILLING IN FIRST ENTRIES OF NCOD
*     ********************************
*
      ITYP = NANA(IA,1)
*
      NCOD(ICOD+2) = 0
      NCOD(ICOD+3) = 0
      NCOD(ICOD+4) = 0
      NCOD(ICOD+5) = 0
      NCOD(ICOD+6) = 0
      NCOD(ICOD+7) = 0
*
      IF(CANA(ITYP:ITYP).EQ.'=') THEN
         NCOD(ICOD+1) = LSYNR+LOPS+1
         NCOD(ICOD+6) = NANA(IA,2)
      ELSEIF(CANA(ITYP:ITYP).EQ.'O') THEN
         NCOD(ICOD+1) = NANA(IA,2) + LSYNR
      ELSEIF(CANA(ITYP:ITYP).EQ.'A') THEN
         NCOD(ICOD+1) = 351
         NCOD(ICOD+7) = NANA(IA,2)
      ELSEIF(CANA(ITYP:ITYP).EQ.'F') THEN
         IF(NANA(IA,2).GT.0) THEN
            NCOD(ICOD+1) = NANA(IA,2) +LSYNR+LOPS
         ELSE
            NCOD(ICOD+1) = 353
            NCOD(ICOD+3) = -NANA(IA,2)
         ENDIF
      ELSE
         AER = 'ERROR, ILLEGAL OPERAND'
         GOTO 2000
      ENDIF
*
      IF(CANA(ITYP:ITYP).NE.'=') THEN
         ISCR = ISCR + 1
         IF(ISCR.GT.LSCR) THEN
            CALL MSG( '### EXPRESSION TOO LONG ')
            LERR = 1
            RETURN
         ENDIF
         MSCR = MAX(MSCR,ISCR)
         NCOD(ICOD+6) = NPRO(IPRO,1) + ISCR
      ENDIF
*
*     FINDING LEFT OPERANDS
*     *********************
*
      JL = NANA(IA,5)
      IF(CANA(NANA(IA,1):NANA(IA,1)).EQ.'='.AND.JL.EQ.1) THEN
         IL = 0
         GOTO 325
      ELSEIF(CANA(NANA(IA,1):NANA(IA,1)).EQ.'A') THEN
         IL = 0
         GOTO 325
      ELSEIF(CANA(NANA(IA,1):NANA(IA,1)).EQ.'F') THEN
         IL = 0
         GOTO 325
      ENDIF
      IL = 1
*
 320  CONTINUE
      IF(INDEX('CVS',CANA(NANA(JL,1):NANA(JL,1))).EQ.0) THEN
         AER = 'ERROR, NO LEFT OPERAND FOUND'
         GOTO 2000
      ENDIF
      IF(CANA(ITYP:ITYP).EQ.'=') THEN
         NCOD(ICOD+7+IL) = NANA(JL,2)
      ELSE
         NCOD(ICOD+7) = NANA(JL,2)
      ENDIF
*
      NANA(IA,5) = NANA(JL,5)
      NANA(NANA(JL,5),6) = IA
*
      NANA(JL,1) = INDEX(CANA,'#')
      NANA(JL,2) = 0
      NANA(JL,3) = 0
      NANA(JL,4) = 0
      NANA(JL,5) = 0
      NANA(JL,6) = 0
*
      JL = NANA(IA,5)
*
      IF(JL.EQ.1) THEN
         CONTINUE
      ELSEIF(CANA(NANA(JL,1):NANA(JL,1)).EQ.','.AND.
     *   NANA(JL,4).GE.NANA(IA,4))  THEN
         NANA(JL,1) = INDEX(CANA,'#')
         IL = IL + 1
         JL = NANA(JL,5)
         GOTO 320
      ENDIF
 325  CONTINUE
*
      IF(IL.GT.1.AND.CANA(ITYP:ITYP).NE.'=') THEN
         AER = 'ERROR, TOO MANY LEFT OPERANDS'
         GOTO 2000
      ENDIF
*
*     FINDING RIGHT OPERANDS
*     **********************
*
      JR = NANA(IA,6)
      IR = 1
 340  CONTINUE
      IF(INDEX('CVS',CANA(NANA(JR,1):NANA(JR,1))).EQ.0) THEN
         AER = 'ERROR, NO RIGHT OPERAND FOUND'
         GOTO 2000
      ENDIF
      IF(CANA(ITYP:ITYP).EQ.'=') THEN
         NCOD(ICOD+7) = NANA(JR,2)
      ELSEIF(CANA(ITYP:ITYP).EQ.'F') THEN
         NCOD(ICOD+6+IR) = NANA(JR,2)
      ELSE
         NCOD(ICOD+7+IR) = NANA(JR,2)
      ENDIF
*
      NANA(IA,6) = NANA(JR,6)
      NANA(NANA(JR,6),5) = IA
*
      NANA(JR,1) = INDEX(CANA,'#')
      NANA(JR,2) = 0
      NANA(JR,3) = 0
      NANA(JR,4) = 0
      NANA(JR,5) = 0
      NANA(JR,6) = 0
*
      JR = NANA(IA,6)
*
      IF(JR.EQ.IANA) THEN
         CONTINUE
      ELSEIF(CANA(NANA(JR,1):NANA(JR,1)).EQ.','.AND.
     *   NANA(JR,4).GE.NANA(IA,4)) THEN
         NANA(JR,1) = INDEX(CANA,'#')
         JR = NANA(JR,6)
         IR = IR + 1
         GOTO 340
      ENDIF
*
*     CHECKS AND UPDATING OF ICOD
*     ***************************
*
      IF(IR.NE.1.AND.CANA(ITYP:ITYP).NE.'A'
     *          .AND.CANA(ITYP:ITYP).NE.'F') THEN
         AER = 'ERROR, TOO MANY RIGHT OPERANDS'
         GOTO 2000
      ENDIF
*
      IF(CANA(NANA(IA,1):NANA(IA,1)).EQ.'=') THEN
         IF(IL.NE.0) NCOD(ICOD+1) = 352
         INCRO = INCR
         INCR  = 7 + IL
      ELSEIF(CANA(NANA(IA,1):NANA(IA,1)).EQ.'F') THEN
         INCRO = INCR
         INCR  = 6 + IR
      ELSE
         INCRO = INCR
         INCR  = 7 + IR
      ENDIF
*
      NCOD(ICOD+2) = ICOD + INCR
      NCOD(ICOD+4) = ICOD + INCR
      NCOD(ICOD+5) = -ILIS
      ICOO = ICOD
      ICOD = NCOD(ICOD+2)
*
*     REPLACING OPERATION JUST PROCESSED BY SCRATCH ADDRESS
*     *****************************************************
*
      NANA(IA,1) = INDEX(CANA,'S')
      NANA(IA,2) = NPRO(IPRO,1) + ISCR
      NANA(IA,3) = IPRO
      NANA(IA,4) = 0
*
      GOTO 310
*
 350  CONTINUE
*
      IF(LCHECK.EQ.1) THEN
         WRITE(6,'(1X,79A)') ('-',J=1,79)
         WRITE(6,'(100(6X,6(A1,I5,I5,A1)/))')
     *   (CANA(NANA(J,1):NANA(J,1)),NANA(J,2),NANA(J,4),'|',J=1,IANA)
      ENDIF
*
*     ALL OPERATIONS EXTRACTED
*     ************************
*
 400  CONTINUE
*
      IF(NANA(1,5).NE.NANA(IANA,4)) THEN
         WRITE(2,'(1X,2A)') '@@@ ERROR, CODE NOT FULLY PROCESSED '
         WRITE(6,'(1X,2A)') '@@@ ERROR, CODE NOT FULLY PROCESSED '
         CALL FOXSTP(1)
      ENDIF
*
      RETURN
*
*     SYNTAX ERROR EXIT
*     *****************
*
 1000 CONTINUE
      WRITE(2,'(1X,2A)') '### ', AER
      WRITE(2,'(1X,10(79A/))') ' ', (A(K:K),K=1,IA2)
      WRITE(2,'(1X,10(79A/))') ' ', (' ',K=1,I-1),'*'
      WRITE(6,'(1X,2A)') '### ', AER
      WRITE(6,'(1X,10(79A/))') ' ', (A(K:K),K=1,IA2)
      WRITE(6,'(1X,10(79A/))') ' ', (' ',K=1,I-1),'*'
      LERR = 1
      RETURN
*
 2000 CONTINUE
      CALL MSG( '@@@ INTERNAL ERROR WHILE EXTRACTING CODE:')
      CALL MSG( AER )
      WRITE(2,'(100(6X,6(A1,I5,I5,A1)/))')
     *     (CANA(NANA(J,1):NANA(J,1)),NANA(J,2),NANA(J,4),'|',J=1,IANA)
      WRITE(2,'(100(6X,6A12/))') ('        ',K=1,IA-1),'*       '
      WRITE(6,'(100(6X,6(A1,I5,I5,A1)/))')
     *     (CANA(NANA(J,1):NANA(J,1)),NANA(J,2),NANA(J,4),'|',J=1,IANA)
      WRITE(6,'(100(6X,6A12/))') ('        ',K=1,IA-1),'*       '
*
      CALL FOXSTP(1)
      END
*
*FACE SUBROUTINE CODEXE(NOP,INC,INA,INB)                                 *FACE
      SUBROUTINE CODEXE(NCOD, DCON,CCON, LERR)                           *NORM
*MPI  SUBROUTINE CODEXE(NCOD, DCON,CCON, LERR)                           *MPI
*     ****************************************
*
*     THIS SUBROUTINE INTERPRETS AND EXECUTES THE INTERIM CODE IN THE
*     ARRAY NCOD
*
      PARAMETER(LADD=100000,LNPRO=30,LMPRO=2000)
      PARAMETER(LOPT=100,LOPV=20,LLOO=500)
      PARAMETER(LSAL=100)
      PARAMETER(LCOMM=50,LPLO=50,LPLV=10,LPEN=1024)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*MPI  INCLUDE 'mpif.h'                                                   *MPI
*
      DIMENSION NCOD(*), DCON(*)                                         *NORM
      CHARACTER CCON(*)*512                                              *NORM
*MPI  DIMENSION NCOD(*), DCON(*)                                         *MPI
*MPI  CHARACTER CCON(*)*512                                              *MPI
*FACE DIMENSION NCOD(10), DCON(10)                                       *FACE
*FACE CHARACTER CCON(10)*512                                             *FACE
*
      DIMENSION NADD(LADD), NPRO(LNPRO,6), MPRO(LMPRO)
*
      DIMENSION NOPT(LOPT,LOPV),XOPT(LOPV,LOPT),XXOPT(LOPV)
      DIMENSION XOBJ(LOPV),XLOO(LLOO,4)
      DIMENSION NNCOM(8,LCOMM),NPLO(4+LPLV,0:LPLO),NNPL(2,0:LPLO)
      DIMENSION MPLO(LNPRO),NMPL(4,LPLO),NRECV(LPEN),NDISP(LPEN)
*
*     NADD: DYNAMICAL ABSOLUTE ADDRESS FOR EACH VARIABLE IN CODE
*
*     NPRO: IPRO DYNAMIC PROCEDURE INFORMATION
*           1: ADDRESS OF FUNCTION BEGINNING
*           2: ADDRESS FROM WHERE CALLED
*           3: IPRO OF LAST PROCEDURE WITH SAME ID (IN CASE OF RECURSION)
*           4: IVAR AT BEGINNING OF ROUTINE
*           5: IMEM AT BEGINNING OF ROUTINE
*           6: IDIM AT BEGINNING OF ROUTINE
*
*     MPRO: JPRO
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
*----------------------------------------------------------------------------
*
      COMMON /CODEX/ ILIS,NSCR
      COMMON /MSTAT/ MMEM,MVAR
      COMMON /PLCOM/ MTASK,ISROOT
*
*     NTYP: INTEGER DESCRIBING THE TYPE
*     NBEG: BEGINNING ADDRESS IN CC, NC
*     NEND: ENDING ADDRESS IN CC, NC (MOMENTARY)
*     NMAX: MAXIMUM ENDING ADDRESS RESERVED
*
*     CC:   CONTAINS IMEM DOUBLE PRECISION NUMBERS
*     NC:   CONTAINS IMEM POSITION INTEGERS
*     IVAR: MOMENTARY NUMBER OF ALLOCATED VARIABLES
*
      PARAMETER(LTYP=11,LSYNR=30,LOPS=20,LFUN=100,LSUB=200)
      CHARACTER CTYID(LTYP)*2
      COMMON /TYCID/ CTYID
*
*     COMPUTED GOTO FOR THE INTERFACES - BYPASSES THE REGULAR GOTO
*     ************************************************************
*
*FACE IC      = 1                                                        *FACE
*FACE NCOD(3) = 1                                                        *FACE
*
*FACE GOTO(1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040,   *FACE
*FACE*     1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050,   *FACE
*FACE*     1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060,   *FACE
*FACE*     1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070,   *FACE
*FACE*     1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080,   *FACE
*FACE*     1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090,   *FACE
*FACE*     1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100,   *FACE
*FACE*     1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110,   *FACE
*FACE*     1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120,   *FACE
*FACE*     1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130,   *FACE
*FACE*     1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140,   *FACE
*FACE*     1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150),  *FACE
*FACE*     NOP                                                           *FACE
*FACE GOTO 2000                                                          *FACE
*
      IF(LERR.NE.0) THEN
         WRITE(6,'(1X,A)') '### ERRORS IN CODE, NO EXECUTION'
         RETURN
      ENDIF
*
*      WRITE(6,'(1X,A)') '--- BEGINNING EXECUTION'                        *NORM
*
      CALL FOXKEY('FOX V10.0')                                             $$$
      CALL FOXPAR(LMEM,LVAR,LDIM)
      CALL GRAPAR(LMEM,LVAR,LDIM)
*
      IOPT = 0
      ILOO = 0
      IPLO = 0
      JPLO = 0
      JCOMM = 0
*
      IC   = 0
*
 1000 CONTINUE
**************
*FACE RETURN                                                             *FACE
*
      IGOTO = NCOD(IC+1)
      ILISN  = NCOD(IC+5)
      IF(ILISN.LT.0) ILIS = ILISN
*
      GOTO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,2
     *4,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,4
     *6,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,6
     *8,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,9
     *0,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,1
     *09,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125
     *,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,1
     *42,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158
     *,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,1
     *75,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191
     *,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,2
     *08,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224
     *,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,2
     *41,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257
     *,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,2
     *74,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290
     *,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,3
     *07,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323
     *,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,3
     *40,341,342,343,344,345,346,347,348,349,350,351,352,353),IGOTO
*
      GOTO 2000
*
***** ASS **
    1 CONTINUE
      IC = NCOD(IC+2)
      GOTO 1000
*
***** B E **
    2 CONTINUE
*
      IVAR = 0
      IMEM = 0
      MMEM = 0
      MVAR = 0
*
      IPRO = 1
      NPRO(IPRO,1) = 0
      NPRO(IPRO,2) = 0
      NPRO(IPRO,3) = 0
      NPRO(IPRO,4) = 0
      NPRO(IPRO,5) = 0
      NPRO(IPRO,6) = 0
      MPRO(1) = 1
*
      MSCR = NCOD(IC+6)
      LSCR = NCOD(IC+7)
      NSCR = NCOD(IC+8)
      MCON = NCOD(IC+9)
      LCON = NCOD(IC+10)
      JPRO = NCOD(IC+11)
*
      IF(JPRO.GT.LMPRO) THEN
         WRITE(6,'(1X,A,I8)')
     *      '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LMPRO BEYOND ',JPRO
         GOTO 3000
      ENDIF
*
      DO 1002 J=1,JPRO
 1002 MPRO(J) = 0
*
      MPICW = 1
      MPICS = 1
      MPICN = -1
      MPIUN = 0
      MPERR = 0
      MGRO = 1
      MTASK = 1
      MRANK = 0
      NGRP = 0
      MEXE = 1
*
*MPI  MPICW = MPI_COMM_WORLD                                             *MPI
*MPI  MPICS = MPI_COMM_SELF                                              *MPI
*MPI  MPICN = MPI_COMM_NULL                                              *MPI
*MPI  MPIUN = MPI_UNDEFINED                                              *MPI
*MPI  CALL MPI_COMM_GROUP(MPI_COMM_WORLD,MGRO,MPERR)                     *MPI
*MPI  CALL MPI_GROUP_SIZE(MGRO,MTASK,MPERR)                              *MPI
*MPI  CALL MPI_GROUP_RANK(MGRO,MRANK,MPERR)                              *MPI
*MPI  NGRP = 1                                                           *MPI
*MPI  MEXE = 2                                                           *MPI
*
      MROOT = 0
      IROOT = MROOT + 1
*
      JCOMM = JCOMM + 1
      NNCOM(1,JCOMM) = MTASK
      NNCOM(2,JCOMM) = MTASK
      NNCOM(3,JCOMM) = MRANK + 1
      NNCOM(4,JCOMM) = 1
      NNCOM(5,JCOMM) = 1
      NNCOM(6,JCOMM) = MPICW
      NNCOM(7,JCOMM) = MPICS
      NNCOM(8,JCOMM) = 1
*
      JCOMM = JCOMM + 1
      NNCOM(1,JCOMM) = 1
      NNCOM(2,JCOMM) = 1
      NNCOM(3,JCOMM) = 1
      NNCOM(4,JCOMM) = MTASK
      NNCOM(5,JCOMM) = MRANK + 1
      NNCOM(6,JCOMM) = MPICS
      IF(NNCOM(5,JCOMM).NE.IROOT) NNCOM(6,JCOMM) = MPICN
      NNCOM(7,JCOMM) = MPICW
      NNCOM(8,JCOMM) = 1
      IF(NNCOM(6,JCOMM).EQ.MPICN) NNCOM(8,JCOMM) = 0
      ISROOT = NNCOM(8,JCOMM)
*
      NPLO(1,JPLO) = IC
      NPLO(2,JPLO) = 0
      NPLO(3,JPLO) = 0
      NPLO(4,JPLO) = 0
      NNPL(1,IPLO) = JPLO
      NNPL(2,IPLO) = JCOMM
*
      MPLO(IPRO) = JPLO
*
*     ALLOCATING CONSTANTS
*     ********************
*
      IMEM = 0
      DO 2002 J=1,MCON
      IF(CCON(J)(1:1).EQ.' ') THEN
         IMEM = IMEM + 1
         NTYP(J) = 1
         NBEG(J) = IMEM
         NEND(J) = IMEM
         NMAX(J) = IMEM
         CC(IMEM) = DCON(J)
      ELSE
         NBEG(J) = IMEM + 1
         NTYP(J) = 2
         DO 3002 K=2,NINT(DCON(J))
         IMEM = IMEM + 1
 3002    NC(IMEM) = ICHAR((CCON(J)(K:K)))
         NEND(J) = IMEM
         NMAX(J) = IMEM
      ENDIF
      NADD(J) = J
 2002 CONTINUE
*
      IVAR = MCON
      MMEM = MAX(MMEM,IMEM)
      CALL VARCHK(IVAR)
      CALL MEMCHK(IMEM)
*
      JADD = LCON + 2
      IC = NCOD(IC+2)
      GOTO 500
*
***** E E **
    3 CONTINUE
      JPLO = MPLO(IPRO)
      IPRO = IPRO - 1
      IF(IPRO.NE.0) THEN
         WRITE(6,'(1X,A,I3)') '@@@ RUNTIME ERROR, IPRO = ',IPRO
      ENDIF
      RETURN
*
***** D V **
    4 CONTINUE
*
      JVAR = NCOD(IC+6)
      IF(JVAR.GE.LADD.OR.JVAR.LE.0) THEN
         WRITE(6,'(1X,A)')
     *   '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LADD'
         GOTO 3000
      ENDIF
      JLEN = NADD(NCOD(IC+7))
      IF(NTYP(JLEN).NE.1) THEN
         WRITE(6,'(1X,A,I6)')
     *      '$$$ VARIABLE LENGTH WRONG TYPE ',NTYP(JLEN)
         GOTO 3000
      ENDIF
*
      JLEN = NINT(CC(NBEG(JLEN)))
      IF(JLEN.LE.0) THEN
         WRITE(6,'(1X,A,I6)')
     *      '$$$ VARIABLE LENGTH MUST BE POSITIVE INTEGER'
         GOTO 3000
      ENDIF
*
      JNUM = 1
      JDIM = NCOD(IC+4) - IC - 7
      IF(IDIM+JDIM.GE.LDIM) THEN
         WRITE(6,'(1X,A)')
     *       '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LDIM'
         GOTO 3000
      ENDIF
      IDIMO = IDIM + 1
      DO 1004 J=1,JDIM
      KDIM = NADD(NCOD(IC+7+J))
      IF(NTYP(KDIM).NE.1) THEN
         WRITE(6,'(1X,A,I6)')
     *      '$$$ ARRAY INDEX HAS WRONG TYPE ',NTYP(KDIM)
         GOTO 3000
      ENDIF
      IDIM = IDIMO + J
      NDIM(IDIM) = NINT(CC(NBEG(KDIM)))
      JNUM = JNUM * NDIM(IDIM)
 1004 CONTINUE
*
      IF(JDIM.EQ.0) THEN
         IVAR = IVAR + 1
         CALL VARCHK(IVAR)
         NTYP(IVAR) = 1
         NBEG(IVAR) = IMEM + 1
         NMAX(IVAR) = IMEM + JLEN
         NEND(IVAR) = IMEM + 1
         NADD(JVAR) = IVAR
         CC(IMEM+1) = 0.D0
         IMEM = IMEM + JLEN
         CALL MEMCHK(IMEM)
      ELSE
         NDIM(IDIMO) = JNUM
         IVAR = IVAR + 1
         NTYP(IVAR) = -JDIM
         NBEG(IVAR) = 0
         NEND(IVAR) = JLEN
         NMAX(IVAR) = IDIMO
         NADD(JVAR) = IVAR
         CALL VARCHK(IVAR+JNUM)
         CALL MEMCHK(IMEM+JNUM*JLEN)
         DO 2004 J=1,JNUM
         IVAR = IVAR + 1
         NTYP(IVAR) = 1
         NBEG(IVAR) = IMEM + 1
         NMAX(IVAR) = IMEM + JLEN
         NEND(IVAR) = IMEM + 1
         CC(IMEM+1) = 0.D0
         IMEM = IMEM + JLEN
 2004    CONTINUE
      ENDIF
*
      IC = NCOD(IC+2)
      GOTO 1000
*
***** B F **
   5  CONTINUE
      IARG = NCOD(ICFU+4) - ICFU - 6
      IV = NCOD(IC+5)
      IF(IARG.NE.IV-IC-6) THEN
         WRITE(6,'(1X,A)') '$$$ WRONG NUMBER OF FUNCTION ARGUMENTS'
         WRITE(6,'(5X,A,I3,A,I3)') 'ARGUMENTS EXPECTED:',IV-IC-6,
     *                           ', SUPPLIED:',IARG
         GOTO 3000
      ENDIF
*
      JADD = NCOD(IC+6)
      DO 1005 J=1,IARG
 1005 NADD(JADD+J) = NADD(NCOD(ICFU+6+J))
*
      IV   = NCOD(IC+5)
      MSCR = NCOD(IV+1)
      JPRO = NCOD(IV+2)
      JADD = NCOD(IV+3)
      NADD(JADD) = NADD(NCOD(ICFU+6))
      NTYP(NADD(JADD)) = 1
      CC(NBEG(NADD(JADD))) = 0.D0
*
      IPRO = IPRO + 1
      IF(IPRO.GT.LNPRO) THEN
         WRITE(6,'(1X,A)')
     *      '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LNPRO'
         GOTO 3000
      ENDIF
*
      NPRO(IPRO,1) = IC
      NPRO(IPRO,2) = ICFU
      NPRO(IPRO,3) = MPRO(JPRO)
      NPRO(IPRO,4) = IVAR
      NPRO(IPRO,5) = IMEM
      NPRO(IPRO,6) = IDIM
      MPRO(JPRO) = IPRO
      MPLO(IPRO) = JPLO
*
      IC = NCOD(IC+2)
      GOTO 500
*
***** E F **
    6 CONTINUE
      JPRO = NPRO(IPRO,3)
      IVAR = NPRO(IPRO,4)
      IMEM = NPRO(IPRO,5)
      IDIM = NPRO(IPRO,6)
      JJPRO = NCOD(NCOD(NPRO(IPRO,1)+5)+2)
      MPRO(JJPRO) = JPRO
      IF(JPRO.NE.0) THEN
         ICO = NPRO(JPRO,1)
         ICFUO = NPRO(JPRO,2)
         JADD = NCOD(ICO+6)
         NADD(JADD) = NCOD(ICFUO+6)
         DO 1006 J=1,IARG
 1006    NADD(JADD+J) = NCOD(ICFUO+7+J)
         IVO  = NCOD(ICO+5)
         MSCR = NCOD(IVO+1)
         JADD = NCOD(IVO+3)
         JALL = NPRO(JPRO,4)
         DO 2006 J=1,MSCR
 2006    NADD(JADD+J) = JALL + J
      ENDIF
      JPLO = MPLO(IPRO)
      IC = NCOD(NPRO(IPRO,2)+2)
      IPRO = IPRO - 1
      GOTO 1000
*
***** B P **
    7 CONTINUE
      IARG = NCOD(ICPR+4) - ICPR - 5
      IV = NCOD(IC+5)
      IF(IARG.NE.IV-IC-6) THEN
         WRITE(6,'(1X,A)') '$$$ WRONG NUMBER OF PROCEDURE ARGUMENTS'
         WRITE(6,'(5X,A,I3,A,I3)') 'ARGUMENTS EXPECTED:',IV-IC-6,
     *                           ', SUPPLIED:',IARG
         GOTO 3000
      ENDIF
*
      JADD = NCOD(IC+6)
      DO 1007 J=1,IARG
 1007 NADD(JADD+J) = NADD(NCOD(ICPR+5+J))
*
      IV   = NCOD(IC+5)
      MSCR = NCOD(IV+1)
      JPRO = NCOD(IV+2)
      JADD = NCOD(IV+3)
      NADD(JADD) = 0
*
      IPRO = IPRO + 1
      IF(IPRO.GT.LNPRO) THEN
         WRITE(6,'(1X,A)')
     *      '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LNPRO'
         GOTO 3000
      ENDIF
*
      NPRO(IPRO,1) = IC
      NPRO(IPRO,2) = ICPR
      NPRO(IPRO,3) = MPRO(JPRO)
      NPRO(IPRO,4) = IVAR
      NPRO(IPRO,5) = IMEM
      NPRO(IPRO,6) = IDIM
      MPRO(JPRO) = IPRO
      MPLO(IPRO) = JPLO
*
      IC = NCOD(IC+2)
      GOTO 500
*
***** E P **
   8  CONTINUE
      JPRO = NPRO(IPRO,3)
      IVAR = NPRO(IPRO,4)
      IMEM = NPRO(IPRO,5)
      IDIM = NPRO(IPRO,6)
      JJPRO = NCOD(NCOD(NPRO(IPRO,1)+5)+2)
      MPRO(JJPRO) = JPRO
      IF(JPRO.NE.0) THEN
         ICO = NPRO(JPRO,1)
         ICPRO = NPRO(JPRO,2)
         JADD = NCOD(ICO+6)
         NADD(JADD) = NCOD(ICPRO+6)
         DO 1008 J=1,IARG
 1008    NADD(JADD+J) = NCOD(ICPRO+7+J)
         IVO  = NCOD(ICO+5)
         MSCR = NCOD(IVO+1)
         JADD = NCOD(IVO+3)
         JALL = NPRO(JPRO,4)
         DO 2008 J=1,MSCR
 2008    NADD(JADD+J) = JALL + J
      ENDIF
      JPLO = MPLO(IPRO)
      IC = NCOD(NPRO(IPRO,2)+2)
      IPRO = IPRO - 1
      GOTO 1000
*
***** C P **
    9 CONTINUE
      ICPR = IC
      IC = NCOD(IC+3)
      GOTO 1000
*
***** C S **
   10 CONTINUE
      IC = NCOD(IC+2)
      GOTO 1000
*
***** B I **
   11 CONTINUE
      MVAR = NADD(NCOD(IC+6))
      IF(NTYP(MVAR).NE.3) THEN
         WRITE(6,'(1X,A,I6)')
     *      '$$$ IF VARIABLE HAS WRONG TYPE ',NTYP(MVAR)
         GOTO 3000
      ENDIF
      NIF = NC(NBEG(MVAR))
      IF(NIF.EQ.1) THEN
         IC = NCOD(IC+2)
      ELSE
         IC = NCOD(IC+3)
      ENDIF
      GOTO 1000
*
***** O I **
   12 CONTINUE
      MVAR = NADD(NCOD(IC+6))
      IF(NTYP(MVAR).NE.3) THEN
         WRITE(6,'(1X,A,I6)')
     *      '$$$ ELSEIF VARIABLE HAS WRONG TYPE ',NTYP(MVAR)
         GOTO 3000
      ENDIF
      NIF = NC(NBEG(MVAR))
      IF(NIF.EQ.1) THEN
         IC = NCOD(IC+2)
      ELSE
         IC = NCOD(IC+3)
      ENDIF
      GOTO 1000
*
***** E I **
   13 CONTINUE
      IC = NCOD(IC+2)
      GOTO 1000
*
***** B W **
   14 CONTINUE
      MVAR = NADD(NCOD(IC+6))
      IF(NTYP(MVAR).NE.3) THEN
         WRITE(6,'(1X,A,I6)')
     *      '$$$ WHILE VARIABLE HAS WRONG TYPE ',NTYP(MVAR)
         GOTO 3000
      ENDIF
      NIF = NC(NBEG(MVAR))
      IF(NIF.EQ.1) THEN
         IC = NCOD(IC+2)
      ELSE
         IC = NCOD(NCOD(IC+3)+2)
      ENDIF
      GOTO 1000
*
***** E W **
   15 CONTINUE
      IC = NCOD(IC+3)
      GOTO 1000
*
***** B L **
   16 CONTINUE
      ILOO = ILOO + 1
      IF(ILOO.GT.LLOO) THEN
         WRITE(6,'(1X,A)')
     *      '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LLOO'
         GOTO 3000
      ENDIF
      ILVA = NADD(NCOD(IC+6))
      XLOO(ILOO,1) = ILVA
*
      IBEG = NADD(NCOD(IC+7))
      IF(NTYP(IBEG).NE.1) THEN
         WRITE(6,'(1X,A)') '$$$ LOOP BEGINNING NOT REAL'
         GOTO 3000
      ELSE
         XBEG = CC(NBEG(IBEG))
         XLOO(ILOO,2) = XBEG
      ENDIF
*
      IEND = NADD(NCOD(IC+8))
      IF(NTYP(IEND).NE.1) THEN
         WRITE(6,'(1X,A)') '$$$ LOOP END NOT REAL'
         GOTO 3000
      ELSE
         XEND = CC(NBEG(IEND))
         XLOO(ILOO,3) = XEND
      ENDIF
*
      IF(NCOD(IC+4)-IC.EQ.9) THEN
         ISTE = NADD(NCOD(IC+9))
         IF(NTYP(ISTE).NE.1) THEN
            WRITE(6,'(1X,A)') '$$$ LOOP STEP NOT REAL'
            GOTO 3000
         ELSE
            XSTE = CC(NBEG(ISTE))
            XLOO(ILOO,4) = XSTE
         ENDIF
      ELSEIF(NCOD(IC+4)-IC.EQ.8) THEN
         XSTE = 1.D0
         XLOO(ILOO,4) = XSTE
      ELSE
         WRITE(6,'(1X,A)') '$$$ LOOP HAS TOO MANY PARAMETERS'
         GOTO 3000
      ENDIF
*
      NTYP(ILVA) = 1
      CC(NBEG(ILVA)) = XBEG
      IF(XSTE.GE.0.D0) THEN
         IF(XBEG.GT.XEND+1.D-12) THEN
            ILOO = ILOO - 1
            IC = NCOD(NCOD(IC+3)+2)
         ELSE
            IC = NCOD(IC+2)
         ENDIF
      ELSE
         IF(XBEG.LT.XEND-1.D-12) THEN
            ILOO = ILOO - 1
            IC = NCOD(NCOD(IC+3)+2)
         ELSE
            IC = NCOD(IC+2)
         ENDIF
      ENDIF
      GOTO 1000
*
***** E L **
   17 CONTINUE
      XSTE = XLOO(ILOO,4)
      XLVA = XLOO(ILOO,2) + XSTE
      XLOO(ILOO,2) = XLVA
      XEND = XLOO(ILOO,3)
      IF(XSTE.GE.0.D0) THEN
         IF(XLVA.GT.XEND+1.D-12) THEN
            ILOO = ILOO - 1
            IC = NCOD(IC+2)
         ELSE
            ILVA = NINT(XLOO(ILOO,1))
            NTYP(ILVA) = 1
            CC(NBEG(ILVA)) = XLVA
            IC = NCOD(NCOD(IC+3)+2)
         ENDIF
      ELSE
         IF(XLVA.LT.XEND-1.D-12) THEN
            ILOO = ILOO - 1
            IC = NCOD(IC+2)
         ELSE
            ILVA = NINT(XLOO(ILOO,1))
            NTYP(ILVA) = 1
            CC(NBEG(ILVA)) = XLVA
            IC = NCOD(NCOD(IC+3)+2)
         ENDIF
      ENDIF
      GOTO 1000
*
***** B O **
   18 CONTINUE
      IF(IOPT.NE.0) THEN
         WRITE(6,'(1X,A)') '$$$ NO NESTED FITS'
         GOTO 3000
      ENDIF
      IOPT = IOPT + 1
      NOPT(IOPT,1) = NCOD(IC+4) - IC - 5
      IF(NOPT(IOPT,1).GT.LOPV-1) THEN
         WRITE(6,'(1X,A)') '$$$ TOO MANY FIT VARIABLES'
         GOTO 3000
      ENDIF
      DO 1018 J=1,NOPT(IOPT,1)
 1018 NOPT(IOPT,1+J) = NCOD(IC+5+J)
      IC = NCOD(IC+2)
      GOTO 1000
*
***** E O **
   19 CONTINUE
*
      IEPS = NADD(NCOD(IC+6))
      IF(NTYP(IEPS).NE.1) THEN
         WRITE(6,'(1X,A)') '$$$ FIT TOLERANCE NOT REAL'
         GOTO 3000
      ELSE
         XEPS = CC(NBEG(IEPS))
      ENDIF
*
      ITER = NADD(NCOD(IC+7))
      IF(NTYP(ITER).NE.1) THEN
         WRITE(6,'(1X,A)') '$$$ ITERATION LIMIT NOT REAL'
         GOTO 3000
      ELSE
         ITER = NINT(CC(NBEG(ITER)))
      ENDIF
*
      IFIT = NADD(NCOD(IC+8))
      IF(NTYP(IFIT).NE.1) THEN
         WRITE(6,'(1X,A)') '$$$ ALGORITHM IDENTIFIER NOT REAL'
         GOTO 3000
      ELSE
         IFIT = NINT(CC(NBEG(IFIT)))
      ENDIF
*
      NXOBJ = NCOD(IC+4) - IC - 8
      DO 1019 J=1,NXOBJ
      IOBJ = NADD(NCOD(IC+8+J))
      IF(NTYP(IOBJ).NE.1) THEN
         WRITE(6,'(1X,A)') '$$$ OBJECTIVE QUANTITY NOT REAL'
         GOTO 3000
      ELSE
         XOBJ(J) = CC(NBEG(IOBJ))
      ENDIF
 1019 CONTINUE
*
      JDIM = NOPT(IOPT,1)
      DO 2019 J=1,JDIM
      JVAR = NADD(NOPT(IOPT,1+J))
      IF(NTYP(JVAR).NE.1) THEN
         WRITE(6,'(1X,A)') '$$$ FIT VARIABLE NOT REAL'
         GOTO 3000
      ENDIF
      XOPT(J,IOPT) = CC(NBEG(JVAR))
      XXOPT(J) = XOPT(J,IOPT)
 2019 CONTINUE
      CALL FIT(IFIT,XXOPT,JDIM,XOBJ,NXOBJ,XEPS,ITER,IEND)
      DO 3019 J=1,JDIM
      XOPT(J,IOPT) = XXOPT(J)
      CC(NBEG(NADD(NOPT(IOPT,1+J)))) = XXOPT(J)
 3019 CONTINUE
      IF(IEND.EQ.1) THEN
         IOPT = IOPT - 1
         IC = NCOD(IC+2)
      ELSE
         IC = NCOD(NCOD(IC+3)+2)
      ENDIF
      GOTO 1000
*
***** D C **
   20 CONTINUE
      WRITE(6,'(1X,A)') 'DEBUG FOXY'
      WRITE(6,'(1X,A)') '**********'
      WRITE(6,'(1X,10I5)') (NADD(JJ),JJ=1,100)
      IC = NCOD(IC+2)
      GOTO 1000
*
***** W F **
   21 CONTINUE
      J1 = IC + 6
      J2 = NCOD(IC+4)
      INUNIT = NADD(NCOD(J1))
      IF(NTYP(INUNIT).NE.1) CALL FOXNTY(INUNIT)
      IUNI = NINT(CC(NBEG(INUNIT)))
      DO 1021 J=1,J2-J1
 1021 CALL FOXPRI(NADD(NCOD(J1+J)),IUNI,J,J2-J1)
      IC = NCOD(IC+2)
      GOTO 1000
*
***** R F **
   22 CONTINUE
      J1 = IC + 6
      J2 = NCOD(IC+4)
      INUNIT = NADD(NCOD(J1))
      IF(NTYP(INUNIT).NE.1) CALL FOXNTY(INUNIT)
      IUNI = NINT(CC(NBEG(INUNIT)))
      DO 1022 J=1,J2-J1
 1022 CALL FOXREA(NADD(NCOD(J1+J)),IUNI)
      IC = NCOD(IC+2)
      GOTO 1000
*
***** W C **
   23 CONTINUE
      IC = NCOD(IC+2)
      GOTO 1000
*
***** I C **
   24 CONTINUE
      IC = NCOD(IC+2)
      GOTO 1000
*
***** B M **
   25 CONTINUE
      IPLO = IPLO + 1
      IF(IPLO.GT.LPLO) THEN
         WRITE(6,'(1X,A)')
     *      '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LPLO'
         GOTO 3000
      ENDIF
      IPLV = NADD(NCOD(IC+6))
*
      IPBE = NADD(NCOD(IC+7))
      IF(NTYP(IPBE).NE.1) THEN
         WRITE(6,'(1X,A)') '$$$ PLOOP INDEX NOT REAL'
         GOTO 3000
      ELSE
         XPBE = CC(NBEG(IPBE))
         NPBE = INT(XPBE)
         IF(XPBE.NE.DBLE(NPBE)) THEN
            WRITE(6,'(1X,A)') '$$$ PLOOP INDEX MUST BE AN INTEGER'
            GOTO 3000
         ELSEIF(NPBE.NE.1) THEN
            WRITE(6,'(1X,A)') '$$$ PLOOP INDEX MUST START FROM 1'
            GOTO 3000
         ENDIF
      ENDIF
*
      IPEN = NADD(NCOD(IC+8))
      IF(NTYP(IPEN).NE.1) THEN
         WRITE(6,'(1X,A)') '$$$ PLOOP INDEX NOT REAL'
         GOTO 3000
      ELSE
         XPEN = CC(NBEG(IPEN))
         NPEN = INT(XPEN)
         NGRP = NGRP*NPEN
         IF(XPEN.NE.DBLE(NPEN)) THEN
            WRITE(6,'(1X,A)') '$$$ PLOOP INDEX MUST BE AN INTEGER'
            GOTO 3000
         ELSEIF(NPEN.LT.1.OR.NGRP.GT.MTASK) THEN
            WRITE(6,'(1X,A)') '$$$ PLOOP INDEX OUT OF ALLOWED RANGE'
            GOTO 3000
         ELSEIF(NPEN.GT.LPEN) THEN
            WRITE(6,'(1X,A)')
     *         '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LPEN'
            GOTO 3000
         ENDIF
      ENDIF
*
      IF(NCOD(IC+4)-IC.EQ.9) THEN
         IPST = NADD(NCOD(IC+9))
         IF(NTYP(IPST).NE.1) THEN
            WRITE(6,'(1X,A)') '$$$ PLOOP STEP NOT REAL'
            GOTO 3000
         ELSE
            XPST = CC(NBEG(IPST))
            NPST = INT(XPST)
            IF(XPST.NE.DBLE(NPST)) THEN
               WRITE(6,'(1X,A)') '$$$ PLOOP STEP MUST BE AN INTEGER'
               GOTO 3000
            ELSEIF(NPST.NE.1) THEN
               WRITE(6,'(1X,A)') '$$$ PLOOP STEP MUST BE EQUAL TO 1'
               GOTO 3000
            ENDIF
         ENDIF
      ELSEIF(NCOD(IC+4)-IC.EQ.8) THEN
         XPST = 1.D0
         NPST = 1
      ELSE
         WRITE(6,'(1X,A)') '$$$ PLOOP HAS TOO MANY PARAMETERS'
         GOTO 3000
      ENDIF
*
      GOTO (1025,2025), MEXE
      WRITE(6,'(1X,A)') '$$$ ERROR, INCORRECT PLOOP EXECUTION MODE'
      GOTO 3000
*
*     SERIAL EXECUTION
*     ****************
 1025 CONTINUE
      NMPL(1,IPLO) = IPLV
      NMPL(2,IPLO) = NPBE
      NMPL(3,IPLO) = NPEN
      NMPL(4,IPLO) = NPST
*
      NTYP(IPLV) = 1
      CC(NBEG(IPLV)) = DBLE(NPBE)
      IC = NCOD(IC+2)
      GOTO 1000
*
*     PARALLEL EXECUTION
*     ******************
 2025 CONTINUE
      DO 3025 KPLO=JPLO,1,-1
         IF(IC.EQ.NPLO(1,KPLO)) GOTO 4025
 3025 CONTINUE
      JPLO = JPLO + 1
      IF(JPLO.GT.LPLO) THEN
         WRITE(6,'(1X,A)')
     *      '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LPLO'
         GOTO 3000
      ENDIF
      NPLO(1,JPLO) = IC
      NPLO(2,JPLO) = IPLV
      NPLO(3,JPLO) = 0
      NPLO(4,JPLO) = -LPLV
      KPLO = JPLO
 4025 NNPL(1,IPLO) = KPLO
*
      DO 5025 ICOMM=JCOMM,1,-1
         IF(NGRP.EQ.NNCOM(1,ICOMM).AND.NPEN.EQ.NNCOM(2,ICOMM)) THEN
            NNPL(2,IPLO) = ICOMM
            NTYP(IPLV) = 1
            CC(NBEG(IPLV)) = DBLE(NNCOM(3,ICOMM))
            ISROOT = NNCOM(8,ICOMM)
            IC = NCOD(IC+2)
            GOTO 1000
         ENDIF
 5025 CONTINUE
      JCOMM = JCOMM + 1
      IF(JCOMM.GT.LCOMM) THEN
         WRITE(6,'(1X,A)')
     *      '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LCOMM'
         GOTO 3000
      ENDIF
      NNPL(2,IPLO) = JCOMM
      NNCOM(1,JCOMM) = NGRP
      NNCOM(2,JCOMM) = NPEN
      ICOMM = NNPL(2,IPLO-1)
      XTASK = DBLE(NNCOM(4,ICOMM))/DBLE(NPEN)
      ITASK = INT(XTASK)
      IF(XTASK.NE.DBLE(ITASK)) THEN
         WRITE(6,'(1X,A)') '$$$ INVALID PLOOP NESTING'
         GOTO 3000
      ELSE
         NNCOM(3,JCOMM) = (NNCOM(5,ICOMM)-1)/ITASK + 1
         NNCOM(4,JCOMM) = ITASK
         NNCOM(5,JCOMM) = NNCOM(5,ICOMM)-((NNCOM(3,JCOMM)-1)*ITASK)
         ICLR = NNCOM(3,JCOMM)
         KCOMM = NNCOM(7,ICOMM)
*MPI     CALL MPI_COMM_SPLIT(KCOMM,ICLR,MRANK,MCOMM,MPERR)               *MPI
         NNCOM(7,JCOMM) = MCOMM
         ICLR = NNCOM(5,JCOMM)
         IF(ICLR.NE.IROOT) ICLR = MPIUN
*MPI     CALL MPI_COMM_SPLIT(KCOMM,ICLR,MRANK,MCOMM,MPERR)               *MPI
         NNCOM(6,JCOMM) = MCOMM
         NNCOM(8,JCOMM) = 1
         IF(NNCOM(6,JCOMM).EQ.MPICN) NNCOM(8,JCOMM) = 0
      ENDIF
*
      NTYP(IPLV) = 1
      CC(NBEG(IPLV)) = DBLE(NNCOM(3,JCOMM))
      ISROOT = NNCOM(8,JCOMM)
      IC = NCOD(IC+2)
      GOTO 1000
*
***** E M **
   26 CONTINUE
      GOTO (1026,2026), MEXE
      WRITE(6,'(1X,A)') '$$$ ERROR, INCORRECT PLOOP EXECUTION MODE'
      GOTO 3000
*
*     SERIAL EXECUTION
*     ****************
 1026 CONTINUE
      NPST = NMPL(4,IPLO)
      NPLV = NMPL(2,IPLO) + NPST
      NMPL(2,IPLO) = NPLV
      NPEN = NMPL(3,IPLO)
      IF(NPLV.GT.NPEN) THEN
         IPLO = IPLO - 1
         IC = NCOD(IC+2)
      ELSE
         IPLV = NMPL(1,IPLO)
         NTYP(IPLV) = 1
         CC(NBEG(IPLV)) = DBLE(NPLV)
         IC = NCOD(NCOD(IC+3)+2)
      ENDIF
      GOTO 1000
*
*     PARALLEL EXECUTION
*     ******************
 2026 CONTINUE
      KPLO = NNPL(1,IPLO)
      ICOMM = NNPL(2,IPLO)
*
      NPEN = NNCOM(2,ICOMM)
      IRANK = NNCOM(3,ICOMM) - 1
      ITASK = NNCOM(4,ICOMM)
      KCOMM = NNCOM(6,ICOMM)
      MCOMM = NNCOM(7,ICOMM)
*
      JPLV = NPLO(3,KPLO)
      KPLV = NPLO(4,KPLO)
      MCOPT = 1
      IF(JPLV.GT.0) THEN
         IPLV = NADD(NCOD(IC+6))
         IF(NTYP(IPLV).EQ.1) THEN
            MCOPT = NINT(CC(NBEG(IPLV)))
         ELSE
            WRITE(6,'(1X,A)') '$$$ PLOOP COMMUNICATION OPTION NOT REAL'
            GOTO 3000
         ENDIF
      ENDIF
      IF(KPLV.GT.0) THEN
         GOTO 4026
      ELSEIF(KPLV.EQ.0) THEN
         GOTO 4526
      ENDIF
*
      KPLV = NCOD(IC+4) - IC - 5
      IF(KPLV.GT.0) THEN
         IPLV = NADD(NCOD(IC+6))
         IF(NTYP(IPLV).EQ.1) THEN
            MCOPT = NINT(CC(NBEG(IPLV)))
            JPLV = 1
            KPLV = KPLV - JPLV
         ELSEIF(NTYP(IPLV).GE.0) THEN
            WRITE(6,'(1X,A)') '$$$ PLOOP COMMUNICATION OPTION NOT REAL'
            GOTO 3000
         ENDIF
      ENDIF
      NPLO(3,KPLO) = JPLV
      NPLO(4,KPLO) = KPLV
*
      IF(KPLV.GT.LPLV) THEN
         WRITE(6,'(1X,A)') '$$$ TOO MANY PLOOP VARIABLES'
         GOTO 3000
      ENDIF
      JJ = IC + 5 + JPLV
      DO 3026 J=1,KPLV
         IPLV = NADD(NCOD(JJ+J))
         NPLO(4+J,KPLO) = IPLV
         IF(NTYP(IPLV).GE.0) THEN
            WRITE(6,'(1X,A)') '$$$ PLOOP VARIABLE MUST BE AN ARRAY'
            GOTO 3000
         ELSEIF(NDIM(NMAX(IPLV)-NTYP(IPLV)).LT.NPEN) THEN
            WRITE(6,'(1X,A)')
     *         '$$$ THE LAST DIMENSION OF A PLOOP VARIABLE MUST BE'
            WRITE(6,'(1X,A)')
     *         '    GREATER THAN OR EQUAL TO THE LENGTH OF PLOOP.'
            GOTO 3000
         ENDIF
 3026 CONTINUE
*
 4026 CONTINUE
      DO 5026 J=1,KPLV
         IPLV = NPLO(4+J,KPLO)
         IPDI = NMAX(IPLV)
         IPJN = NDIM(IPDI)
         IPJD = -NTYP(IPLV)
         IPJL = NEND(IPLV)
*
         IPLV = IPLV + 1
         IPBE = NBEG(IPLV)
         JCOUNT = IPJN/NDIM(IPDI+IPJD)
         ICOUNT = IPJL*JCOUNT
         IPVA = IPBE + ICOUNT*IRANK
         JPVA = IPLV + JCOUNT*IRANK
*
         GOTO (6026,6026,1426,1426,5426,5426,9426,9426), MCOPT
         WRITE(6,'(1X,A)')
     *      '$$$ INVALID PLOOP COMMUNICATION OPTION NUMBER'
         GOTO 3000
*
 6026    IF(KCOMM.NE.MPICN.AND.NPEN.GT.1) THEN
*MPI        CALL MPI_ALLGATHER(NTYP(JPVA),JCOUNT,MPI_INTEGER,NTYP(IPLV), *MPI
*MPI *         JCOUNT,MPI_INTEGER,KCOMM,MPERR)                           *MPI
*MPI        CALL MPI_ALLGATHER(NBEG(JPVA),JCOUNT,MPI_INTEGER,NBEG(IPLV), *MPI
*MPI *         JCOUNT,MPI_INTEGER,KCOMM,MPERR)                           *MPI
*MPI        CALL MPI_ALLGATHER(NEND(JPVA),JCOUNT,MPI_INTEGER,NEND(IPLV), *MPI
*MPI *         JCOUNT,MPI_INTEGER,KCOMM,MPERR)                           *MPI
*
            DO 7026 JJ=0,JCOUNT-1
               IDISP = IPLV + JJ
               JDISP = JJ*IPJL
               DO 8026 K=0,NPEN-1
                  KPVA = IDISP + K*JCOUNT
                  NRECV(K+1) = NEND(KPVA) - NBEG(KPVA) + 1
                  NDISP(K+1) = JDISP + K*ICOUNT
 8026          CONTINUE
               IPSB = IPVA + JDISP
               KPVA = JPVA + JJ
               KCOUNT = NEND(KPVA) - NBEG(KPVA) + 1
*MPI           CALL MPI_ALLGATHERV(CC(IPSB),KCOUNT,MPI_DOUBLE_PRECISION, *MPI
*MPI *            CC(IPBE),NRECV,NDISP,MPI_DOUBLE_PRECISION,KCOMM,MPERR) *MPI
*MPI           CALL MPI_ALLGATHERV(NC(IPSB),KCOUNT,MPI_INTEGER,          *MPI
*MPI *            NC(IPBE),NRECV,NDISP,MPI_INTEGER,KCOMM,MPERR)          *MPI
 7026       CONTINUE
         ENDIF
         IF(MCOPT.EQ.1.AND.ITASK.GT.1) THEN
*MPI        CALL MPI_BCAST(NTYP(IPLV),IPJN,MPI_INTEGER,                  *MPI
*MPI *         MROOT,MCOMM,MPERR)                                        *MPI
*MPI        CALL MPI_BCAST(NBEG(IPLV),IPJN,MPI_INTEGER,                  *MPI
*MPI *         MROOT,MCOMM,MPERR)                                        *MPI
*MPI        CALL MPI_BCAST(NEND(IPLV),IPJN,MPI_INTEGER,                  *MPI
*MPI *         MROOT,MCOMM,MPERR)                                        *MPI
*
            DO 9026 JJ=0,IPJN-1
               IPSB = IPBE + JJ*IPJL
               KPVA = IPLV + JJ
               KCOUNT = NEND(KPVA) - NBEG(KPVA) + 1
*MPI           CALL MPI_BCAST(CC(IPSB),KCOUNT,MPI_DOUBLE_PRECISION,      *MPI
*MPI *            MROOT,MCOMM,MPERR)                                     *MPI
*MPI           CALL MPI_BCAST(NC(IPSB),KCOUNT,MPI_INTEGER,               *MPI
*MPI *            MROOT,MCOMM,MPERR)                                     *MPI
 9026       CONTINUE
         ENDIF
         GOTO 5026
*
 1426    IF(KCOMM.NE.MPICN.AND.NPEN.GT.1) THEN
*MPI        CALL MPI_GATHER(NTYP(JPVA),JCOUNT,MPI_INTEGER,NTYP(IPLV),    *MPI
*MPI *         JCOUNT,MPI_INTEGER,MROOT,KCOMM,MPERR)                     *MPI
*MPI        CALL MPI_GATHER(NBEG(JPVA),JCOUNT,MPI_INTEGER,NBEG(IPLV),    *MPI
*MPI *         JCOUNT,MPI_INTEGER,MROOT,KCOMM,MPERR)                     *MPI
*MPI        CALL MPI_GATHER(NEND(JPVA),JCOUNT,MPI_INTEGER,NEND(IPLV),    *MPI
*MPI *         JCOUNT,MPI_INTEGER,MROOT,KCOMM,MPERR)                     *MPI
*
            DO 2426 JJ=0,JCOUNT-1
               IDISP = IPLV + JJ
               JDISP = JJ*IPJL
               DO 3426 K=0,NPEN-1
                  KPVA = IDISP + K*JCOUNT
                  NRECV(K+1) = NEND(KPVA) - NBEG(KPVA) + 1
                  NDISP(K+1) = JDISP + K*ICOUNT
 3426          CONTINUE
               IPSB = IPVA + JDISP
               KPVA = JPVA + JJ
               KCOUNT = NEND(KPVA) - NBEG(KPVA) + 1
*MPI           CALL MPI_GATHERV(CC(IPSB),KCOUNT,MPI_DOUBLE_PRECISION,    *MPI
*MPI *            CC(IPBE),NRECV,NDISP,MPI_DOUBLE_PRECISION,             *MPI
*MPI *            MROOT,KCOMM,MPERR)                                     *MPI
*MPI           CALL MPI_GATHERV(NC(IPSB),KCOUNT,MPI_INTEGER,             *MPI
*MPI *            NC(IPBE),NRECV,NDISP,MPI_INTEGER,                      *MPI
*MPI *            MROOT,KCOMM,MPERR)                                     *MPI
 2426       CONTINUE
         ENDIF
         IF(MCOPT.EQ.3.AND.IRANK.EQ.MROOT.AND.ITASK.GT.1) THEN
*MPI        CALL MPI_BCAST(NTYP(IPLV),IPJN,MPI_INTEGER,                  *MPI
*MPI *         MROOT,MCOMM,MPERR)                                        *MPI
*MPI        CALL MPI_BCAST(NBEG(IPLV),IPJN,MPI_INTEGER,                  *MPI
*MPI *         MROOT,MCOMM,MPERR)                                        *MPI
*MPI        CALL MPI_BCAST(NEND(IPLV),IPJN,MPI_INTEGER,                  *MPI
*MPI *         MROOT,MCOMM,MPERR)                                        *MPI
*
            DO 4426 JJ=0,IPJN-1
               IPSB = IPBE + JJ*IPJL
               KPVA = IPLV + JJ
               KCOUNT = NEND(KPVA) - NBEG(KPVA) + 1
*MPI           CALL MPI_BCAST(CC(IPSB),KCOUNT,MPI_DOUBLE_PRECISION,      *MPI
*MPI *            MROOT,MCOMM,MPERR)                                     *MPI
*MPI           CALL MPI_BCAST(NC(IPSB),KCOUNT,MPI_INTEGER,               *MPI
*MPI *            MROOT,MCOMM,MPERR)                                     *MPI
 4426       CONTINUE
         ENDIF
         GOTO 5026
*
 5426    IF(KCOMM.NE.MPICN.AND.NPEN.GT.1) THEN
*MPI        CALL MPI_SCATTER(NTYP(IPLV),JCOUNT,MPI_INTEGER,NTYP(JPVA),   *MPI
*MPI *         JCOUNT,MPI_INTEGER,MROOT,KCOMM,MPERR)                     *MPI
*MPI        CALL MPI_SCATTER(NBEG(IPLV),JCOUNT,MPI_INTEGER,NBEG(JPVA),   *MPI
*MPI *         JCOUNT,MPI_INTEGER,MROOT,KCOMM,MPERR)                     *MPI
*MPI        CALL MPI_SCATTER(NEND(IPLV),JCOUNT,MPI_INTEGER,NEND(JPVA),   *MPI
*MPI *         JCOUNT,MPI_INTEGER,MROOT,KCOMM,MPERR)                     *MPI
*
            DO 6426 JJ=0,JCOUNT-1
               IDISP = IPLV + JJ
               JDISP = JJ*IPJL
               DO 7426 K=0,NPEN-1
                  KPVA = IDISP + K*JCOUNT
                  NRECV(K+1) = NEND(KPVA) - NBEG(KPVA) + 1
                  NDISP(K+1) = JDISP + K*ICOUNT
 7426          CONTINUE
               IPSB = IPVA + JDISP
               KPVA = JPVA + JJ
               KCOUNT = NEND(KPVA) - NBEG(KPVA) + 1
*MPI           CALL MPI_SCATTERV(CC(IPBE),NRECV,NDISP,                   *MPI
*MPI *            MPI_DOUBLE_PRECISION,CC(IPSB),KCOUNT,                  *MPI
*MPI *            MPI_DOUBLE_PRECISION,MROOT,KCOMM,MPERR)                *MPI
*MPI           CALL MPI_SCATTERV(NC(IPBE),NRECV,NDISP,                   *MPI
*MPI *            MPI_INTEGER,NC(IPSB),KCOUNT,                           *MPI
*MPI *            MPI_INTEGER,MROOT,KCOMM,MPERR)                         *MPI
 6426       CONTINUE
         ENDIF
         IF(MCOPT.EQ.5.AND.ITASK.GT.1) THEN
*MPI        CALL MPI_BCAST(NTYP(JPVA),JCOUNT,MPI_INTEGER,MROOT,          *MPI
*MPI *         MCOMM,MPERR)                                              *MPI
*MPI        CALL MPI_BCAST(NBEG(JPVA),JCOUNT,MPI_INTEGER,MROOT,          *MPI
*MPI *         MCOMM,MPERR)                                              *MPI
*MPI        CALL MPI_BCAST(NEND(JPVA),JCOUNT,MPI_INTEGER,MROOT,          *MPI
*MPI *         MCOMM,MPERR)                                              *MPI
*
            DO 8426 JJ=0,JCOUNT-1
               IPSB = IPVA + JJ*IPJL
               KPVA = JPVA + JJ
               KCOUNT = NEND(KPVA) - NBEG(KPVA) + 1
*MPI           CALL MPI_BCAST(CC(IPSB),KCOUNT,MPI_DOUBLE_PRECISION,      *MPI
*MPI *            MROOT,MCOMM,MPERR)                                     *MPI
*MPI           CALL MPI_BCAST(NC(IPSB),KCOUNT,MPI_INTEGER,               *MPI
*MPI *            MROOT,MCOMM,MPERR)                                     *MPI
 8426       CONTINUE
         ENDIF
         GOTO 5026
*
 9426    IF(KCOMM.NE.MPICN.AND.NPEN.GT.1) THEN
*MPI        CALL MPI_BCAST(NTYP(IPLV),JCOUNT,MPI_INTEGER,MROOT,          *MPI
*MPI *         KCOMM,MPERR)                                              *MPI
*MPI        CALL MPI_BCAST(NBEG(IPLV),JCOUNT,MPI_INTEGER,MROOT,          *MPI
*MPI *         KCOMM,MPERR)                                              *MPI
*MPI        CALL MPI_BCAST(NEND(IPLV),JCOUNT,MPI_INTEGER,MROOT,          *MPI
*MPI *         KCOMM,MPERR)                                              *MPI
*
            DO 1526 JJ=0,JCOUNT-1
               JDISP = JJ*IPJL
               IPSB = IPBE + JDISP
               KPVA = IPLV + JJ
               KCOUNT = NEND(KPVA) - NBEG(KPVA) + 1
*MPI           CALL MPI_BCAST(CC(IPSB),KCOUNT,MPI_DOUBLE_PRECISION,      *MPI
*MPI *            MROOT,KCOMM,MPERR)                                     *MPI
*MPI           CALL MPI_BCAST(NC(IPSB),KCOUNT,MPI_INTEGER,               *MPI
*MPI *            MROOT,KCOMM,MPERR)                                     *MPI
*
               NTYP(JPVA+JJ) = NTYP(KPVA)
               NBEG(JPVA+JJ) = NBEG(KPVA)
               NEND(JPVA+JJ) = NEND(KPVA)
*
               KPVA = IPVA + JDISP
               DO 2526 K=0,KCOUNT-1
                  CC(KPVA+K) = CC(IPSB+K)
                  NC(KPVA+K) = NC(IPSB+K)
 2526          CONTINUE
 1526       CONTINUE
         ENDIF
         IF(MCOPT.EQ.7.AND.ITASK.GT.1) THEN
*MPI        CALL MPI_BCAST(NTYP(JPVA),JCOUNT,MPI_INTEGER,MROOT,          *MPI
*MPI *         MCOMM,MPERR)                                              *MPI
*MPI        CALL MPI_BCAST(NBEG(JPVA),JCOUNT,MPI_INTEGER,MROOT,          *MPI
*MPI *         MCOMM,MPERR)                                              *MPI
*MPI        CALL MPI_BCAST(NEND(JPVA),JCOUNT,MPI_INTEGER,MROOT,          *MPI
*MPI *         MCOMM,MPERR)                                              *MPI
*
            DO 3526 JJ=0,JCOUNT-1
               IPSB = IPVA + JJ*IPJL
               KPVA = JPVA + JJ
               KCOUNT = NEND(KPVA) - NBEG(KPVA) + 1
*MPI           CALL MPI_BCAST(CC(IPSB),KCOUNT,MPI_DOUBLE_PRECISION,      *MPI
*MPI *            MROOT,MCOMM,MPERR)                                     *MPI
*MPI           CALL MPI_BCAST(NC(IPSB),KCOUNT,MPI_INTEGER,               *MPI
*MPI *            MROOT,MCOMM,MPERR)                                     *MPI
 3526       CONTINUE
         ENDIF
 5026 CONTINUE
*
 4526 IPLO = IPLO - 1
      ICOMM = NNPL(2,IPLO)
      NGRP = NNCOM(1,ICOMM)
      ISROOT = NNCOM(8,ICOMM)
      IC = NCOD(IC+2)
      GOTO 1000
*
***** G I **
   27 CONTINUE
      J1 = IC + 6
      J2 = NCOD(IC+4)
      INUNIT = NADD(NCOD(J1))
      IF(NTYP(INUNIT).NE.1) CALL FOXNTY(INUNIT)
      IUNI = NINT(CC(NBEG(INUNIT)))
      IF(J1.EQ.J2) THEN
         CALL GSHOW(IUNI)
      ELSEIF((J1+1).EQ.J2) THEN
         CALL GREAD(IUNI,NADD(NCOD(J1+1)))
      ELSE
         WRITE(6,'(1X,A)') '$$$ TOO MANY GUIIO PARAMETERS'
         GOTO 3000
      ENDIF
      IC = NCOD(IC+2)
      GOTO 1000
*
*****
   28 GOTO 2000
***** JUMP
   29 IC = NCOD(IC+2)
      GOTO 1000
   30 GOTO 2000                                                           $$$$
*
*GEN  ARIEXE INSERT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
*     OPERATOR +
*     **********
*
   31 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1031 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   1
            CC(NBEG(INC)) = CC(NBEG(INA)) + CC(NBEG(INB))
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   4
            CALL REACM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL REAVE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   6
            CALL READA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL REACD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR +,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  3) THEN
         IF(ITB.EQ.  3) THEN
            NTYP(INC) =   3
            CALL LOALO (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR +,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  4) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   4
            CALL CMARE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   4
            CALL CMACM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   7
            CALL CMADA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL CMACD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR +,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  5) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   5
            CALL VEARE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL VEAVE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR +,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  6) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   6
            CALL DAARE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   7
            CALL DAACM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   6
            CALL DAADA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL DAACD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR +,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  7) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   7
            CALL CDARE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   7
            CALL CDACM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   7
            CALL CDADA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL CDACD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR +,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR +,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR -
*     **********
*
   32 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1032 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   1
            CC(NBEG(INC)) = CC(NBEG(INA)) - CC(NBEG(INB))
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   4
            CALL RESCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL RESVE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   6
            CALL RESDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL RESCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR -,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  4) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   4
            CALL CMSRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   4
            CALL CMSCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   7
            CALL CMSDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL CMSCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR -,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  5) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   5
            CALL VESRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL VESVE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR -,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  6) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   6
            CALL DASRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   7
            CALL DASCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   6
            CALL DASDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL DASCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR -,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  7) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   7
            CALL CDSRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   7
            CALL CDSCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   7
            CALL CDSDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL CDSCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR -,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR -,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR *
*     **********
*
   33 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1033 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   1
            CC(NBEG(INC)) = CC(NBEG(INA)) * CC(NBEG(INB))
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   4
            CALL REMCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL REMVE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   6
            CALL REMDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL REMCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR *,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  3) THEN
         IF(ITB.EQ.  3) THEN
            NTYP(INC) =   3
            CALL LOMLO (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR *,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  4) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   4
            CALL CMMRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   4
            CALL CMMCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   7
            CALL CMMDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL CMMCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR *,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  5) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   5
            CALL VEMRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL VEMVE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR *,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  6) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   6
            CALL DAMRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   7
            CALL DAMCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   6
            CALL DAMDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL DAMCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR *,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  7) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   7
            CALL CDMRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   7
            CALL CDMCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   7
            CALL CDMDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL CDMCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR *,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR *,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR /
*     **********
*
   34 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1034 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   1
            CC(NBEG(INC)) = CC(NBEG(INA)) / CC(NBEG(INB))
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   4
            CALL REDCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL REDVE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   6
            CALL REDDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL REDCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR /,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  4) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   4
            CALL CMDRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   4
            CALL CMDCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   7
            CALL CMDDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL CMDCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR /,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  5) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   5
            CALL VEDRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL VEDVE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR /,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  6) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   6
            CALL DADRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   7
            CALL DADCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   6
            CALL DADDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL DADCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR /,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  7) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   7
            CALL CDDRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  4) THEN
            NTYP(INC) =   7
            CALL CDDCM (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  6) THEN
            NTYP(INC) =   7
            CALL CDDDA (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  7) THEN
            NTYP(INC) =   7
            CALL CDDCD (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR /,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR /,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR ^
*     **********
*
   35 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1035 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   1
            CC(NBEG(INC)) = CC(NBEG(INA)) ** CC(NBEG(INB))
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR ^,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  5) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   5
            CALL VEPOW (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR ^,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR ^,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR <
*     **********
*
   36 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1036 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   3
            CALL RELRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR <,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  2) THEN
         IF(ITB.EQ.  2) THEN
            NTYP(INC) =   3
            CALL STLST (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR <,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR <,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR >
*     **********
*
   37 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1037 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   3
            CALL REGRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR >,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  2) THEN
         IF(ITB.EQ.  2) THEN
            NTYP(INC) =   3
            CALL STGST (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR >,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR >,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR =
*     **********
*
   38 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1038 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   3
            CALL REERE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR =,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  2) THEN
         IF(ITB.EQ.  2) THEN
            NTYP(INC) =   3
            CALL STEST (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR =,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR =,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR #
*     **********
*
   39 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1039 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   3
            CALL RENRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR #,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  2) THEN
         IF(ITB.EQ.  2) THEN
            NTYP(INC) =   3
            CALL STNST (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR #,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR #,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR &
*     **********
*
   40 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1040 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   5
            CALL REURE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL REUVE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR &,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  2) THEN
         IF(ITB.EQ.  2) THEN
            NTYP(INC) =   2
            CALL STUST (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR &,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  5) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   5
            CALL VEURE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL VEUVE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR &,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  8) THEN
         IF(ITB.EQ.  8) THEN
            NTYP(INC) =   8
            CALL GRUGR (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR &,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR &,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR |
*     **********
*
   41 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1041 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   1
            CALL VELGET(INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   1
            CALL VESUB (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR |,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  2) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   2
            CALL STSUB (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   2
            CALL STSUB (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR |,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  4) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   1
            CALL CMPRE (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR |,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  5) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   1
            CALL VELGET(INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   5
            CALL VESUB (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR |,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  6) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   1
            CALL DAPIC (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   1
            CALL DAPIC (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR |,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  7) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   4
            CALL DAPIC (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSEIF(ITB.EQ.  5) THEN
            NTYP(INC) =   4
            CALL DAPIC (INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR |,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR |,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     OPERATOR %
*     **********
*
   42 CONTINUE
      INA = NADD(NCOD(IC+7))
      INB = NADD(NCOD(IC+8))
      INC = NADD(NCOD(IC+6))
 1042 ITA = NTYP(INA)
      ITB = NTYP(INB)
      IF(ITA.EQ.  1) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   6
            CALL DADIDA(INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR %,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  4) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   7
            CALL DADIDA(INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR %,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  6) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   6
            CALL DADIDA(INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR %,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSEIF(ITA.EQ.  7) THEN
         IF(ITB.EQ.  1) THEN
            NTYP(INC) =   7
            CALL DADIDA(INA,INB,INC)
            IC = NCOD(IC+2)
            GOTO 1000
         ELSE
            WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR %,'//
     *      ' RIGHT ARGUMENT HAS WRONG TYPE: '//CTYID(ITB)
            GOTO 3000
         ENDIF
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN OPERATOR %,'//
     *   ' LEFT ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
   43 GOTO 2000
 1043 GOTO 2000
   44 GOTO 2000
 1044 GOTO 2000
   45 GOTO 2000
 1045 GOTO 2000
   46 GOTO 2000
 1046 GOTO 2000
   47 GOTO 2000
 1047 GOTO 2000
   48 GOTO 2000
 1048 GOTO 2000
   49 GOTO 2000
 1049 GOTO 2000
   50 GOTO 2000
 1050 GOTO 2000
*
*     INTRINSIC FUNCTION COPY
*     *************************
*
   51 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1051 ITA = NTYP(INA)
      IF(NTYP(INC).LT.0) THEN
         WRITE(6,'(1X,A)') '$$$ ERROR, COPY TARGET IS INDEXED'
         GOTO 3000
      ENDIF
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) =       (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  2) THEN
         NTYP(INC) =   2
         CALL STCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  3) THEN
         NTYP(INC) =   3
         CALL LOCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VECOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DACOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   7
         CALL CDCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  8) THEN
         NTYP(INC) =   8
         CALL GRCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN COPY, '//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION RE
*     *************************
*
   52 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1052 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL RECOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  2) THEN
         NTYP(INC) =   1
         CALL STCRE (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   1
         CALL CMRE  (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL VERE  (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   1
         CALL DACNST(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION RE    ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION ST
*     *************************
*
   53 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1053 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   2
         CALL RECSTC(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  2) THEN
         NTYP(INC) =   2
         CALL STCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  3) THEN
         NTYP(INC) =   2
         CALL LOCST (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   2
         CALL RECSTC(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION ST    ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LO
*     *************************
*
   54 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1054 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   3
         CALL LO    (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  3) THEN
         NTYP(INC) =   3
         CALL LOCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LO    ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION CM
*     *************************
*
   55 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1055 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   4
         CALL CMCMPL(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   4
         CALL CMCMPL(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   4
         CALL CDCNST(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION CM    ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION VE
*     *************************
*
   56 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1056 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL RECOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   5
         CALL CMPVE (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VECOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION VE    ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION DA
*     *************************
*
   57 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1057 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   6
         CALL DA    (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DACOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   6
         CALL CDRE  (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION DA    ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION CD
*     *************************
*
   58 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1058 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   7
         CALL RECD  (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   7
         CALL CDCMPL(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   7
         CALL CDCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION CD    ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LRE
*     *************************
*
   59 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1059 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL LRE   (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LRE   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LST
*     *************************
*
   60 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1060 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL LST   (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LST   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LLO
*     *************************
*
   61 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1061 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL LLO   (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LLO   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LCM
*     *************************
*
   62 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1062 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL LCM   (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LCM   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LVE
*     *************************
*
   63 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1063 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL LVE   (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LVE   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LDA
*     *************************
*
   64 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1064 ITA = NTYP(INA)
      IF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL LDA   (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LDA   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LCD
*     *************************
*
   65 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1065 ITA = NTYP(INA)
      IF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL LCD   (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LCD   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LGR
*     *************************
*
   66 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1066 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL LGR   (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LGR   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION TYPE
*     *************************
*
   67 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1067 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL FOXTYP(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  2) THEN
         NTYP(INC) =   1
         CALL FOXTYP(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  3) THEN
         NTYP(INC) =   1
         CALL FOXTYP(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   1
         CALL FOXTYP(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL FOXTYP(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   1
         CALL FOXTYP(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   1
         CALL FOXTYP(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  8) THEN
         NTYP(INC) =   1
         CALL FOXTYP(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION TYPE  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LENGTH
*     *************************
*
   68 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1068 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL FOXLEN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  2) THEN
         NTYP(INC) =   1
         CALL FOXLEN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  3) THEN
         NTYP(INC) =   1
         CALL FOXLEN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   1
         CALL FOXLEN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL FOXLEN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   1
         CALL FOXLEN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   1
         CALL FOXLEN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  8) THEN
         NTYP(INC) =   1
         CALL FOXLEN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LENGTH,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION VARMEM
*     *************************
*
   69 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1069 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL FOXMEM(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  2) THEN
         NTYP(INC) =   1
         CALL FOXMEM(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  3) THEN
         NTYP(INC) =   1
         CALL FOXMEM(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   1
         CALL FOXMEM(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL FOXMEM(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   1
         CALL FOXMEM(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   1
         CALL FOXMEM(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  8) THEN
         NTYP(INC) =   1
         CALL FOXMEM(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION VARMEM,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION VARPOI
*     *************************
*
   70 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1070 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL FOXPOI(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  2) THEN
         NTYP(INC) =   1
         CALL FOXPOI(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  3) THEN
         NTYP(INC) =   1
         CALL FOXPOI(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   1
         CALL FOXPOI(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL FOXPOI(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   1
         CALL FOXPOI(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   1
         CALL FOXPOI(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  8) THEN
         NTYP(INC) =   1
         CALL FOXPOI(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION VARPOI,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION EXP
*     *************************
*
   71 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1071 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = EXP   (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMEXP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VEEXP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DAEXP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION EXP   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LOG
*     *************************
*
   72 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1072 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = LOG   (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMLOG (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VELOG (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DALOG (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LOG   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION SIN
*     *************************
*
   73 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1073 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = SIN   (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMSINE(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VESINE(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DASINE(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION SIN   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION COS
*     *************************
*
   74 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1074 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = COS   (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMCOSE(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VECOSE(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DACOSE(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION COS   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION TAN
*     *************************
*
   75 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1075 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = TAN   (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VEDTAN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DADTAN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION TAN   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION ASIN
*     *************************
*
   76 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1076 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = ASIN  (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VEASIN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DAASIN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION ASIN  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION ACOS
*     *************************
*
   77 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1077 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = ACOS  (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VEACOS(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DAACOS(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION ACOS  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION ATAN
*     *************************
*
   78 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1078 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = ATAN  (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VEATAN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DAATAN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION ATAN  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION SINH
*     *************************
*
   79 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1079 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = SINH  (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMSINH(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VESINH(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DASINH(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION SINH  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION COSH
*     *************************
*
   80 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1080 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = COSH  (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMCOSH(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VECOSH(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DACOSH(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION COSH  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION TANH
*     *************************
*
   81 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1081 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = TANH  (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VETANH(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DATANH(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION TANH  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION SQRT
*     *************************
*
   82 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1082 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = SQRT  (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMSQRT(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VESQRT(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DASQRT(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION SQRT  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION ISRT
*     *************************
*
   83 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1083 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL REISRT(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VEISRT(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DAISRT(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION ISRT  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION ISRT3
*     *************************
*
   84 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1084 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL REISR3(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VEISR3(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DAISR3(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION ISRT3 ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION SQR
*     *************************
*
   85 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1085 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL RESQR (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMSQR (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VESQR (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DASQR (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   7
         CALL CDSQR (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION SQR   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION ERF
*     *************************
*
   86 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1086 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL REERF (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DAERF (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION ERF   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION WERF
*     *************************
*
   87 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1087 ITA = NTYP(INA)
      IF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMWERF(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   7
         CALL CDWERF(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION WERF  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION VMIN
*     *************************
*
   88 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1088 ITA = NTYP(INA)
      IF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL VELMIN(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION VMIN  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION VMAX
*     *************************
*
   89 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1089 ITA = NTYP(INA)
      IF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL VELMAX(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION VMAX  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION ABS
*     *************************
*
   90 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1090 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = ABS   (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   1
         CALL CMABS (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL VEABS (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   1
         CALL DANOR (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   1
         CALL CDNOR (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION ABS   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION NORM
*     *************************
*
   91 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1091 ITA = NTYP(INA)
      IF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VEABS (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   1
         CALL DANOR (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   1
         CALL CDNOR (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION NORM  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION CONS
*     *************************
*
   92 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1092 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL RECOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   1
         CALL VECNST(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   1
         CALL DACNST(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   4
         CALL CDCNST(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION CONS  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION REAL
*     *************************
*
   93 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1093 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL RECOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   1
         CALL CMRE  (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DACOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   6
         CALL CDRE  (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION REAL  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION IMAG
*     *************************
*
   94 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1094 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL REZERO(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   1
         CALL CMIM  (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DAZERO(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   6
         CALL CDIM  (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION IMAG  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION CMPLX
*     *************************
*
   95 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1095 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   4
         CALL CMCMPL(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   7
         CALL CDCMPL(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   7
         CALL CDCOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION CMPLX ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION CONJ
*     *************************
*
   96 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1096 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL RECOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  4) THEN
         NTYP(INC) =   4
         CALL CMCONJ(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  6) THEN
         NTYP(INC) =   6
         CALL DACOP (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  7) THEN
         NTYP(INC) =   7
         CALL CDCONJ(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION CONJ  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION INT
*     *************************
*
   97 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1097 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = INT   (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VEINT (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION INT   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION NINT
*     *************************
*
   98 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1098 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CC(NBEG(INC)) = NINT  (CC(NBEG(INA)))
         IC = NCOD(IC+2)
         GOTO 1000
      ELSEIF(ITA.EQ.  5) THEN
         NTYP(INC) =   5
         CALL VENINT(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION NINT  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION NOT
*     *************************
*
   99 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1099 ITA = NTYP(INA)
      IF(ITA.EQ.  3) THEN
         NTYP(INC) =   3
         CALL LONOT (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION NOT   ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION TRIM
*     *************************
*
  100 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1100 ITA = NTYP(INA)
      IF(ITA.EQ.  2) THEN
         NTYP(INC) =   2
         CALL STRIM (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION TRIM  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION LTRIM
*     *************************
*
  101 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1101 ITA = NTYP(INA)
      IF(ITA.EQ.  2) THEN
         NTYP(INC) =   2
         CALL SLTRIM(INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION LTRIM ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
*
*     INTRINSIC FUNCTION GRIU
*     *************************
*
  102 CONTINUE
      IF(NCOD(IC+4)-IC-6.NE.1) GOTO 2800
      INA = NADD(NCOD(IC+7))
      INC = NADD(NCOD(IC+6))
 1102 ITA = NTYP(INA)
      IF(ITA.EQ.  1) THEN
         NTYP(INC) =   1
         CALL GRIU  (INA,INC)
         IC = NCOD(IC+2)
         GOTO 1000
      ELSE
         WRITE(6,'(1X,A,I6)') ' $$$ ERROR IN FUNCTION GRIU  ,'//
     *   ' ARGUMENT HAS WRONG TYPE: '//CTYID(ITA)
         GOTO 3000
      ENDIF
  103 GOTO 2000
 1103 GOTO 2000
  104 GOTO 2000
 1104 GOTO 2000
  105 GOTO 2000
 1105 GOTO 2000
  106 GOTO 2000
 1106 GOTO 2000
  107 GOTO 2000
 1107 GOTO 2000
  108 GOTO 2000
 1108 GOTO 2000
  109 GOTO 2000
 1109 GOTO 2000
  110 GOTO 2000
 1110 GOTO 2000
  111 GOTO 2000
 1111 GOTO 2000
  112 GOTO 2000
 1112 GOTO 2000
  113 GOTO 2000
 1113 GOTO 2000
  114 GOTO 2000
 1114 GOTO 2000
  115 GOTO 2000
 1115 GOTO 2000
  116 GOTO 2000
 1116 GOTO 2000
  117 GOTO 2000
 1117 GOTO 2000
  118 GOTO 2000
 1118 GOTO 2000
  119 GOTO 2000
 1119 GOTO 2000
  120 GOTO 2000
 1120 GOTO 2000
  121 GOTO 2000
 1121 GOTO 2000
  122 GOTO 2000
 1122 GOTO 2000
  123 GOTO 2000
 1123 GOTO 2000
  124 GOTO 2000
 1124 GOTO 2000
  125 GOTO 2000
 1125 GOTO 2000
  126 GOTO 2000
 1126 GOTO 2000
  127 GOTO 2000
 1127 GOTO 2000
  128 GOTO 2000
 1128 GOTO 2000
  129 GOTO 2000
 1129 GOTO 2000
  130 GOTO 2000
 1130 GOTO 2000
  131 GOTO 2000
 1131 GOTO 2000
  132 GOTO 2000
 1132 GOTO 2000
  133 GOTO 2000
 1133 GOTO 2000
  134 GOTO 2000
 1134 GOTO 2000
  135 GOTO 2000
 1135 GOTO 2000
  136 GOTO 2000
 1136 GOTO 2000
  137 GOTO 2000
 1137 GOTO 2000
  138 GOTO 2000
 1138 GOTO 2000
  139 GOTO 2000
 1139 GOTO 2000
  140 GOTO 2000
 1140 GOTO 2000
  141 GOTO 2000
 1141 GOTO 2000
  142 GOTO 2000
 1142 GOTO 2000
  143 GOTO 2000
 1143 GOTO 2000
  144 GOTO 2000
 1144 GOTO 2000
  145 GOTO 2000
 1145 GOTO 2000
  146 GOTO 2000
 1146 GOTO 2000
  147 GOTO 2000
 1147 GOTO 2000
  148 GOTO 2000
 1148 GOTO 2000
  149 GOTO 2000
 1149 GOTO 2000
  150 GOTO 2000
 1150 GOTO 2000
*
*     INTRINSIC PROCEDURES
*     ********************
*
  151 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL MEMALL( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  152 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL MEMFRE( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  153 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL MEMDPV( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  154 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL MEMWRT( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  155 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL SCRLEN( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  156 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL CPUSEC( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  157 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL PWTIME( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  158 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL PNPRO ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  159 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL PROOT ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  160 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL QUIT  ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  161 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL SLEEPM( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  162 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL OS    ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  163 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL ARGGET( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  164 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL OPENF ( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  165 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL OPENFB( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  166 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL CLOSEF( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  167 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL REWF  ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  168 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL BACKF ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  169 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL READS ( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  170 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL READB ( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  171 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL WRITEB( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  172 CONTINUE
      NARG =  6
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      CALL READM ( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6)
      IC = NCOD(IC+2)
      GOTO 1000
*
  173 CONTINUE
      NARG =  6
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      CALL WRITEM( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6)
      IC = NCOD(IC+2)
      GOTO 1000
*
  174 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DAINI ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  175 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL DANOT ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  176 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL DANOTW( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  177 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL DAEPS ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  178 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL DAEPSM( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  179 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL EPSMIN( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  180 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL DAFSET( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  181 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL DAFILT( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  182 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DAPEW ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  183 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL DAREA ( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  184 CONTINUE
      NARG =  5
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      CALL DAPRV ( IA 1,IA 2,IA 3,IA 4,IA 5)
      IC = NCOD(IC+2)
      GOTO 1000
*
  185 CONTINUE
      NARG =  5
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      CALL DAREV ( IA 1,IA 2,IA 3,IA 4,IA 5)
      IC = NCOD(IC+2)
      GOTO 1000
*
  186 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DAFLO ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  187 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL CDFLO ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  188 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DAGMD ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  189 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL RERAN ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  190 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL DARAN ( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  191 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL DADIU ( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  192 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DADMU ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  193 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL DADER ( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  194 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL DAINT ( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  195 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DAPLU ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  196 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DASCL ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  197 CONTINUE
      NARG =  6
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      CALL DATRN ( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6)
      IC = NCOD(IC+2)
      GOTO 1000
*
  198 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DASGN ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  199 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL DAPEE ( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  200 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DAPEA ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  201 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL DACODE( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  202 CONTINUE
      NARG =  5
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      CALL DANORO( IA 1,IA 2,IA 3,IA 4,IA 5)
      IC = NCOD(IC+2)
      GOTO 1000
*
  203 CONTINUE
      NARG =  5
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      CALL DANORS( IA 1,IA 2,IA 3,IA 4,IA 5)
      IC = NCOD(IC+2)
      GOTO 1000
*
  204 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL DACLIW( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  205 CONTINUE
      NARG =  5
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      CALL DACQLC( IA 1,IA 2,IA 3,IA 4,IA 5)
      IC = NCOD(IC+2)
      GOTO 1000
*
  206 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DAPEP ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  207 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL DANOW ( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  208 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL DAEST ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  209 CONTINUE
      NARG =  7
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      IA 7 = NADD(NCOD(IC+12))
      CALL MTREE ( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6,IA 7)
      IC = NCOD(IC+2)
      GOTO 1000
*
  210 CONTINUE
      NARG =  5
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      CALL CDF2  ( IA 1,IA 2,IA 3,IA 4,IA 5)
      IC = NCOD(IC+2)
      GOTO 1000
*
  211 CONTINUE
      NARG =  8
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      IA 7 = NADD(NCOD(IC+12))
      IA 8 = NADD(NCOD(IC+13))
      CALL CDNF  ( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6,IA 7,IA 8)
      IC = NCOD(IC+2)
      GOTO 1000
*
  212 CONTINUE
      NARG =  7
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      IA 7 = NADD(NCOD(IC+12))
      CALL CDNFDA( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6,IA 7)
      IC = NCOD(IC+2)
      GOTO 1000
*
  213 CONTINUE
      NARG =  7
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      IA 7 = NADD(NCOD(IC+12))
      CALL CDNFDS( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6,IA 7)
      IC = NCOD(IC+2)
      GOTO 1000
*
  214 CONTINUE
      NARG =  5
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      CALL LINV  ( IA 1,IA 2,IA 3,IA 4,IA 5)
      IC = NCOD(IC+2)
      GOTO 1000
*
  215 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL LDET  ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  216 CONTINUE
      NARG =  6
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      CALL LEV   ( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6)
      IC = NCOD(IC+2)
      GOTO 1000
*
  217 CONTINUE
      NARG =  5
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      CALL MBLOCK( IA 1,IA 2,IA 3,IA 4,IA 5)
      IC = NCOD(IC+2)
      GOTO 1000
*
  218 CONTINUE
      NARG =  5
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      CALL LSLINE( IA 1,IA 2,IA 3,IA 4,IA 5)
      IC = NCOD(IC+2)
      GOTO 1000
*
  219 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL SUBSTR( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  220 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL STCRE ( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  221 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL RECST ( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  222 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL VELSET( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  223 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL VELGET( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  224 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL VEDOT ( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  225 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL VEUNIT( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  226 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL VEZERO( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  227 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL IMUNIT( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  228 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL LTRUE ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  229 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL LFALSE( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  230 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL INTPOL( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  231 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL CLEAR ( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  232 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL GRMOVE( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  233 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL GRDRAW( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  234 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL GRDOT ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  235 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL GRTRI ( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  236 CONTINUE
      NARG =  4
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      CALL GRPOLY( IA 1,IA 2,IA 3,IA 4)
      IC = NCOD(IC+2)
      GOTO 1000
*
  237 CONTINUE
      NARG = 10
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      IA 7 = NADD(NCOD(IC+12))
      IA 8 = NADD(NCOD(IC+13))
      IA 9 = NADD(NCOD(IC+14))
      IA10 = NADD(NCOD(IC+15))
      CALL GRCURV( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6,IA 7,IA 8,IA 9,IA10)
      IC = NCOD(IC+2)
      GOTO 1000
*
  238 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL GRCHAR( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  239 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL GRCOLR( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  240 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL GRWDTH( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  241 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL GRPROJ( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  242 CONTINUE
      NARG =  7
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      IA 7 = NADD(NCOD(IC+12))
      CALL GRZOOM( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6,IA 7)
      IC = NCOD(IC+2)
      GOTO 1000
*
  243 CONTINUE
      NARG =  7
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      IA 7 = NADD(NCOD(IC+12))
      CALL GRMIMA( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6,IA 7)
      IC = NCOD(IC+2)
      GOTO 1000
*
  244 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL GREPS ( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  245 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL GRSTYL( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  246 CONTINUE
      NARG =  2
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      CALL GROUTF( IA 1,IA 2)
      IC = NCOD(IC+2)
      GOTO 1000
*
  247 CONTINUE
      NARG =  3
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      CALL GUISET( IA 1,IA 2,IA 3)
      IC = NCOD(IC+2)
      GOTO 1000
*
  248 CONTINUE
      NARG =  5
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      CALL RKCO  ( IA 1,IA 2,IA 3,IA 4,IA 5)
      IC = NCOD(IC+2)
      GOTO 1000
*
  249 CONTINUE
      NARG =  1
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      CALL POLSET( IA 1)
      IC = NCOD(IC+2)
      GOTO 1000
*
  250 CONTINUE
      NARG =  7
      IF(NCOD(IC+4)-IC-5.NE.NARG) GOTO 2900
      IA 1 = NADD(NCOD(IC+ 6))
      IA 2 = NADD(NCOD(IC+ 7))
      IA 3 = NADD(NCOD(IC+ 8))
      IA 4 = NADD(NCOD(IC+ 9))
      IA 5 = NADD(NCOD(IC+10))
      IA 6 = NADD(NCOD(IC+11))
      IA 7 = NADD(NCOD(IC+12))
      CALL POLVAL( IA 1,IA 2,IA 3,IA 4,IA 5,IA 6,IA 7)
      IC = NCOD(IC+2)
      GOTO 1000
*
  251 GOTO 2000
  252 GOTO 2000
  253 GOTO 2000
  254 GOTO 2000
  255 GOTO 2000
  256 GOTO 2000
  257 GOTO 2000
  258 GOTO 2000
  259 GOTO 2000
  260 GOTO 2000
  261 GOTO 2000
  262 GOTO 2000
  263 GOTO 2000
  264 GOTO 2000
  265 GOTO 2000
  266 GOTO 2000
  267 GOTO 2000
  268 GOTO 2000
  269 GOTO 2000
  270 GOTO 2000
  271 GOTO 2000
  272 GOTO 2000
  273 GOTO 2000
  274 GOTO 2000
  275 GOTO 2000
  276 GOTO 2000
  277 GOTO 2000
  278 GOTO 2000
  279 GOTO 2000
  280 GOTO 2000
  281 GOTO 2000
  282 GOTO 2000
  283 GOTO 2000
  284 GOTO 2000
  285 GOTO 2000
  286 GOTO 2000
  287 GOTO 2000
  288 GOTO 2000
  289 GOTO 2000
  290 GOTO 2000
  291 GOTO 2000
  292 GOTO 2000
  293 GOTO 2000
  294 GOTO 2000
  295 GOTO 2000
  296 GOTO 2000
  297 GOTO 2000
  298 GOTO 2000
  299 GOTO 2000
  300 GOTO 2000
  301 GOTO 2000
  302 GOTO 2000
  303 GOTO 2000
  304 GOTO 2000
  305 GOTO 2000
  306 GOTO 2000
  307 GOTO 2000
  308 GOTO 2000
  309 GOTO 2000
  310 GOTO 2000
  311 GOTO 2000
  312 GOTO 2000
  313 GOTO 2000
  314 GOTO 2000
  315 GOTO 2000
  316 GOTO 2000
  317 GOTO 2000
  318 GOTO 2000
  319 GOTO 2000
  320 GOTO 2000
  321 GOTO 2000
  322 GOTO 2000
  323 GOTO 2000
  324 GOTO 2000
  325 GOTO 2000
  326 GOTO 2000
  327 GOTO 2000
  328 GOTO 2000
  329 GOTO 2000
  330 GOTO 2000
  331 GOTO 2000
  332 GOTO 2000
  333 GOTO 2000
  334 GOTO 2000
  335 GOTO 2000
  336 GOTO 2000
  337 GOTO 2000
  338 GOTO 2000
  339 GOTO 2000
  340 GOTO 2000
  341 GOTO 2000
  342 GOTO 2000
  343 GOTO 2000
  344 GOTO 2000
  345 GOTO 2000
  346 GOTO 2000
  347 GOTO 2000
  348 GOTO 2000
  349 GOTO 2000
  350 GOTO 2000
*GEN  ARIEXE INSERT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*
*     ARRAY LOOKUP
*     ************
*
  351 CONTINUE
      INDA = NCOD(IC+4) - IC - 7
      IARR = NCOD(IC+7)
      IF(-NTYP(NADD(IARR)).NE.INDA) THEN
         WRITE(6,'(1X,A)') '$$$ INDEXED VARIABLE IS NOT DECLARED ARRAY'
         GOTO 3000
      ENDIF
      IIDIM = NMAX(NADD(IARR))
      JNUM = NDIM(IIDIM)
      IAAD = 0
*
      DO 1351 J=INDA,2,-1
      JDIM = NADD(NCOD(IC+7+J))
      IF(NTYP(JDIM).NE.1) THEN
         WRITE(6,'(1X,A,I2,A)') '$$$ ARRAY INDEX ',J,' NOT REAL'
         GOTO 3000
      ENDIF
      IJ = NINT(CC(NBEG(JDIM)))
      NJ = NDIM(IIDIM+J)
      KJ = NDIM(IIDIM+J-1)
      IF(IJ.LT.1.OR.IJ.GT.NJ) THEN
         WRITE(6,'(1X,A,I2,A)') '$$$ ARRAY INDEX ',J,' OUT OF BOUND'
         WRITE(6,'(1X,A,I6,A,I6)') '    INDEX, BOUND = ',IJ,',',NJ
         GOTO 3000
      ENDIF
      IAAD = IAAD + IJ - 1
      IAAD = IAAD*KJ
 1351 CONTINUE
*
      JDIM = NADD(NCOD(IC+8))
      IF(NTYP(JDIM).NE.1) THEN
         WRITE(6,'(1X,A,I2,A)') '$$$ ARRAY INDEX ',1,' NOT REAL'
         GOTO 3000
      ENDIF
      IJ = NINT(CC(NBEG(JDIM)))
      NJ = NDIM(IIDIM+1)
      IF(IJ.LT.1.OR.IJ.GT.NJ) THEN
         WRITE(6,'(1X,A,I2,A)') '$$$ ARRAY INDEX ',1,' OUT OF BOUND'
         WRITE(6,'(1X,A,I6,A,I6)') '    INDEX, BOUND = ',IJ,',',NJ
         GOTO 3000
      ENDIF
      IAAD = IAAD + IJ + NADD(IARR)
*
      INC = NADD(NCOD(IC+6))
      INA = IAAD
      GOTO 1051
*
*     COPY INTO ARRAY POSITION
*     ************************
*
  352 CONTINUE
      INDA = NCOD(IC+4) - IC - 7
      IARR = NCOD(IC+6)
      IF(-NTYP(NADD(IARR)).NE.INDA) THEN
         WRITE(6,'(1X,A)') '$$$ INDEXED VARIABLE IS NOT DECLARED ARRAY'
         GOTO 3000
      ENDIF
      IIDIM = NMAX(NADD(IARR))
      JNUM = NDIM(IIDIM)
      IAAD = 0
*
      DO 1352 J=INDA,2,-1
      JDIM = NADD(NCOD(IC+8+INDA-J))
      IF(NTYP(JDIM).NE.1) THEN
         WRITE(6,'(1X,A,I2,A)') '$$$ ARRAY INDEX ',J,' NOT REAL'
         GOTO 3000
      ENDIF
      IJ = NINT(CC(NBEG(JDIM)))
      NJ = NDIM(IIDIM+J)
      KJ = NDIM(IIDIM+J-1)
      IF(IJ.LT.1.OR.IJ.GT.NJ) THEN
         WRITE(6,'(1X,A,I2,A)') '$$$ ARRAY INDEX ',J,' OUT OF BOUND'
         WRITE(6,'(1X,A,I6,A,I6)') '    INDEX, BOUND = ',IJ,',',NJ
         GOTO 3000
      ENDIF
      IAAD = IAAD + IJ - 1
      IAAD = IAAD*KJ
 1352 CONTINUE
*
      JDIM = NADD(NCOD(IC+7+INDA))
      IF(NTYP(JDIM).NE.1) THEN
         WRITE(6,'(1X,A,I2,A)') '$$$ ARRAY INDEX ',1,' NOT REAL'
         GOTO 3000
      ENDIF
      IJ = NINT(CC(NBEG(JDIM)))
      NJ = NDIM(IIDIM+1)
      IF(IJ.LT.1.OR.IJ.GT.NJ) THEN
         WRITE(6,'(1X,A,I2,A)') '$$$ ARRAY INDEX ',1,' OUT OF BOUND'
         WRITE(6,'(1X,A,I6,A,I6)') '    INDEX, BOUND = ',IJ,',',NJ
         GOTO 3000
      ENDIF
      IAAD = IAAD + IJ + NADD(IARR)
*
      INA = NADD(NCOD(IC+7))
      INC = IAAD
      GOTO 1051
*
*     FUNCTION BRANCHING
*     ******************
*
  353 CONTINUE
      ICFU = IC
      IC = NCOD(IC+3)
      GOTO 1000
*
*
************************************************************************
*
*     ALLOCATION OF SCRATCH VARIABLES
*     *******************************
*
  500 CONTINUE
*
      IF(JADD+LSCR.GT.LADD) THEN
         WRITE(6,'(1X,A)')
     *      '!!! RUNTIME MEMORY EXHAUSTION, INCREASE LADD'
         GOTO 3000
      ENDIF
      CALL VARCHK(IVAR+MSCR)
      CALL MEMCHK(IMEM+MSCR*NSCR)
*
      DO 1500 J=1,MSCR
      IVAR = IVAR + 1
      NTYP(IVAR) = 1
      NBEG(IVAR) = IMEM + 1
      NMAX(IVAR) = IMEM + NSCR
      NEND(IVAR) = IMEM + 1
      IMEM = IMEM + NSCR
      NADD(JADD+J) = IVAR
 1500 CONTINUE
      MMEM = MAX(MMEM,IMEM)
*
      GOTO 1000
*
*     CODE ERROR EXIT
*     ***************
*
 2000 CONTINUE
      WRITE(6,'(1X,A,I6)') '@@@ CODE ERROR IN CODEXE: ',NCOD(IC+1)
      CALL FOXSTP(1)
*
*     EXECUTION ERROR EXIT
*     ********************
*
 2800 CONTINUE
      WRITE(6,'(1X,A/5X,A,I3)')
     *   '$$$ WRONG NUMBER OF INTRINSIC FUNCTION ARGUMENTS',
     *   'ARGUMENTS EXPECTED:  1, SUPPLIED:',NCOD(IC+4)-IC-6
      GOTO 3000
*
 2900 CONTINUE
      WRITE(6,'(1X,A/5X,2(A,I3))')
     *   '$$$ WRONG NUMBER OF INTRINSIC PROCEDURE ARGUMENTS',
     *   'ARGUMENTS EXPECTED:',NARG,', SUPPLIED:',NCOD(IC+4)-IC-5
      GOTO 3000
*
 3000 CONTINUE
      CALL FOXSTL
      END
*
      SUBROUTINE GETCOM(A,LAIN,LA, NA,LNA,IA, IASS,IEND,LERR)
*     *******************************************************
*
*     THIS SUBROUTINE GETS THE NEXT COMMAND FROM THE INPUT FILE AND STORES
*     IT IN THE CHARACTER A.
*     IT CAPITALIZES THE LETTERS IN THE COMMAND. IF THE COMMAND CONTAINS
*     AN ASSIGNMENT (I.E. := ), IASS = 1, ELSE IASS = 0.
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
      CHARACTER A*(*),ALIN*(2*512),AA*1,AA1*1
      INTEGER NA(*)
*
      SAVE KA,ICALL,ICOM,IPAR,ICHA
      DATA KA, ICOM, IPAR, ICHA, ICALL / 0, 0, 0, 0, 0 /
*
      IF(ICALL.EQ.0) THEN
         ICALL = 1
         CALL DFCHAR
      ENDIF
*----------------------------------------------------------------------------
*
      IASS = 0
      IEND = 0
      IA   = 0
*
      DO 10 I=LA+1,KA
  10  A(I-LA:I-LA) = A(I:I)
*
      LAO = 0
      KA  = KA - LA
*
*     EXTRACTING, CAPITALIZING AND PRE-ANALYZING FOX COMMAND
*     ******************************************************
*
  20  CONTINUE
*
      AA = A(LAO+1:LAO+1)
      IF(KA.GE.LAO+1.AND.AA.NE.SCSP.AND.AA.NE.SCHT.AND.
     *   IPAR.EQ.0.AND.AA.NE.'{'.AND.ICOM.EQ.0.AND.ICHA.EQ.0) THEN
         IA = IA + 1
         NA(IA) = LAO+1
      ENDIF
*
      DO 30 LA=LAO+1,KA
*     *****************
*
      AA = A(LA:LA)
      AA1 = '*'
      IF(LA.LT.KA) AA1 = A(LA+1:LA+1)
      IF(ICOM.NE.0) THEN
         IF(AA.EQ.'{') THEN
            ICOM = ICOM + 1
         ELSEIF(AA.EQ.'}') THEN
            ICOM = ICOM - 1
            IF(ICOM.EQ.0.AND.LA.NE.KA.AND.AA1.NE.SCSP.AND.AA1.NE.SCHT
     *         .AND.AA1.NE.'{'.AND.IPAR.EQ.0) THEN
               IA = IA + 1
               IF(IA.GT.LNA) THEN
                  WRITE(6,'(1X,A)')'!!! MEMORY EXHAUSTION, INCREASE LNA'
                  CALL FOXSTP(1)
               ENDIF
               NA(IA) = LA + 1
            ENDIF
         ENDIF
         A(LA:LA) = ' '
      ELSEIF(ICHA.EQ.0) THEN
         IF(AA.EQ.'''') THEN
            ICHA = 1
         ELSEIF(AA.EQ.'{') THEN
            ICOM = 1
            A(LA:LA) = ' '
         ELSEIF(AA.EQ.'}') THEN
            WRITE(2,'(1X,A)') '### ERROR, NO MATCHING "{" FOUND FOR "}"'
            WRITE(6,'(1X,A)') '### ERROR, NO MATCHING "{" FOUND FOR "}"'
            LERR = 1
         ELSEIF(AA.EQ.'(') THEN
            IPAR = IPAR + 1
         ELSEIF(AA.EQ.')') THEN
            IPAR = IPAR - 1
         ELSEIF(AA.EQ.':')THEN
            IASS = 1
         ELSEIF(AA.EQ.SCSP.OR.AA.EQ.SCHT) THEN
            IF(AA1.NE.SCSP.AND.AA1.NE.SCHT
     *         .AND.AA1.NE.'{'.AND.IPAR.EQ.0) THEN
               IA = IA + 1
               IF(IA.GT.LNA) THEN
                  WRITE(6,'(1X,A)')'!!! MEMORY EXHAUSTION, INCREASE LNA'
                  CALL FOXSTP(1)
               ENDIF
               NA(IA) = LA + 1
            ENDIF
            A(LA:LA) = ' '
         ELSEIF(AA.EQ.';')  THEN
            IF(NA(IA).NE.LA) THEN
               IA = IA + 1
               IF(IA.GT.LNA) THEN
                  WRITE(6,'(1X,A)')'!!! MEMORY EXHAUSTION, INCREASE LNA'
                  CALL FOXSTP(1)
               ENDIF
               NA(IA) = LA
            ENDIF
            RETURN
         ELSE
            ICI = ICHAR(AA)
            IF(ICI.GE.ICSA.AND.ICI.LE.ICSZ) THEN
               ICI = ICI + IDIF
               A(LA:LA) = CHAR(ICI)
            ENDIF
         ENDIF
      ELSE
         IF(AA.EQ.'''') THEN
            IF(ICHA.EQ.1) THEN
               IF(AA1.EQ.'''') THEN
                  ICHA = 2
               ELSE
                  ICHA = 0
               ENDIF
            ELSEIF(ICHA.EQ.2) THEN
               ICHA = 1
            ENDIF
         ENDIF
      ENDIF
*
  30  CONTINUE
*
      IF(KA.NE.0) KA = ILAST(A,1,KA)
      LA = KA
*
      IF(LA+512.GE.LAIN) THEN
         WRITE(2,'(1X,A)') '### ERROR, COMMAND TOO LONG'
         WRITE(6,'(1X,A)') '### ERROR, COMMAND TOO LONG'
         LA = 0
         LERR = 1
      ENDIF
*
      READ(1,'(A)',END=100) ALIN
      IF(ILAST(ALIN,1,2*512).GT.512) THEN
         WRITE(2,'(1X,A)') '### ERROR, LINE EXCEEDS 512 COLUMNS'
         WRITE(6,'(1X,A)') '### ERROR, LINE EXCEEDS 512 COLUMNS'
         WRITE(2,'(A)') ALIN(1:ILAST(ALIN,1,2*512))
         WRITE(6,'(A)') ALIN(1:ILAST(ALIN,1,2*512))
         CALL FOXSTP(1)
      ENDIF
      A(LA+1:LA+512) = ALIN(1:512)
      LAO = LA
      KA  = ILAST(A,1,LA+512)
*
      GOTO 20
*
  100 CONTINUE
      IF(LA.NE.1) THEN
         IF(ICOM.NE.0) THEN
            WRITE(2,'(1X,A)') '### ERROR, EXPECTED CURLY BRACKET "}"'
            WRITE(6,'(1X,A)') '### ERROR, EXPECTED CURLY BRACKET "}"'
         ELSEIF(ICHA.NE.0) THEN
            WRITE(2,'(1X,A)') '### ERROR, EXPECTED END OF STRING'
            WRITE(6,'(1X,A)') '### ERROR, EXPECTED END OF STRING'
         ENDIF
*        WRITE(2,'(1X,A)') '--- CODE INCOMPLETE'
*        WRITE(6,'(1X,A)') '--- CODE INCOMPLETE'
         CALL FOXSTP(0)
      ENDIF
      IEND = 1
*
      RETURN
      END
*
      SUBROUTINE HEADER(IU,C1,C2,C3,C4)
*     *********************************
*
*     THIS SUBROUTINE PRODUCES A HEADER FILLED WITH THE INFORMATION IN
*     THE FOUR CHARACTERS
*
      CHARACTER*29 C1,C2,C3,C4
*
      WRITE(IU,'(23X,A)') ' '
      WRITE(IU,'(23X,A)') ' '
      WRITE(IU,'(23X,A)') ' '
      WRITE(IU,'(23X,A)') ' '
      WRITE(IU,'(23X,A)') '************************************'
      WRITE(IU,'(23X,A)') '*                                  *'
      WRITE(IU,'(23X,A)') '*   '            //C1//         '  *'
      WRITE(IU,'(23X,A)') '*                                  *'
      WRITE(IU,'(23X,A)') '*   '            //C2//         '  *'
      WRITE(IU,'(23X,A)') '*                                  *'
      WRITE(IU,'(23X,A)') '*   '            //C3//         '  *'
      WRITE(IU,'(23X,A)') '*                                  *'
      WRITE(IU,'(23X,A)') '*   '            //C4//         '  *'
      WRITE(IU,'(23X,A)') '*                                  *'
      WRITE(IU,'(23X,A)') '************************************'
      WRITE(IU,'(23X,A)') ' '
      WRITE(IU,'(23X,A)') ' '
      WRITE(IU,'(23X,A)') ' '
      WRITE(IU,'(23X,A)') ' '
*
      RETURN
      END
*
      SUBROUTINE INPUT
*     ****************
*
*     THIS SUBROUTINE INPUTS DIRECTIVES AND FILENAMES FROM THE USER AND
*     OPENS SOURCE AND OUTPUT FILES
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      CHARACTER A*160
*
      CALL DFCHAR
*
*     PROCESS THE INPUT .FOX FILE NAME ETC VIA COMMAND LINE
*
*     NOTE: IARGC(), GETARG ARE NOT INCLUDED IN THE FORTRAN 95 STANDARD,
*           THOUGH MOST OF FORTRAN COMPILERS SUPPORT THEM (AS OF 2010).
*           THE FORTRAN 2003 STANDARD NEITHER INCLUDE THEM, BUT IT DOES
*           INCLUDE INTRINSIC ROUTINES WITH SIMILAR FUNCTIONALITY:
*           IARGC  --> COMMAND_ARGUMENT_COUNT
*           GETARG --> GET_COMMAND_ARGUMENT
*
*      IF(IARGC().GE.1) THEN
*         CALL GETARG(1,A(1:79))
*         IF(A(1:4).EQ.'-gui') THEN                                       *NORM
*            LGUI = 1                                                     *NORM
*            IF(IARGC().GE.2) THEN                                        *NORM
*               CALL GETARG(2,A(1:79))                                    *NORM
*            ELSE                                                         *NORM
*               GOTO 5                                                    *NORM
*            ENDIF                                                        *NORM
*         ENDIF                                                           *NORM
*         IL = ILAST(A,1,79)
*         IF((A(IL-3:IL-3).EQ.'.').AND.
*     *      (A(IL-2:IL-2).EQ.'F'.OR.A(IL-2:IL-2).EQ.'f').AND.
*     *      (A(IL-1:IL-1).EQ.'O'.OR.A(IL-1:IL-1).EQ.'o').AND.
*     *      (A(IL  :IL  ).EQ.'X'.OR.A(IL  :IL  ).EQ.'x'))
*     *   A(IL-3:IL) = '    '
*         WRITE(6,'(1X,A)') A(1:79)                                       *NORM
*         GOTO 25
*      ENDIF
*
*   5  CONTINUE
*      WRITE(6,'(1X,A)') 'GIVE SOURCE FILE NAME WITHOUT EXTENSION .FOX'   *NORM
*      DO 6 J=1,79
*   6  A(J:J) = ' '
*
      OPEN(12,FILE='foxyinp.dat',STATUS='OLD',ERR=10)
      READ(12,'(A)') A(1:79)
      CLOSE(12)
*      WRITE(6,'(1X,A)') A(1:79)                                          *NORM
*MPI  GOTO 25                                                            *MPI
  10  CONTINUE
*      DO 11 J=81,159
*  11  A(J:J) = ' '
*      READ(5,'(A)',END=12) A(81:159)
*      GOTO 13
*  12  CONTINUE
*      REWIND(5)
*  13  CONTINUE
*      DO 20 I=81,159
*      IF(A(I:I).NE.' ') THEN
*         DO 15 J=1,79
*         A(J:J) = A(J+80:J+80)
*  15     CONTINUE
*         OPEN(12,FILE='foxyinp.dat',STATUS='UNKNOWN')
*         WRITE(12,'(A)') A(1:79)
*         CLOSE(12)
*         GOTO 25
*      ENDIF
  20  CONTINUE
  25  CONTINUE
*
      IL = ILAST(A,1,79)
      OPEN(1,FILE=A(1:IL)//'.fox',STATUS='OLD')
      OPEN(2,FILE=A(1:IL)//'.lis',STATUS='UNKNOWN')
*
      RETURN
      END
*
      SUBROUTINE MSG(A)
*     *****************
*
*     THIS SUBROUTINE WRITES A MESSAGE TO ALL APPROPRIATE UNITS
*
      CHARACTER A*(*)
*
      WRITE(2,'(1X,A)') A
      WRITE(6,'(1X,A)') A
*
      RETURN
      END
*
      SUBROUTINE FOXSTL
*     *****************
*
*     THIS SUBROUTINE HANDLES STOPS OF THE PROGRAM WITH LIS LINE INFORMATION
*
      COMMON /CODEX/ ILIS,NSCR
*
      IF(ILIS.LT.0) THEN
         WRITE(6,'(/1X,A,I6)') '    ERROR OCCURED IN .LIS LINE ',-ILIS
      ENDIF
      CALL FOXSTP(1)
*
      END
*
      SUBROUTINE FOXSTP(I)
*     ********************
*
*     THIS SUBROUTINE HANDLES STOPS OF THE PROGRAM
*
*MPI  INCLUDE 'mpif.h'                                                   *MPI
*MPI  IF(I.EQ.0) THEN                                                    *MPI
*MPI     CALL MPI_FINALIZE(MPERR)                                        *MPI
*MPI  ELSEIF(I.EQ.1) THEN                                                *MPI
*MPI     CALL MPI_ABORT(MPI_COMM_WORLD,MPERC,MPERR)                      *MPI
*MPI  ENDIF                                                              *MPI
      STOP
*
      END
