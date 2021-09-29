*******************************************************************************
*                                                                             *
*                                                                             *
*                              GRAPHICS PACKAGE                               *
*                                                                             *
*                           PART OF THE COSY SYSTEM                           *
*                                                                             *
*                                 VERSION 10.0                                *
*                                                                             *
*                           UPDATED IN AUGUST, 2017                           *
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
*                                                                             *
*PGP                                                                     *PGP *
*GRW                                                                     *GRW *
*AQT                                                                     *AQT *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE GRPRI(INA,IUNITI)
*     ****************************
*
*     THIS SUBROUTINE DISPLAYS A GRAPHICS OBJECT AND ACTS AS AN INTERFACE
*     TO ALL AVAILABLE GRAPHICS DRIVERS. IT CALLS THE FOLLOWING
*     GENERIC PRINT ROUTINES:
*
*     ...BEG         : BEGINS THE PICTURE
*     ...MOV(X,Y,Z)  : PERFORMS A MOVE TO GIVEN POSITION
*     ...DRA(X,Y,Z)  : DRAWS A LINE TO GIVEN POSITION
*     ...DOT(X,Y,Z)  : PERFORMS A MOVE TO GIVEN POSITION AND PRINTS A DOT.
*                      THIS IS USED BY COSY TRACKING, THUS IT HAS TO BE FAST.
*     ...TRI(X,Y,Z)  : DRAWS A TRIANGLE WITH PREVIOUS TWO POSITIONS
*     ...PLY(IA,IPST): DRAW POLYNOMIAL CURVE OR SURFACE (SEE E.G. TTYPLY FOR
*                      A SIMPLE IMPLEMENTATION USING ...DRA AND ...TRI)
*                      IA IS THE STARTING POSITION OF THE GRPOLY GRAPHICS
*                      OBJECT IN COSY MEMORY.
*                      IPST SHOWS THE STATUS OF THE XYZ POLYNOMIALS, AND IF
*                      ANY ONE OF THEM IS TM ARITHMETIC FAILURE, IPST=-1.
*     ...CHA(STR,L)  : PRINTS ASCII STRING STR WITH LENGTH L
*     ...COL(CLR)    : SETS AN RGBA COLOR IF POSSIBLE
*                      CLR(I),I=1,4 FOR R,G,B,A WITH VALUES BETWEEN 0. AND 1.
*                      SEE GRCOLR FOR THE CONVENTIONAL COSY INDEXED COLOR.
*     ...WID(W)      : SETS A WIDTH IF POSSIBLE (1:DEFAULT)
*     ...END         : CONCLUDES THE PICTURE
*
*     SPECIFIC GRAPHIC STANDARDS ARE IDENTIFIED BY THE UNIT NUMBER IUNIT,
*     CLASSIFIED BY 3 LETTERS. THE FOLLOWING GRAPHICS STANDARDS ARE AVAILABLE:
*
*     IUNIT     ID    DESCRIPTION
*
*      >0       TTY   LOW-RESOLUTION 80 BY 24 ASCII OUTPUT TO FILE IUNIT
*
*      -1->-9         STANDARD INTERACTIVE OUTPUT (PLATFORM DEPENDENT)
*                        LINUX:      PGPLOT
*                        MS WINDOWS: GRWIN
*                        MAC:        AQUATERM
*
*     -10       POS   POSTSCRIPT OUTPUT         (FILE: pic001.ps  ...)
*     -11       ASC   METAFILE OUTPUT IN ASCII  (FILE: pic001.dat ...)
*     -12       PDF   PDF OUTPUT                (FILE: pic001.pdf ...)
*     -13       SVG   SVG OUTPUT                (FILE: pic001.svg ...)
*     -14       STL   STL OUTPUT (ASCII)        (FILE: pic001.stl ...)
*     -15       STL   STL OUTPUT (BINARY)       (FILE: pic001.stl ...)
*     -20       PGP   POSTSCRIPT VIA PGPLOT     (FILE: pic001.ps  ...)
*     -22       PGP   LATEX VIA PGPLOT          (FILE: pic001.tex ...)
*
*    -101->-110 PGP   PGPLOT   (X11 AND MS WINDOWS WITH THE PGPLOT LIBRARY)
*    -111->-120 GRW   GRWIN    (MS WINDOWS)
*    -121->-130 AQT   AQUATERM (APPLE MAC OS X)
*
*    -201->-210 GUI   OUTPUT IN COSY GUI WINDOW WITH GIVEN NUMBER
*
*   -1000             DUMMY UNIT (NO OUTPUT)
*
*     TO ADD A NEW GRAPHICS DRIVER, THINK OF A THREE LETTER ID, ADD THE TEN
*     ELEMENTARY ROUTINES, AND CALL THEM AT THE APPROPRIATE PLACES.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      CHARACTER STR*1024
*
      IF(NTYP(INA).NE.NGR) CALL FOXNTY(INA)
*
      CALL GRIUF(IUNITI,IUNIT,IUNITA)
*
      CALL GRPRE(INA,NELE,NARI,0)
      I = NBEG(INA)+1
*
      IF((IUNIT.LE.-101.AND.IUNIT.GE.-110).OR.
     *      (IUNIT.EQ.-20).OR.(IUNIT.EQ.-22)) THEN
         IF(IUNITA.NE.IUNIT) THEN
            PRINT*,'$$$ WARNING IN GRPRI, PGPLOT IS NOT LINKED.'
            RETURN
         ENDIF
         CALL PGPBEG(IUNIT)
         DO 110 IG=1,NELE
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         GOTO(101,102,103,104,105,106,107,108,109,109), NCI
         GOTO 110
  101    CALL PGPMOV(CC(I),CC(I+1),CC(I+2))
         GOTO 110
  102    CALL PGPDRA(CC(I),CC(I+1),CC(I+2))
         GOTO 110
  103    CALL PGPDOT(CC(I),CC(I+1),CC(I+2))
         GOTO 110
  104    CALL PGPTRI(CC(I),CC(I+1),CC(I+2))
         GOTO 110
  105    CALL PGPPLY(I,NINT(CC(I)))
         GOTO 110
  106    IS = I
         DO 100 J=1,NLEN
         STR(J:J) = CHAR(NINT(CC(IS)))
         IS = IS+1
  100    CONTINUE
         CALL PGPCHA(STR,NLEN)
         GOTO 110
  107    CALL PGPCOL(CC(I))
         GOTO 110
  108    CALL PGPWID(CC(I))
         GOTO 110
  109    CALL GENSTL(NCI,NLEN,I)
  110    I = I+NLEN
         CALL PGPEND(IUNIT,NARI)
         RETURN
*
      ELSEIF(IUNIT.LE.-111.AND.IUNIT.GE.-120) THEN
         IF(IUNITA.NE.IUNIT) THEN
            PRINT*,'$$$ WARNING IN GRPRI, GRWIN IS NOT LINKED.'
            RETURN
         ENDIF
         CALL GRWBEG(IUNIT)
         DO 210 IG=1,NELE
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         GOTO(201,202,203,204,205,206,207,208,209,209), NCI
         GOTO 210
  201    CALL GRWMOV(CC(I),CC(I+1),CC(I+2))
         GOTO 210
  202    CALL GRWDRA(CC(I),CC(I+1),CC(I+2))
         GOTO 210
  203    CALL GRWDOT(CC(I),CC(I+1),CC(I+2))
         GOTO 210
  204    CALL GRWTRI(CC(I),CC(I+1),CC(I+2))
         GOTO 210
  205    CALL GRWPLY(I,NINT(CC(I)))
         GOTO 210
  206    IS = I
         DO 200 J=1,NLEN
         STR(J:J) = CHAR(NINT(CC(IS)))
         IS = IS+1
  200    CONTINUE
         CALL GRWCHA(STR,NLEN)
         GOTO 210
  207    CALL GRWCOL(CC(I))
         GOTO 210
  208    CALL GRWWID(CC(I))
         GOTO 210
  209    CALL GENSTL(NCI,NLEN,I)
  210    I = I+NLEN
         CALL GRWEND(NARI)
         RETURN
*
      ELSEIF(IUNIT.LE.-121.AND.IUNIT.GE.-130) THEN
         IF(IUNITA.NE.IUNIT) THEN
            PRINT*,'$$$ WARNING IN GRPRI, AQUATERM IS NOT LINKED.'
            RETURN
         ENDIF
         CALL AQTBEG(IUNIT)
         DO 310 IG=1,NELE
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         GOTO(301,302,303,304,305,306,307,308,309,309), NCI
         GOTO 310
  301    CALL AQTMOV(CC(I),CC(I+1),CC(I+2))
         GOTO 310
  302    CALL AQTDRA(CC(I),CC(I+1),CC(I+2))
         GOTO 310
  303    CALL AQTDOT(CC(I),CC(I+1),CC(I+2))
         GOTO 310
  304    CALL AQTTRI(CC(I),CC(I+1),CC(I+2))
         GOTO 310
  305    CALL AQTPLY(I,NINT(CC(I)))
         GOTO 310
  306    IS = I
         DO 300 J=1,NLEN
         STR(J:J) = CHAR(NINT(CC(IS)))
         IS = IS+1
  300    CONTINUE
         CALL AQTCHA(STR,NLEN)
         GOTO 310
  307    CALL AQTCOL(CC(I))
         GOTO 310
  308    CALL AQTWID(CC(I))
         GOTO 310
  309    CALL GENSTL(NCI,NLEN,I)
  310    I = I+NLEN
         CALL AQTEND(NARI)
         RETURN
*
      ELSEIF(IUNIT.LE.-201.AND.IUNIT.GE.-210) THEN
         CALL GUIBEG(IUNIT,NELE)
         DO 410 IG=1,NELE
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         GOTO(401,402,403,404,405,406,407,408,409,409), NCI
         GOTO 410
  401    CALL GUIMOV(CC(I),CC(I+1),CC(I+2))
         GOTO 410
  402    CALL GUIDRA(CC(I),CC(I+1),CC(I+2))
         GOTO 410
  403    CALL GUIDOT(CC(I),CC(I+1),CC(I+2))
         GOTO 410
  404    CALL GUITRI(CC(I),CC(I+1),CC(I+2))
         GOTO 410
  405    CALL GUIPLY(I,NINT(CC(I)))
         GOTO 410
  406    IS = I
         DO 400 J=1,NLEN
         STR(J:J) = CHAR(NINT(CC(IS)))
         IS = IS+1
  400    CONTINUE
         CALL GUICHA(STR,NLEN)
         GOTO 410
  407    CALL GUICOL(CC(I))
         GOTO 410
  408    CALL GUIWID(CC(I))
         GOTO 410
  409    CALL GUISTL(NCI,NLEN,I)
  410    I = I+NLEN
         CALL GUIEND(IUNIT,NARI)
         RETURN
*
      ELSEIF(IUNIT.EQ.-10) THEN
         CALL POSBEG
         DO 510 IG=1,NELE
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         GOTO(501,502,503,504,505,506,507,508,509,509), NCI
         GOTO 510
  501    CALL POSMOV(CC(I),CC(I+1),CC(I+2))
         GOTO 510
  502    CALL POSDRA(CC(I),CC(I+1),CC(I+2))
         GOTO 510
  503    CALL POSDOT(CC(I),CC(I+1),CC(I+2))
         GOTO 510
  504    CALL POSTRI(CC(I),CC(I+1),CC(I+2))
         GOTO 510
  505    CALL POSPLY(I,NINT(CC(I)))
         GOTO 510
  506    IS = I
         DO 500 J=1,NLEN
         STR(J:J) = CHAR(NINT(CC(IS)))
         IS = IS+1
  500    CONTINUE
         CALL POSCHA(STR,NLEN)
         GOTO 510
  507    CALL POSCOL(CC(I))
         GOTO 510
  508    CALL POSWID(CC(I))
         GOTO 510
  509    CALL GENSTL(NCI,NLEN,I)
  510    I = I+NLEN
         CALL POSEND(NARI)
         RETURN
*
      ELSEIF(IUNIT.EQ.-12) THEN
         CALL PDFBEG(INA)
         DO 610 IG=1,NELE
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         GOTO(601,602,603,604,605,606,607,608,609,609), NCI
         GOTO 610
  601    CALL PDFMOV(CC(I),CC(I+1),CC(I+2))
         GOTO 610
  602    CALL PDFDRA(CC(I),CC(I+1),CC(I+2))
         GOTO 610
  603    CALL PDFDOT(CC(I),CC(I+1),CC(I+2))
         GOTO 610
  604    CALL PDFTRI(CC(I),CC(I+1),CC(I+2))
         GOTO 610
  605    CALL PDFPLY(I,NINT(CC(I)))
         GOTO 610
  606    IS = I
         DO 600 J=1,NLEN
         STR(J:J) = CHAR(NINT(CC(IS)))
         IS = IS+1
  600    CONTINUE
         CALL PDFCHA(STR,NLEN)
         GOTO 610
  607    CALL PDFCOL(CC(I))
         GOTO 610
  608    CALL PDFWID(CC(I))
         GOTO 610
  609    CALL GENSTL(NCI,NLEN,I)
  610    I = I+NLEN
         CALL PDFEND(NARI)
         RETURN
*
      ELSEIF(IUNIT.EQ.-13) THEN
         CALL SVGBEG
         DO 710 IG=1,NELE
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         GOTO(701,702,703,704,705,706,707,708,709,709), NCI
         GOTO 710
  701    CALL SVGMOV(CC(I),CC(I+1),CC(I+2))
         GOTO 710
  702    CALL SVGDRA(CC(I),CC(I+1),CC(I+2))
         GOTO 710
  703    CALL SVGDOT(CC(I),CC(I+1),CC(I+2))
         GOTO 710
  704    CALL SVGTRI(CC(I),CC(I+1),CC(I+2))
         GOTO 710
  705    CALL SVGPLY(I,NINT(CC(I)))
         GOTO 710
  706    IS = I
         DO 700 J=1,NLEN
         STR(J:J) = CHAR(NINT(CC(IS)))
         IS = IS+1
  700    CONTINUE
         CALL SVGCHA(STR,NLEN)
         GOTO 710
  707    CALL SVGCOL(CC(I))
         GOTO 710
  708    CALL SVGWID(CC(I))
         GOTO 710
  709    CALL GENSTL(NCI,NLEN,I)
  710    I = I+NLEN
         CALL SVGEND(NARI)
         RETURN
*
      ELSEIF((IUNIT.EQ.-14).OR.(IUNIT.EQ.-15)) THEN
         CALL STLBEG(IUNIT)
         DO 810 IG=1,NELE
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         GOTO(801,802,803,804,805,806,807,808,809,809), NCI
         GOTO 810
  801    CALL STLMOV(CC(I),CC(I+1),CC(I+2))
         GOTO 810
  802    CALL STLDRA(CC(I),CC(I+1),CC(I+2))
         GOTO 810
  803    CALL STLDOT(CC(I),CC(I+1),CC(I+2))
         GOTO 810
  804    CALL STLTRI(CC(I),CC(I+1),CC(I+2))
         GOTO 810
  805    CALL STLPLY(I,NINT(CC(I)))
         GOTO 810
  806    IS = I
         DO 800 J=1,NLEN
         STR(J:J) = CHAR(NINT(CC(IS)))
         IS = IS+1
  800    CONTINUE
         CALL STLCHA(STR,NLEN)
         GOTO 810
  807    CALL STLCOL(CC(I))
         GOTO 810
  808    CALL STLWID(CC(I))
         GOTO 810
  809    CALL GENSTL(NCI,NLEN,I)
  810    I = I+NLEN
         CALL STLEND(NARI)
         RETURN
*
      ELSEIF(IUNIT.EQ.-11) THEN
         CALL ASCBEG
         DO 910 IG=1,NELE
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         GOTO(901,902,903,904,905,906,907,908,909,909), NCI
         GOTO 910
  901    CALL ASCMOV(CC(I),CC(I+1),CC(I+2))
         GOTO 910
  902    CALL ASCDRA(CC(I),CC(I+1),CC(I+2))
         GOTO 910
  903    CALL ASCDOT(CC(I),CC(I+1),CC(I+2))
         GOTO 910
  904    CALL ASCTRI(CC(I),CC(I+1),CC(I+2))
         GOTO 910
  905    CALL ASCPLY(I,NINT(CC(I)))
         GOTO 910
  906    IS = I
         DO 900 J=1,NLEN
         STR(J:J) = CHAR(NINT(CC(IS)))
         IS = IS+1
  900    CONTINUE
         CALL ASCCHA(STR,NLEN)
         GOTO 910
  907    CALL ASCCOL(CC(I))
         GOTO 910
  908    CALL ASCWID(CC(I))
         GOTO 910
  909    CALL ASCSTL(NCI,NLEN,I)
  910    I = I+NLEN
         CALL ASCEND
         RETURN
*
      ELSEIF(IUNIT.GT.0) THEN
         CALL TTYBEG
         DO 1010 IG=1,NELE
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         GOTO(1001,1002,1003,1004,1005,1006,1007,1008,1009,1009), NCI
         GOTO 1010
 1001    CALL TTYMOV(CC(I),CC(I+1),CC(I+2))
         GOTO 1010
 1002    CALL TTYDRA(CC(I),CC(I+1),CC(I+2))
         GOTO 1010
 1003    CALL TTYDOT(CC(I),CC(I+1),CC(I+2))
         GOTO 1010
 1004    CALL TTYTRI(CC(I),CC(I+1),CC(I+2))
         GOTO 1010
 1005    CALL TTYPLY(I,NINT(CC(I)))
         GOTO 1010
 1006    IS = I
         DO 1000 J=1,NLEN
         STR(J:J) = CHAR(NINT(CC(IS)))
         IS = IS+1
 1000    CONTINUE
         CALL TTYCHA(STR,NLEN)
         GOTO 1010
 1007    CALL TTYCOL(CC(I))
         GOTO 1010
 1008    CALL TTYWID(CC(I))
         GOTO 1010
 1009    CALL GENSTL(NCI,NLEN,I)
 1010    I = I+NLEN
         CALL TTYEND(IUNIT,NARI)
         RETURN
*
      ELSEIF(IUNIT.EQ.-1000) THEN
         RETURN
*
      ELSE
         PRINT*,'### ERROR IN GRPRI, UNKNOWN GRAPHICS UNIT: ',IUNIT
         CALL FOXSTP(1)
      ENDIF
*
      END
*
      SUBROUTINE GRIUF(IUNITI,IUNIT,IUNITA)
*     *************************************
*
*     FOR THE INCOMING GRAPHICS UNIT NUMBER IUNITI FOR GRPRI,
*     THIS SUBROUTINE RETURNS THE ALLOCATED GRAPHICS UNIT NUMBER IUNIT,
*     AND THE ACTUALLY ALLOCATED GRAPHICS UNIT NUMBER IUNITA.
*     IF NOTHING IS ACTUALLY ALLOCATED, IUNITA=0 IS RETURNED.
*
      IUNIT = IUNITI
      IF(IUNIT.LE.-1.AND.IUNIT.GE.-9) THEN
         IUNIT = IUNITI - 200
*PGP     IUNIT = IUNITI - 100                                            *PGP
*GRW     IUNIT = IUNITI - 110                                            *GRW
*AQT     IUNIT = IUNITI - 120                                            *AQT
      ENDIF
*
      IF((IUNIT.LE.-101.AND.IUNIT.GE.-110).OR.
     *      (IUNIT.EQ.-20).OR.(IUNIT.EQ.-22)) THEN
         IUNITA = 0
*PGP     IUNITA = IUNIT                                                  *PGP
      ELSEIF(IUNIT.LE.-111.AND.IUNIT.GE.-120) THEN
         IUNITA = 0
*GRW     IUNITA = IUNIT                                                  *GRW
      ELSEIF(IUNIT.LE.-121.AND.IUNIT.GE.-130) THEN
         IUNITA = 0
*AQT     IUNITA = IUNIT                                                  *AQT
      ELSEIF(IUNIT.LE.-201.AND.IUNIT.GE.-210) THEN
         IUNITA = IUNIT
      ELSEIF(IUNIT.LE.-10.AND.IUNIT.GE.-15) THEN
         IUNITA = IUNIT
      ELSEIF(IUNIT.GT.0) THEN
         IUNITA = IUNIT
      ELSE
         IUNITA = 0
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GRIU(INUNIT,INIUA)
*     *****************************
*
*     THIS SUBROUTINE INTERFACES FROM COSY PROGRAM TO SUBROUTINE GRIUF
*
*     FOR THE GRAPHICS UNIT NUMBER INUNIT,
*     THE ACTUALLY ALLOCATED GRAPHICS UNIT NUMBER IS RETURNED TO INIUA.
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
      IUNITI = NINT(CC(NBEG(INUNIT)))
*
      CALL GRIUF(IUNITI,IUNIT,IUNITA)
      NTYP(INIUA) = NRE
      CC(NBEG(INIUA)) = DBLE(IUNITA)
*
      RETURN
      END
*
      SUBROUTINE GRAPAR(KMEM,KVAR,KDIM)
*     *********************************
*
*     THIS SUBROUTINE VERIFIES THAT THE PARAMETERS IN A "MEMORY MANAGEMENT"
*     COMMON BLOCK OUTSIDE THIS PACKAGE AGREE WITH THOSE HERE
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
         PRINT*,'@@@ ERROR IN GRAPAR, WRONG MEMORY MANAGEMENT '//
     *          'PARAMETERS IN FOXGRAF '
         CALL FOXSTP(1)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GRPRE(INA,NELE,NARI,IM)
*     **********************************
*
*     THIS SUBROUTINE PREPROCESSES THE GRAPHICS OBJECT INA FOR DISPLAY.
*     NELE IS THE NUMBER OF COSY GRAPHICS PRIMITIVES CONTAINED IN INA.
*     NARI IS THE NUMBER OF TM ARITHMETIC FAILURE GRPOLY GRAPHICS OBJECTS.
*     IF IM=1, IT RETURNS WHEN THE OBJECT BOUNDS BDG(,) COMPUTATION IS DONE.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*-----PDF GRPOLY RELATED DATA -----------------------------------------------
      DOUBLE PRECISION PLYS(2,0:1),PLYDEC(2,2)
      INTEGER IPLYS,IPLYE
      COMMON /GRPDAT/ PLYS,PLYDEC,IPLYS,IPLYE
*----------------------------------------------------------------------------
*
      PARAMETER(PID=.17453292519943295474371680597869271878153085708D-1)
      PARAMETER(RATIO=1.D-5)
      DOUBLE PRECISION BDP(3,2),BDZ(3,2),BD2(2,2),BD(3,2),X(2),XG(2)
      LOGICAL LCR
*
      IZOOM = 0
      DO 10 J=1,3
      PRANG(J) = 0.D0
      DO 10 LM=1,2
      BDP(J,LM) = 0.D0
 10   CONTINUE
      EPSXYZ = 0.D0
      EPSCOL = 0.D0
      LWR = 0
*
      NELE = 0
      NCO = 0
      NPLY = 0
      IPLYE = 0
      IPB = 0
      IEB = 0
      NARI = 0
      IBEG = 0
      IBEGP = 0
*
*     INITIAL SWEEPING
*     ****************
*
      I = NBEG(INA)+1
      IAE = NEND(INA)
      NLEN = 0
 100  I = I+NLEN
      IF(I.GT.IAE) GOTO 200
      NELE = NELE+1
      NCI = NC(I)/100000000
      NLEN = MOD(NC(I),100000000)
      GOTO(110,110,110,110,120,100,100,100,100,130,140,150), NCI
      PRINT*,'@@@ CODE ERROR IN GRPRE'
      CALL FOXSTP(1)
*
*     MOVE, DRAW, DOT, TRI
*
 110  NCO = NCO+1
      IF(IBEG.EQ.0) THEN
         IBEG = 1
         DO 111 J=1,3
         DO 111 LM=1,2
         BDG(J,LM) = CC(I+J-1)
 111     CONTINUE
      ELSE
         DO 112 J=1,3
         BDG(J,1) = MIN(BDG(J,1),CC(I+J-1))
         BDG(J,2) = MAX(BDG(J,2),CC(I+J-1))
 112     CONTINUE
      ENDIF
      GOTO 100
*
*     POLY
*
 120  CONTINUE
      IF(CC(I).GE.0.D0) THEN
         NPLY = NPLY+1
         IF(NPLY.EQ.1) IPB = I
         IPLYE = I
         CALL PLYPRE(I,IOBJ,NP1,NP2,LCR,-2)
         CALL PLYBD(I,IOBJ,2,2,BD)
         DO 121 J=1,3
         BDP(J,1) = MIN(BDP(J,1),BD(J,1))
         BDP(J,2) = MAX(BDP(J,2),BD(J,2))
 121     CONTINUE
      ELSE
         NARI = NARI+1
      ENDIF
      GOTO 100
*
*     EPS FOR POLY
*
 130  IF(IEB.EQ.0) IEB = I
      GOTO 100
*
*     PROJ ANGLES
*
 140  IS = I
      DO 141 J=1,3
      PRANG(J) = CC(IS)
      IS = IS+1
 141  CONTINUE
      GOTO 100
*
*     ZOOM
*
 150  IZOOM = 1
      IS = I
      DO 151 L=1,2
      DO 151 J=1,3
      BDZ(J,L) = CC(IS)
      IS = IS+1
 151  CONTINUE
      GOTO 100
*
 200  CONTINUE
      IF(NCO+NPLY.EQ.0) THEN
         DO 201 J=1,3
         DO 201 L=1,2
         BDG(J,L) = 0.D0
 201     CONTINUE
      ELSE
         DO 202 J=1,3
         BD3(J,1) = MIN(BDG(J,1),BDP(J,1))
         BD3(J,2) = MAX(BDG(J,2),BDP(J,2))
 202     CONTINUE
      ENDIF
*
*     REFINE BOUNDS FOR POLYNOMIALS
*     *****************************
*
      IF(NPLY.GE.1) THEN
         DO 210 J=1,3
         DO 210 LM=1,2
         BDP(J,LM) = 0.D0
 210     CONTINUE
*
         IBEG = 0
         I = IPB
         IF(IEB.NE.0.AND.IEB.LT.IPB) I = IEB
 300     IF(I.GT.IPLYE) GOTO 400
         NCI = NC(I)/100000000
         NLEN = MOD(NC(I),100000000)
         IF(NCI.EQ.5) THEN
            IF(CC(I).GE.0.D0) THEN
               CALL PLYPRE(I,IOBJ,NP1,NP2,LCR,-1)
               CALL PLYBD(I,IOBJ,NP1,NP2,BD)
               DO 301 J=1,3
               BDP(J,1) = MIN(BDP(J,1),BD(J,1))
               BDP(J,2) = MAX(BDP(J,2),BD(J,2))
 301           CONTINUE
            ENDIF
         ELSEIF(NCI.EQ.10) THEN
            CALL GENSTL(NCI,NLEN,I)
         ENDIF
         I = I+NLEN
         GOTO 300
*
 400     CONTINUE
         DO 401 J=1,3
         BDG(J,1) = MIN(BDG(J,1),BDP(J,1))
         BDG(J,2) = MAX(BDG(J,2),BDP(J,2))
 401     CONTINUE
      ENDIF
*
      EPSXYZ = 0.D0
      EPSCOL = 0.D0
      LWR = 0
      IF(IM.EQ.1) RETURN
*
*     STORE ANGLES, ZOOM, BOUNDS INFO IN COMMON GRCOM
*     ***********************************************
*     BDG(J,L)  : 3D BOUNDS OF THE GRAPHICS OBJECT INA
*                 [0,0]^3 WHEN NO POSITION DETERMINING GR PRIMITIVE EXISTS.
*     BD3(J,L)  : 3D BOUNDS OF THE GRAPHICS OBJECT INA WITH ZOOMING INFO.
*                 IT IS THE ZOOMING BOX IF GRZOOM IS USED.
*                 IT EQUALS TO BDG(J,L) IF GRZOOM IS NOT USED.
*     BDG2(J,L) : 2D PROJECTED BOUNDS OF BDG(J,L)
*     BD2(J,L)  : 2D PROJECTED BOUNDS OF BD3(J,L)
*     BD2F(J,L) : THE PICTURE FRAME BY INFLATING BD2(J,L) BY +-4%.   HOWEVER,
*                 WHEN ZOOMING, THE SIZE IS SET AT LEAST RATIO*{SIZE OF BDG2}.
*                 WHEN SIZE IS ZERO, THE WIDTH 2 IS ENFORCED.
*
      PHI   = PRANG(1)*PID
      THETA = PRANG(2)*PID
      SP = SIN(PHI)
      CP = COS(PHI)
      ST = SIN(THETA)
      CT = COS(THETA)
*
      IF(IZOOM.EQ.0) THEN
         DO 510 J=1,3
         DO 510 L=1,2
         BD3(J,L) = BDG(J,L)
 510     CONTINUE
      ELSE
         DO 520 J=1,3
         DO 520 L=1,2
         BD3(J,L) = BDZ(J,L)
 520     CONTINUE
      ENDIF
*
*     PROJECTED 2D BOUNDS AND FRAME
*
      IBEG = 0
      DO 550 K=1,2
      DO 550 L=1,2
      DO 550 M=1,2
      XB = BD3(1,K)
      YB = BD3(2,L)
      ZB = BD3(3,M)
      X(1) = CP*XB+SP*YB
      X(2) = CT*(-SP*XB+CP*YB)+ST*ZB
      XB = BDG(1,K)
      YB = BDG(2,L)
      ZB = BDG(3,M)
      XG(1) = CP*XB+SP*YB
      XG(2) = CT*(-SP*XB+CP*YB)+ST*ZB
      IF(IBEG.EQ.0) THEN
         IBEG = 1
         DO 530 J=1,2
         DO 530 LM=1,2
         BD2(J,LM) = X(J)
         BDG2(J,LM) = XG(J)
 530     CONTINUE
      ELSE
         DO 540 J=1,2
         BD2(J,1) = MIN(BD2(J,1),X(J))
         BD2(J,2) = MAX(BD2(J,2),X(J))
         BDG2(J,1) = MIN(BDG2(J,1),XG(J))
         BDG2(J,2) = MAX(BDG2(J,2),XG(J))
 540     CONTINUE
      ENDIF
 550  CONTINUE
*
      DO 560 J=1,2
      X(J) = BD2(J,2)-BD2(J,1)
      XG(J) = BDG2(J,2)-BDG2(J,1)
      IF(XG(J).EQ.0.D0) XG(J) = 2.D0
      IF(X(J).LT.RATIO*XG(J)) THEN
         IF(X(J).EQ.0D0) THEN
            X(J) = XG(J)
         ELSE
            X(J) = RATIO*XG(J)
         ENDIF
         BD2M = 0.5D0*(BD2(J,1)+BD2(J,2))
         BD2(J,1) = BD2M-0.5D0*X(J)
         BD2(J,2) = BD2M+0.5D0*X(J)
      ENDIF
      BD2F(J,1) = BD2(J,1)-0.04D0*X(J)
      BD2F(J,2) = BD2(J,2)+0.04D0*X(J)
 560  CONTINUE
*
      RETURN
      END
*
      SUBROUTINE PROJEC(X,Y,Z,XP,YP)
*     ******************************
*
*     THIS SUBROUTINE PERFORMS A PROJECTION FROM 3D COORDINATES INTO
*     2D COORDINATES AS SEEN FROM SPHERICAL COORDINATES PHI, THETA
*     AND SCALES THE RESULT INTO COSY'S 2D COORDINATES [0,1]^2.
*     NOTE: WHEN ZOOMING IS ENABLED, THE COORDINATES MAY LIE OUTSIDE
*     [0,1]^2. IT IS THE GRAPHICS DRIVER'S RESPONSIBILITY TO CLIP
*     THE DRAWING REGION ACCORDINGLY.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*
      XP = CP*X+SP*Y
      YP = CT*(-SP*X+CP*Y)+ST*Z
*
      XP = (XP-BD2F(1,1))/(BD2F(1,2)-BD2F(1,1))
      YP = (YP-BD2F(2,1))/(BD2F(2,2)-BD2F(2,1))
*
      RETURN
      END
*
      SUBROUTINE PROJPNT(NP)
*     **********************
*
*     THIS SUBROUTINE PERFORMS 3D to 2D COORDINATE PROJECTIONS FOR
*     THE GRPOLY DISCRETIZED POINTS PNT(J,I), I=0,..,NP, J=1,2,3 FOR X,Y,Z.
*
*     CAUTION: THE ORIGINAL DATA PNT(J,I) IS OVERWRITTEN.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      DO 10 I=0,NP
      CALL PROJEC(PNT(1,I),PNT(2,I),PNT(3,I),X,Y)
      PNT(1,I) = X
      PNT(2,I) = Y
 10   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE PLYLEN(IAG,IPLEN,IPRER)
*     **********************************
*
*     THIS SUBROUTINE PROVIDES THE LENGTH INFORMATION OF POLYNOMIALS OF THE
*     GRPOLY GRAPHICS OBJECT WITH THE STARTING POSITION IAG IN COSY MEMORY.
*
*     IPLEN : LENGTH OF EACH OF THE SEVEN POLYNOMIALS
*     IPRER : FLAG IF THE REMAINDER ERROR IS GIVEN (1: YES, 0: NO)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      INTEGER IPLEN(7),IPRER(7)
*
      DO 10 J=1,7
      IPLEN(J) = 0
      IPRER(J) = 0
 10   CONTINUE
*
      IF(CC(IAG).EQ.-1.D0) RETURN
*
*     SPACIAL POLYNOMIALS
      C = CC(IAG)
      DO 20 J=1,3
      B = AINT(C/1.D3)
      IPLEN(4-J) = INT(C-B*1.D3)
      C = B
 20   CONTINUE
*
*     COLOR POLYNOMIALS
      C = CC(IAG+1)
      DO 30 J=1,4
      B = AINT(C/1.D3)
      IPLEN(8-J) = INT(C-B*1.D3)
      C = B
 30   CONTINUE
*
      IPOS = IAG+2
      DO 40 J=1,7
      IPOS = IPOS+IPLEN(J)
      IF(IPLEN(J).GE.1) THEN
         IF(NC(IPOS-1).EQ.-1) IPRER(J) = 1
      ENDIF
 40   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE PLYPRE(IAG,IOBJ,NP1,NP2,LCR,ID)
*     ******************************************
*
*     THIS SUBROUTINE PREPARES SOME PROPERTIES OF THE GRPOLY GRAPHICS OBJECT
*     WITH THE STARTING POSITION IAG IN COSY MEMORY.
*     IOBJ IDENTIFIES -- 1: CURVE, 2: SURFACE.
*     ID CONTROLS HOW TO OBTAIN THE DISCRETIZATION NUMBERS NP1 AND NP2.
*     EXCEPT FOR ID=-2, NON-ZERO (NP1,NP2) PAIR IS RETURNED; NP2==0 FOR A CURVE.
*     SEE GRPOLY FOR THE USER INPUT, AND SEE PLYDSC FOR AUTOMATIC ESTIMATES.
*     LCR IDENTIFIES IF ANY COLOR POLYNOMIAL IS GIVEN.  .TRUE.: GIVEN.
*     FOR LCR, THE REMAINDER ERROR IS NOT CONSIDERED FOR THE JUDGEMENT.
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
*-----GRPOLY GRAPHICS OBJECT PARAMETERS  ------------------------------------
      INTEGER IPLEN(7),IPRER(7),MP1(7),MP2(7)
      COMMON /GRPCOM/ IPLEN,IPRER,MP1,MP2
*----------------------------------------------------------------------------
*
      LOGICAL LCR
*
      CALL PLYLEN(IAG,IPLEN,IPRER)
      NCP = NC(IAG+1)
      IOBJ = ABS(NCP)/1000000
*
      IF(ID.GE.-1) THEN
         INP = ABS(NCP)-IOBJ*1000000
         IF(INP.EQ.0) THEN
            NP1 = 0
            NP2 = 0
            CALL PLYDSC(IAG,1,NP1)
            IF(IOBJ.EQ.2) CALL PLYDSC(IAG,2,NP2)
            IF(ID.EQ.0.OR.(ID.EQ.-1.AND.IZOOM.EQ.0))
     *         NC(IAG+1) = -(IOBJ*1000000+NP1*1000+NP2)
         ELSE
            NP1 = INP/1000
            NP2 = MOD(INP,1000)
         ENDIF
      ELSEIF(ID.EQ.-2) THEN
         IF(NCP.LT.0) NC(IAG+1) = IOBJ*1000000
      ENDIF
*
      DO 10 J=1,7
      MP1(J) = 0
      MP2(J) = 0
 10   CONTINUE
*
      IPOS = IAG+2
      IF(IOBJ.EQ.1) THEN
         DO 30 J=1,7
         DO 20 I=IPOS,IPOS+IPLEN(J)-IPRER(J)-1
         MP1(J) = MAX(MP1(J),NC(I))
 20      CONTINUE
         IPOS = IPOS+IPLEN(J)
 30      CONTINUE
      ELSE
         DO 50 J=1,7
         DO 40 I=IPOS,IPOS+IPLEN(J)-IPRER(J)-1
         IP2 = INT(NC(I)/100)
         IP1 = NC(I)-IP2*100
         MP1(J) = MAX(MP1(J),IP1)
         MP2(J) = MAX(MP2(J),IP2)
 40      CONTINUE
         IPOS = IPOS+IPLEN(J)
 50      CONTINUE
      ENDIF
*
      LCR = .FALSE.
      DO 60 J=4,7
      IF((IPLEN(J)-IPRER(J)).GE.1) LCR = .TRUE.
 60   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE PLYDSC(IAG,IV,NP)
*     ****************************
*
*     THIS SUBROUTINE DETERMINES A DISCRETIZATION NUMBER NP ALONG THE PARAMETER
*     VARIABLE IV FOR THE GRPOLY GRAPHICS OBJECT STARTING AT IAG IN COSY MEMORY.
*     THE RESULTING NP IS AT LEAST 1.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRPOLY GRAPHICS OBJECT PARAMETERS  ------------------------------------
      INTEGER IPLEN(7),IPRER(7),MP1(7),MP2(7)
      COMMON /GRPCOM/ IPLEN,IPRER,MP1,MP2
*----------------------------------------------------------------------------
*
      PARAMETER(ERTOL=0.005D0,NGP=999)
*
      NP = 0
      IPOS = IAG+2
      DO 10 J=1,7
      IF(IPLEN(J).EQ.0) GOTO 10
*
      IF(J.LE.3) THEN
         ER = EPSXYZ
         IF(ER.EQ.0.D0) ER = ERTOL*(BD3(J,2)-BD3(J,1))
         IF(ER.EQ.0.D0) GOTO 10
      ELSE
         ER = EPSCOL
         IF(ER.EQ.0.D0) ER = ERTOL
      ENDIF
*
      CALL PLYDDB(IPOS,IPLEN(J)-IPRER(J),IV,VB)
      DP = SQRT(VB/ER)
*
C     FOR THE 1ST ORDER METHOD, USE THE FOLLOWING INSTEAD.
C     HOWEVER, ERTOL NEEDS TO BE INCREASED A LOT LIKE 0.02D0.
C     CALL PLYDBW(IPOS,IPLEN(J)-IPRER(J),IV,WB)
C     DP = WB/ER
*
      NP = MAX(NP,INT(2.D0*DP)+1)
      IPOS = IPOS+IPLEN(J)
 10   CONTINUE
*
      NP = MIN(NGP,MAX(1,NP))
*
      RETURN
      END
*
      SUBROUTINE PLYEV(IAG,ND,P2,NP1)
*     *******************************
*
*     THIS SUBROUTINE COMPUTES THE ND-DIM COORDINATES OF DISCRETIZED POINTS
*     PNT(J,I) WITH RESPECT TO THE PARAMETER P1 IN [-1,1], I=0,..,NP1,
*     OF THE GRPOLY GRAPHICS OBJECT (STARTING AT IAG IN COSY MEMORY)
*     WHEN THE SECOND PARAMETER VALUE IS P2.  J=1,..,ND.
*     FOR THE COLOR PART (J=4,5,6,7), ANY OUT OF BOUND VALUE IS CORRECTED
*     TO STAY INSIDE THE BOUND [0,1].
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
*-----GRPOLY GRAPHICS OBJECT PARAMETERS  ------------------------------------
      INTEGER IPLEN(7),IPRER(7),MP1(7),MP2(7)
      COMMON /GRPCOM/ IPLEN,IPRER,MP1,MP2
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CP1(7,0:LNO),POW(0:LNO)
      CHARACTER RGBA*4
      DATA RGBA / 'RGBA' /
*
      MP1M = 0
      MP2M = 0
      DO 10 J=1,ND
      MP1M = MAX(MP1M,MP1(J))
      MP2M = MAX(MP2M,MP2(J))
      DO 10 M=0,MP1(J)
      CP1(J,M) = 0.D0
 10   CONTINUE
*
      IPOS = IAG+2
*
      IF(MP2M.EQ.0) THEN
         DO 30 J=1,ND
         DO 20 I=IPOS,IPOS+IPLEN(J)-IPRER(J)-1
         CP1(J,NC(I)) = CC(I)
 20      CONTINUE
         IPOS = IPOS+IPLEN(J)
 30      CONTINUE
      ELSE
         POW(0) = 1.D0
         DO 40 M=1,MP2M
         POW(M) = POW(M-1)*P2
 40      CONTINUE
*
         DO 60 J=1,ND
         DO 50 I=IPOS,IPOS+IPLEN(J)-IPRER(J)-1
         IP2 = INT(NC(I)/100)
         IP1 = NC(I)-IP2*100
         CP1(J,IP1) = CP1(J,IP1)+CC(I)*POW(IP2)
 50      CONTINUE
         IPOS = IPOS+IPLEN(J)
 60      CONTINUE
      ENDIF
*
      DO 80 I=0,NP1
      P1 = (2.D0*I)/NP1-1.D0
      POW(1) = P1
      DO 70 M=2,MP1M
      POW(M) = POW(M-1)*P1
 70   CONTINUE
*
      DO 80 J=1,ND
      PNT(J,I) = CP1(J,0)
      DO 80 M=1,MP1(J)
      PNT(J,I) = PNT(J,I)+CP1(J,M)*POW(M)
 80   CONTINUE
*
      IF(ND.GE.4) THEN
         DO 100 J=4,ND
         IWC = 0
         DO 90 I=0,NP1
         IF(PNT(J,I).LT.0.D0) THEN
            IWC = IWC+1
            PNT(J,I) = 0.D0
         ELSEIF(PNT(J,I).GT.1.D0) THEN
            IWC = IWC+1
            PNT(J,I) = 1.D0
         ENDIF
 90      CONTINUE
         IF(IWC.NE.0) THEN
            PRINT*,'$$$ WARNING IN PLYEV, COLOR VALUE OF THE '//
     *             RGBA(J-3:J-3)//' POLYNOMIAL OUT OF BOUND'
         ENDIF
 100  CONTINUE
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE PLYPOS(IAG,XC)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE LAST POSITION XC(J), J=1,2,3 FOR X,Y,Z,
*     OF THE GRPOLY GRAPHICS OBJECT (STARTING AT IAG IN COSY MEMORY)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLE PRECISION XC(3)
      INTEGER IPLEN(7),IPRER(7)
*
      CALL PLYLEN(IAG,IPLEN,IPRER)
      IS = IAG+2
      DO 10 J=1,3
      XC(J) = 0.D0
      DO 10 K=1,IPLEN(J)
      IF(NC(IS).GE.0) XC(J) = XC(J)+CC(IS)
      IS = IS+1
 10   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE PLYBD(IAG,IOBJ,NP1,NP2,BD)
*     *************************************
*
*     THIS SUBROUTINE COMPUTES THE 3D BOUNDS [BD(J,1),BD(J,2)] OF THE GRPOLY
*     GRAPHICS OBJECT WITH THE STARTING POSITION IAG IN COSY MEMORY
*     USING NP1 X NP2 RASTERING POINTS. J=1,2,3 IS FOR X,Y,Z.
*     IOBJ IDENTIFIES -- 1: CURVE, 2: SURFACE.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION BD(3,2)
*
      IBEG = 0
*
      IF(IOBJ.EQ.1) THEN
         CALL PLYEV(IAG,3,0.D0,NP1)
         DO 30 I1=0,NP1
         IF(IBEG.EQ.0) THEN
            IBEG = 1
            DO 10 J=1,3
            DO 10 LM=1,2
            BD(J,LM) = PNT(J,I1)
 10         CONTINUE
         ELSE
            DO 20 J=1,3
            BD(J,1) = MIN(BD(J,1),PNT(J,I1))
            BD(J,2) = MAX(BD(J,2),PNT(J,I1))
 20         CONTINUE
         ENDIF
 30      CONTINUE
      ELSEIF(IOBJ.EQ.2) THEN
         DO 60 I2=0,NP2
         CALL PLYEV(IAG,3,(2.D0*I2)/NP2-1.D0,NP1)
         DO 60 I1=0,NP1
         IF(IBEG.EQ.0) THEN
            IBEG = 1
            DO 40 J=1,3
            DO 40 LM=1,2
            BD(J,LM) = PNT(J,I1)
 40         CONTINUE
         ELSE
            DO 50 J=1,3
            BD(J,1) = MIN(BD(J,1),PNT(J,I1))
            BD(J,2) = MAX(BD(J,2),PNT(J,I1))
 50         CONTINUE
         ENDIF
 60      CONTINUE
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE PLYBND(IA,L,B)
*     *************************
*
*     THIS SUBROUTINE COMPUTES THE BOUND OF ONE CURVE OR SURFACE POLYNOMIAL
*     NOTE: THE POLYNOMIAL BOUNDING PERFORMED HERE IS VERY ROUGH.
*
*     IA : STARTING POSITION OF THE POLYNOMIAL COEFFICIENTS IN COSY MEMORY
*     L  : NUMBER OF COEFFICIENTS IN POLYNOMIAL
*     B  : THE UPPER AND LOWER BOUNDS
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLE PRECISION B(2)
*
      DO 10 I=1,2
      B(I) = 0.D0
 10   CONTINUE
      IF(L.LE.0) RETURN
*
      DO 20 I=IA,IA+L-1
      N = NC(I)
      C = CC(I)
      IF(N.EQ.0) THEN
         B(1) = B(1)+C
         B(2) = B(2)+C
      ELSEIF(N.EQ.-1) THEN
         B(1) = B(1)-C
         B(2) = B(2)+C
      ELSEIF(MOD(N,2).EQ.0.AND.MOD(N/100,2).EQ.0) THEN
         IF(C.GT.0.D0) THEN
            B(2) = B(2)+C
         ELSE
            B(1) = B(1)+C
         ENDIF
      ELSE
         B(1) = B(1)-ABS(C)
         B(2) = B(2)+ABS(C)
      ENDIF
 20   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE PLYDER(IA,LA,IV,IB,LB)
*     *********************************
*
*     THIS SUBROUTINE COMPUTES THE DERIVATIVE OF A POLYNOMIAL FOR GRPOLY
*     STARTING AT THE COSY MEMORY ADDRESS IA WITH LENGTH LA
*     WITH RESPECT TO VARIABLE IV (IV IS LIMITED TO 1 OR 2),
*     RETURNING THE RESULTING POLYNOMIAL STARTING AT IB WITH LENGTH LB.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      INTEGER IP(2)
*
      IS = IB
      DO 20 I=IA,IA+LA-1
      IP(2) = INT(NC(I)/100)
      IP(1) = NC(I)-IP(2)*100
      IFAC = IP(IV)
      DO 10 J=1,2
      IF(IV.EQ.J) THEN
         IP(J) = IP(J)-1
         IF(IP(J).LT.0) GOTO 20
      ENDIF
 10   CONTINUE
      NC(IS) = IP(1)+100*IP(2)
      CC(IS) = CC(I)*IFAC
      IS = IS+1
 20   CONTINUE
*
      LB = IS-IB
*
      RETURN
      END
*
      SUBROUTINE PLYDBW(IA,L,IV,WB)
*     *****************************
*
*     THIS SUBROUTINE COMPUTES WB, A WIDTH OF BOUNDS OF THE DERIVATIVE OF
*     A GRPOLY POLYNOMIAL STARTING AT THE COSY MEMORY ADDRESS IA WITH LENGTH L
*     WITH RESPECT TO VARIABLE IV (IV IS LIMITED TO 1 OR 2).
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLE PRECISION B(2)
*
      WB = 0.D0
      IF(L.LE.0) RETURN
*
      CALL FOXALL(ISC,1,L)
      CALL PLYDER(IA,L,IV,NBEG(ISC),LD)
      CALL PLYBND(NBEG(ISC),LD,B)
      WB = B(2)-B(1)
      CALL FOXDAL(ISC,1)
*
      RETURN
      END
*
      SUBROUTINE PLYDDB(IA,L,IV,VB)
*     *****************************
*
*     THIS SUBROUTINE COMPUTES VB, AN ABSOLUTE VALUE OF BOUNDS OF
*     THE 2ND DERIVATIVE OF A GRPOLY POLYNOMIAL
*     STARTING AT THE COSY MEMORY ADDRESS IA WITH LENGTH L
*     WITH RESPECT TO VARIABLE IV (IV IS LIMITED TO 1 OR 2).
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DOUBLE PRECISION B(2)
*
      VB = 0.D0
      IF(L.LE.0) RETURN
*
      CALL FOXALL(ISC,1,L)
      CALL PLYDER(IA,L,IV,NBEG(ISC),LD)
      CALL PLYDER(NBEG(ISC),LD,IV,NBEG(ISC),LDD)
      CALL PLYBND(NBEG(ISC),LDD,B)
      VB = MAX(ABS(B(1)),ABS(B(2)))
      CALL FOXDAL(ISC,1)
*
      RETURN
      END
*
      SUBROUTINE CURVDT0(ND,NP,N)
*     ***************************
*
*     THIS SUBROUTINE PREPROCESSES PNT(J,I) FOR CURVDT. I=0,..,NP, J=1,..,ND.
*     J=1,2,3 ARE FOR X,Y,Z, AND J=4,5,6,7 ARE FOR R,G,B,A. ND IS 6 OR 7.
*     AFTER 2D PROJECTION, IF THERE ARE ANY CONSECUTIVE 2D SPACIALLY IDENTICAL
*     POINTS, THE REDUNDANT POINTS ARE DISCARDED, KEEPING ONLY THE LAST ONE
*     WITH ITS COLOR PROPERTY, RESULTING I=0,..,N.
*     PNT(1,I),PNT(2,I) : 2D PROJECTED X,Y COORDINATES IN [0.1]^2.
*     PNT(4,I)-PNT(7,I) : R,G,B,A.
*
*     CAUTION: THE ORIGINAL DATA PNT(J,I) IS OVERWRITTEN.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      CALL PROJPNT(NP)
*
      N = 0
      X0 = PNT(1,0)
      Y0 = PNT(2,0)
*
      DO 20 I=1,NP
      X = PNT(1,I)
      Y = PNT(2,I)
      IF(((X-X0)*(X-X0)+(Y-Y0)*(Y-Y0)).LT.1.D-20) GOTO 20
      N = N+1
      DO 10 J=1,ND
      PNT(J,N) = PNT(J,I)
 10   CONTINUE
      X0 = X
      Y0 = Y
 20   CONTINUE
*
      IF(N.EQ.0) THEN
         DO 30 J=1,ND
         PNT(J,0) = PNT(J,NP)
 30      CONTINUE
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE CURVDT(N,NCV,JX,KAP)
*     *******************************
*
*     THIS SUBROUTINE PREPARES THE GRPOLY GRADIENT COLOR CURVE DRAWING DATA.
*     CALL CURVDT0 AHEAD OF TIME.
*     OUT OF N+1 CURVDT0 PROCESSED POINTS, ALONG THE CURVE PROGRESSION,
*     NCV+1 PAIRS OF POINTS, THE LEFT AND RIGHT SIDES, ARE PREPARED IN
*     THE COSY MEMORY, FORMING A BAND WITH THE USER SPECIFIED LINE WIDTH.
*
*     THE LEFT SIDE POINTS:  K=0,..,NCV
*        (X1,Y1)(K) = ( CC(JX(1,1)+K) , CC(JX(2,1)+K) )
*     THE RIGHT SIDE POINTS: K=0,..,NCV
*        (X2,Y2)(K) = ( CC(JX(1,2)+K) , CC(JX(2,2)+K) )
*     THE COLOR VALUES OF THE PAIR (X1,Y1)(K) AND (X2,Y2)(K):
*        PNT(J,M), J=4,5,6,7 FOR R,G,B,A WHERE M = NC(JX(1,1)+K)
*
*     KAP: THE NUMBER OF AZIMUTHAL DIVISION OF A QUADRANT FOR A CAP FORMATION
*          KAP=4: EVERY 25 DEGREES,  KAP=3: EVERY 30 DEGREES
*
*     CAUTION: THE ORIGINAL DATA PNT(J,I) AS WELL AS PNTO(J,I) ARE OVERWRITTEN.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      PARAMETER(PI=3.141592653589793238462643383279502884197169399375D0)
      DOUBLE PRECISION SCALE(2),S(2)
      INTEGER JX(2,2)
*
*     THE SCALING RATIO IS 1.5 IN MOST OF THE COSY 2D GRAPHICS DRIVERS.
      SCALE(1) = 1.5D0
      SCALE(2) = 1.0D0
*
*     COMPUTE ANGLES  -- ANGLES IN THE SCALED SPACE
*     **************
*     PNT(3,I) : THE ANGLE OF THE 2D LINE SEGMENT FROM THE (I-1)-TH POINT TO
*                THE I-TH POINT. THE VALUE IS FROM -PI TO PI; COUNTERCLOCKWISE;
*                THE POSITIVE X DIRECTION IS 0.
*     PNTO(3,I): TPI=1/TAN((PI+PHI)/2) AT THE JOINT OF THE TWO LINE SEGMENTS
*                FORMED BY THE (I-1),I,(I+1)-TH POINTS,
*                WHERE PHI IS THE ANGLE FORMED BY THE TWO LINE SEGMENTS.
*                SPECIAL CASES: FOR PHI=0, TPI=0. FOR PHI=PI, TPI=-1E6.
*
      PNT(3,0) = 0.D0
      DO 20 I=1,N
      DO 10 J=1,2
      S(J) = (PNT(J,I)-PNT(J,I-1))*SCALE(J)
 10   CONTINUE
      PNT(3,I) = ATAN2(S(2),S(1))
 20   CONTINUE
      IF(N.GE.1) PNT(3,0) = PNT(3,1)
*
      PNTO(3,0) = 0.D0
      PNTO(3,N) = 0.D0
*
      DO 30 I=1,N-1
      CX = COS(PNT(3,I))
      CY = SIN(PNT(3,I))
      CXN = COS(PNT(3,I+1))
      CYN = SIN(PNT(3,I+1))
      PROD = CX*CXN+CY*CYN
      TPI = 0.D0
      IF((1-ABS(PROD)).GT.1.D-15) THEN
         PHI = ACOS(PROD)
         IF((CX*CYN-CY*CXN).LT.0.D0) PHI = -PHI
         TPI = 1.D0/TAN(0.5D0*(PI+PHI))
      ELSEIF(PROD.LT.0.D0) THEN
*        THE DIRECTION IS OPPOSITE, SO PHI=PI AND 0.5*(PI+PHI)=PI.
*        GIVE A NON-ZERO, FAIRLY LARGE NEGATIVE VALUE TO TPI.
         TPI = -1.D6
      ENDIF
      PNTO(3,I) = TPI
 30   CONTINUE
*
*     PREPARE POINTS FOR COLOR GRADIENT SHADING
*     *****************************************
*
      PM = PI/(2.D0*KAP)
      W = 0.00075D0*ACW
*
      NCV = -1
*
*   * THE STARTING POINT WITH THE HALF END CAP *
      I = 0
      ANG = PNT(3,I)
      DO 40 K=2*KAP,KAP,-1
      CALL CURVPS(NCV,I,ANG,0.D0,K,PM,W,SCALE,JX)
 40   CONTINUE
*
*   * THE POINTS BETWEEN *
      DO 60 I=1,N-1
      ANG = PNT(3,I)
      TPI = PNTO(3,I)
      IF(ABS(TPI).LE.0.2D0) THEN
*        THE BENDING ANGLE IS SHALLOW (~< PI/8 ) -- A MITER CONNECTION
*        [ .LE.1.0D0 SETS THE THRESHOLD PI/2.
*        [ .LE.0.5D0 SETS THE THRESHOLD AROUND 0.3*PI.
         CALL CURVPS(NCV,I,ANG,TPI,KAP,PM,W,SCALE,JX)
      ELSE
*        THE BENDING ANGLE IS SHARP   (~> PI/8 ) -- INSERT A HALF CAP
         DO 50 K=KAP,0,-1
         CALL CURVPS(NCV,I,ANG,0.D0,K,PM,W,SCALE,JX)
 50      CONTINUE
         ANG = PNT(3,I+1)
         CALL CURVPS(NCV,I,ANG,0.D0,KAP,PM,W,SCALE,JX)
      ENDIF
 60   CONTINUE
*
*   * THE ENDING POINT WITH THE HALF END CAP *
      I = N
      ANG = PNT(3,I)
      DO 70 K=KAP,0,-1
      CALL CURVPS(NCV,I,ANG,0.D0,K,PM,W,SCALE,JX)
 70   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE CURVPS(NCV,I,ANG,TPI,K,PM,W,SCALE,JX)
*     ************************************************
*
*     THIS SUBROUTINE COMPUTES AN ADDITIONAL, MOMENTARILY THE NCV-TH, PAIR OF
*     POINTS FOR CURVDT, STORING IN THE RESPECTIVE COSY MEMORY SLOTS.
*     NCV:      THE MOMENTARY NUMBER OF POINT PAIRS FOR CURDVT.
*               NCV IS INCREMENTALLY UPDATED.
*     I:        THE POINTER TO THE CURVDT0 PROCESSED DATA PNT(J,I).
*     ANG:      THE ANGLE OF THE LINE SEGMENT.
*     TPI:      THE TPI VALUE AT THE JOINT OF TWO LINE SEGMENTS.
*     K:        THE IDENTIFIER OF THE LOCATION IN A CAP.
*     PM:       THE AZIMUTHAL SPACING ANGLE TO FORM A CAP.
*     W:        THE HALF LINE WIDTH.
*     SCALE(J): THE 2D GRAPHICS SCALING FACTORS. J=1,2 FOR X,Y.
*     JX(J,LR): THE STARTING ADDRESSES OF THE COSY MEMORY SLOTS.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION C(2),S(2),SCALE(2)
      INTEGER JX(2,2)
*
      NCV = NCV+1
      NC(JX(1,1)+NCV) = I
      DO 10 LR=1,2
      ANGP = ANG+PM*K*(3-2*LR)
      TPIP = TPI*(3-2*LR)
      C(1) = COS(ANGP)+TPIP*COS(ANG)
      C(2) = SIN(ANGP)+TPIP*SIN(ANG)
      DO 10 J=1,2
      S(J) = C(J)*W/SCALE(J)
      CC(JX(J,LR)+NCV) = PNT(J,I)+S(J)
 10   CONTINUE
*
      RETURN
      END
*
      LOGICAL FUNCTION LGRDOT(X1,Y1,Z1)
*     *********************************
*
*     THIS FUNCTION DETERMINES IF THE GIVEN POINT IS WITHIN
*     THE COSY VIEWING AREA [0,1]^2 AFTER PROJECTION.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*
      LGRDOT = .TRUE.
      IF(IZOOM.EQ.0) RETURN
*
      CALL PROJEC(X1,Y1,Z1,XP1,YP1)
      IF(XP1.LT.0.D0.OR.XP1.GT.1.D0.OR.YP1.LT.0.D0.OR.YP1.GT.1.D0)
     *   LGRDOT = .FALSE.
*
      RETURN
      END
*
      LOGICAL FUNCTION LGRDRA(X1,Y1,Z1,X2,Y2,Z2)
*     ******************************************
*
*     THIS FUNCTION DETERMINES IF THE LINE BETWEEN GIVEN POINTS
*     MIGHT INTERSECT THE COSY VIEWING AREA [0,1]^2 AFTER PROJECTION.
*
*     NOTE: WHEN USING A POLYGON CONCEPT IN THE RESPECTIVE "...DRA" ROUTINE,
*     REMOVING THE LINE SEGMENT MAY MESS UP THE POLYGON.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*
      LGRDRA = .TRUE.
      IF(IZOOM.EQ.0) RETURN
*
      CALL PROJEC(X1,Y1,Z1,XP1,YP1)
      CALL PROJEC(X2,Y2,Z2,XP2,YP2)
      IF(MAX(XP1,XP2).LT.0.D0.OR.MIN(XP1,XP2).GT.1.D0.OR.
     *   MAX(YP1,YP2).LT.0.D0.OR.MIN(YP1,YP2).GT.1.D0)
     *   LGRDRA = .FALSE.
*
      RETURN
      END
*
      LOGICAL FUNCTION LGRTRI(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)
*     ***************************************************
*
*     THIS FUNCTION DETERMINES IF THE TRIANGLE SPANNED BY GIVEN POINTS
*     MIGHT INTERSECT THE COSY VIEWING AREA [0,1]^2 AFTER PROJECTION.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*
      LGRTRI = .TRUE.
      IF(IZOOM.EQ.0) RETURN
*
      CALL PROJEC(X1,Y1,Z1,XP1,YP1)
      CALL PROJEC(X2,Y2,Z2,XP2,YP2)
      CALL PROJEC(X3,Y3,Z3,XP3,YP3)
      IF(MAX(XP1,XP2,XP3).LT.0.D0.OR.MIN(XP1,XP2,XP3).GT.1.D0.OR.
     *   MAX(YP1,YP2,YP3).LT.0.D0.OR.MIN(YP1,YP2,YP3).GT.1.D0)
     *   LGRTRI = .FALSE.
*
      RETURN
      END
*
      LOGICAL FUNCTION LGRPLY(IAG)
*     ****************************
*
*     THIS FUNCTION DETERMINES IF THE GRPOLY POLYNOMIAL GRAPHICS OBJECT
*     WITH THE STARTING POSITION IAG IN COSY MEMORY MIGHT INTERSECT
*     THE COSY VIEWING AREA [0,1]^2 AFTER PROJECTION.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*
      INTEGER IPLEN(7),IPRER(7)
      DOUBLE PRECISION B3(3,2),B2(2,2),B(2),X(2)
      LOGICAL LCR
      PARAMETER(IBOUND=1)
*
      LGRPLY = .TRUE.
      IF(IZOOM.EQ.0) RETURN
*
*     COMPUTE 3D BOUNDS
*
      IF(IBOUND.EQ.1) THEN
*
*        ACCURATE BOUNDING
*
         CALL PLYPRE(IAG,IOBJ,NP1,NP2,LCR,0)
         CALL PLYBD(IAG,IOBJ,NP1,NP2,B3)
      ELSE
*
*        VERY ROUGH BOUNDING: QUICKLY ESTIMATED, BUT IT MAY NOT BE ABLE TO
*        EXCLUDE UNNECESSARY POLYNOMIAL GRAPHICS OBJECTS.
*
         CALL PLYLEN(IAG,IPLEN,IPRER)
         IPOS = IAG+2
         DO 20 J=1,3
         CALL PLYBND(IPOS,IPLEN(J)-IPRER(J),B)
         DO 10 L=1,2
         B3(J,L) = B(L)
 10      CONTINUE
         IPOS = IPOS+IPLEN(J)
 20      CONTINUE
      ENDIF
*
*     PROJECT 8 CORNER POINTS OF THE 3D BOUNDS AND COMPUTE 2D BOUNDS
*
      IBEG = 0
      DO 50 K=1,2
      DO 50 L=1,2
      DO 50 M=1,2
      CALL PROJEC(B3(1,K),B3(2,L),B3(3,M),X(1),X(2))
      IF(IBEG.EQ.0) THEN
         IBEG = 1
         DO 30 J=1,2
         DO 30 LM=1,2
         B2(J,LM) = X(J)
 30      CONTINUE
      ELSE
         DO 40 J=1,2
         B2(J,1) = MIN(B2(J,1),X(J))
         B2(J,2) = MAX(B2(J,2),X(J))
 40      CONTINUE
      ENDIF
 50   CONTINUE
*
*     TEST 2D BOUNDS
*
      DO 60 J=1,2
      IF(B2(J,1).GT.1.D0.OR.B2(J,2).LT.0.D0) LGRPLY = .FALSE.
 60   CONTINUE
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     GENERIC DRIVER ROUTINES                                                 *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE GENBEG
*     *****************
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      DO 10 J=1,3
      XC(J) = 0.D0
      XP(J) = 0.D0
      ACC(J) = 0.D0
 10   CONTINUE
      ACC(4) = 1.D0
      ACW = 1.D0
      EPSXYZ = 0.D0
      EPSCOL = 0.D0
      LCM = 0
      LWR = 0
      LDOTW = 1
      LDOTC = 0
      IPLY = 0
*
      RETURN
      END
*
      SUBROUTINE GENMOV(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      XP(1) = XC(1)
      XP(2) = XC(2)
      XP(3) = XC(3)
      XC(1) = XXX
      XC(2) = YYY
      XC(3) = ZZZ
      LCM = 0
*
      RETURN
      END
*
      SUBROUTINE GENDRA(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      XP(1) = XC(1)
      XP(2) = XC(2)
      XP(3) = XC(3)
      XC(1) = XXX
      XC(2) = YYY
      XC(3) = ZZZ
      LCM = 1
*
      RETURN
      END
*
      SUBROUTINE GENDOT(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      XP(1) = XC(1)
      XP(2) = XC(2)
      XP(3) = XC(3)
      XC(1) = XXX
      XC(2) = YYY
      XC(3) = ZZZ
      LCM = 2
*
      RETURN
      END
*
      SUBROUTINE GENTRI(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      XP(1) = XC(1)
      XP(2) = XC(2)
      XP(3) = XC(3)
      XC(1) = XXX
      XC(2) = YYY
      XC(3) = ZZZ
      LCM = 3
*
      RETURN
      END
*
      SUBROUTINE GENPLY(IA)
*     *********************
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      IPLY = IPLY+1
      CALL PLYPOS(IA,XC)
      LCM = 0
*
      RETURN
      END
*
      SUBROUTINE GENCHA
*     *****************
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      LCM = 0
      LDOTC = 0
*
      RETURN
      END
*
      SUBROUTINE GENCOL(CLR)
*     **********************
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CLR(4)
*
      LCM = 0
      DO 10 I=1,4
      ACC(I) = CLR(I)
 10   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE GENWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      LCM = 0
      ACW = AJ
      LDOTW = 1
*
      RETURN
      END
*
      SUBROUTINE GENSTL(ID,NLEN,I)
*     ****************************
*
*     THIS SUBROUTINE REGISTERS USER INPUT OPTIONAL INFORMATION
*     DURING THE GRPRI OUTPUT.
*
*     IT ALSO CAN BE USED TO HANDLE SOME OTHER SIMILAR USER SERVICE PROCEDURES.
*     ADJUST THE FLOWS CONTROLLED BY COMPUTED GOTO STATEMENTS AND SIMILAR LINES
*     IN GRPRI AND GEPRE.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      IF(ID.EQ.9) THEN
         LWR = NINT(CC(I))
      ELSEIF(ID.EQ.10) THEN
         EPSXYZ = CC(I)
         EPSCOL = CC(I+1)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GENEND
*     *****************
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     TTY (LOW RESOLUTION ASCII) DRIVER ROUTINES                              *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE TTYBEG
*     *****************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CHARACTER A*1
      COMMON /TTYCHX/ A(79,23)
*
      DO 10 I=1,79
      DO 10 J=1,23
      A(I,J) = ' '
 10   CONTINUE
*
      CALL GENBEG
*
      RETURN
      END
*
      SUBROUTINE TTYMOV(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL GENMOV(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE TTYDRA(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER A*1
      COMMON /TTYCHX/ A(79,23)
*
      CALL PROJEC(XC(1),XC(2),XC(3),X,Y)
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      NS = 1+NINT(MAX(ABS(X-XX)*79.D0,ABS(Y-YY)*23.D0))
      DS = 1.D0/NS
*
      DO 10 I=0,NS
      S = I*DS
      PX = X+S*(XX-X)
      PY = Y+S*(YY-Y)
      IX = 1+INT(PX*78.999D0)
      IY = 1+INT(PY*22.999D0)
      IF(IX.GE.1.AND.IX.LE.79.AND.IY.GE.1.AND.IY.LE.23) A(IX,IY) = '*'
 10   CONTINUE
*
      CALL GENDRA(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE TTYDOT(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER A*1
      COMMON /TTYCHX/ A(79,23)
*
      CALL PROJEC(XXX,YYY,ZZZ,X,Y)
      IX = 1+INT(X*78.999D0)
      IY = 1+INT(Y*22.999D0)
      IF(IX.GE.1.AND.IX.LE.79.AND.IY.GE.1.AND.IY.LE.23) A(IX,IY) = '.'
*
      CALL GENDOT(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE TTYTRI(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER A*1
      COMMON /TTYCHX/ A(79,23)
*
      PX = XP(1)
      PY = XP(2)
      PZ = XP(3)
      CX = XC(1)
      CY = XC(2)
      CZ = XC(3)
      CALL TTYDRA(XXX,YYY,ZZZ)
      CALL TTYDRA(PX,PY,PZ)
      IF(LCM.NE.3) THEN
         CALL TTYDRA(CX,CY,CZ)
      ELSE
         CALL TTYMOV(CX,CY,CZ)
      ENDIF
*
      CALL GENTRI(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE TTYPLY(IA,IPST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      LOGICAL LGRPLY,LCR
*
      IF(IPST.EQ.-1) GOTO 90
      IF(.NOT.LGRPLY(IA)) GOTO 90
*
      CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
      ND = 3
*
*     CURVE
*     *****
*
      IF(IOBJ.EQ.1) THEN
         CALL PLYEV(IA,ND,0.D0,NP1)
         CALL TTYMOV(PNT(1,0),PNT(2,0),PNT(3,0))
         DO 10 I=1,NP1
         CALL TTYDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 10      CONTINUE
*
*     SURFACE
*     *******
*
      ELSEIF(IOBJ.EQ.2) THEN
         CALL PLYEV(IA,ND,-1.D0,NP1)
         CALL TTYMOV(PNT(1,0),PNT(2,0),PNT(3,0))
         DO 20 I=1,NP1
         CALL TTYDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 20      CONTINUE
*
         DO 60 I2=1,NP2
         DO 30 I=0,NP1
         DO 30 J=1,ND
 30      PNTO(J,I) = PNT(J,I)
         CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
         DO 40 I=0,NP1
         CALL TTYMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
         CALL TTYDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 40      CONTINUE
*
         CALL TTYMOV(PNT(1,0),PNT(2,0),PNT(3,0))
         DO 50 I=1,NP1
         CALL TTYDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 50      CONTINUE
 60      CONTINUE
      ENDIF
*
 90   CALL GENPLY(IA)
*
      RETURN
      END
*
      SUBROUTINE TTYCHA(STR,LST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER A*1
      COMMON /TTYCHX/ A(79,23)
      CHARACTER STR*1024
*
      CALL PROJEC(XC(1),XC(2),XC(3),X,Y)
      X1 = X
      IY = 1+INT(Y*22.999D0)
      IF(IY.GE.1.AND.IY.LE.23) THEN
         DO 10 I=1,LST
         IX = 1+INT(X1*78.999D0)
         IF(IX.GE.1.AND.IX.LE.79) A(IX,IY) = STR(I:I)
         X1 = X1 + 1.D0/79.D0
 10      CONTINUE
      ENDIF
*
      CALL GENCHA
*
      RETURN
      END
*
      SUBROUTINE TTYCOL(CLR)
*     **********************
*
      DOUBLE PRECISION CLR(4)
*
      CALL GENCOL(CLR)
*
      RETURN
      END
*
      SUBROUTINE TTYWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL GENWID(AJ)
*
      RETURN
      END
*
      SUBROUTINE TTYEND(IUNIT,NARI)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CHARACTER A*1
      COMMON /TTYCHX/ A(79,23)
      CHARACTER STR*100
*
      IF(NARI.GT.0) THEN
         WRITE(STR,'(A,I0)') '$$$ WARNING, GRPOLY ARITHMETIC FAILURES: '
     *             ,NARI
         LST = ILAST(STR,1,100)
         X = 0.007D0
         Y = 0.010D0
         X1 = X
         IY = 1+INT(Y*22.999D0)
         IF(IY.GE.1.AND.IY.LE.23) THEN
            DO 10 I=1,LST
            IX = 1+INT(X1*78.999D0)
            IF(IX.GE.1.AND.IX.LE.79) A(IX,IY) = STR(I:I)
            X1 = X1 + 1.D0/79.D0
 10         CONTINUE
         ENDIF
      ENDIF
*
      WRITE(IUNIT,'(1X,79A)') ((A(I,J),I=1,79),J=23,1,-1)
*
      CALL GENEND
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     METAFILE DRIVER ROUTINES                                                *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE ASCBEG
*     *****************
*
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*-----GRAPH FILE OUTPUT -----------------------------------------------------
      INTEGER ICNT
      CHARACTER FBASE*20
      COMMON /GRFILE/ ICNT,FBASE
*----------------------------------------------------------------------------
*
      CHARACTER NUMBER*3
*
      WRITE(NUMBER,'(I3.3)') MOD(ICNT,1000)
      IL = ILAST(FBASE,1,20)
      OPEN(4,FILE=FBASE(1:IL)//NUMBER//'.dat',STATUS='UNKNOWN')
      ICNT = ICNT+1
*
      WRITE(4,'(A)') 'BEG   GRAPHICS METAFILE CREATED BY COSY INFINITY'
      WRITE(4,'(A,2E24.16)') 'PRJ ',PHI,THETA
      IF(IZOOM.EQ.1) WRITE(4,'(A,6E24.16)') 'ZOO ',BD3(1,1),BD3(1,2),
     *    BD3(2,1),BD3(2,2),BD3(3,1),BD3(3,2)
*
      CALL GENBEG
*
      RETURN
      END
*
      SUBROUTINE ASCMOV(X,Y,Z)
*     ************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      WRITE(4,'(A,3E24.16)') 'MOV ',X,Y,Z
*
      CALL GENMOV(X,Y,Z)
*
      RETURN
      END
*
      SUBROUTINE ASCDRA(X,Y,Z)
*     ************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      WRITE(4,'(A,3E24.16)') 'DRA ',X,Y,Z
*
      CALL GENDRA(X,Y,Z)
*
      RETURN
      END
*
      SUBROUTINE ASCDOT(X,Y,Z)
*     ************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      WRITE(4,'(A,3E24.16)') 'DOT ',X,Y,Z
*
      CALL GENDOT(X,Y,Z)
*
      RETURN
      END
*
      SUBROUTINE ASCTRI(X,Y,Z)
*     ************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      WRITE(4,'(A,3E24.16)') 'TRI ',X,Y,Z
*
      CALL GENTRI(X,Y,Z)
*
      RETURN
      END
*
      SUBROUTINE ASCPLY(IA,IPST)
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
      INTEGER IPLEN(7),IPRER(7)
*
      IF(IPST.NE.-1) THEN
         WRITE(4,'(A)') 'PLY'
         CALL PLYLEN(IA,IPLEN,IPRER)
*
         IPOS = IA+2
         DO 20 I=1,7
         DO 10 J=IPOS,IPOS+IPLEN(I)-IPRER(I)-1
         IP2 = INT(NC(J)/100)
         IP1 = NC(J)-IP2*100
         WRITE(4,'(A,I1,2(X,I2),X,E24.16)') 'PL',I,IP1,IP2,CC(J)
 10      CONTINUE
         IF(IPRER(I).EQ.1) THEN
            WRITE(4,'(A,I1,5X,A,X,E24.16)') 'PL',I,'R',CC(J)
         ENDIF
         IPOS = IPOS+IPLEN(I)
 20      CONTINUE
         WRITE(4,'(A)') 'PLE'
      ELSE
         WRITE(4,'(A)') 'PLY   -ARITHMETIC FAILURE-'
      ENDIF
*
      CALL GENPLY(IA)
*
      RETURN
      END
*
      SUBROUTINE ASCCHA(STR,LST)
*     **************************
*
      CHARACTER STR*1024
*
      WRITE(4,'(A)') 'CHA   '//STR(1:LST)
*
      CALL GENCHA
*
      RETURN
      END
*
      SUBROUTINE ASCCOL(CLR)
*     **********************
*
      DOUBLE PRECISION CLR(4)
*
      WRITE(4,'(A,4F12.8)') 'RGB ',(CLR(J),J=1,4)
*
      CALL GENCOL(CLR)
*
      RETURN
      END
*
      SUBROUTINE ASCWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      WRITE(4,'(A,F9.4)') 'WID',AJ
*
      CALL GENWID(AJ)
*
      RETURN
      END
*
      SUBROUTINE ASCSTL(ID,NLEN,I)
*     ****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CALL GENSTL(ID,NLEN,I)
*
      IF(ID.EQ.9) THEN
         WRITE(4,'(A,I2)') 'LWR ',LWR
      ELSEIF(ID.EQ.10) THEN
         WRITE(4,'(A,2E24.16)') 'EPS ',EPSXYZ,EPSCOL
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE ASCEND
*     *****************
*
      WRITE(4,'(A)') 'END   COSY PICTURE'
      CLOSE(4)
*
      CALL GENEND
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     PGPLOT DRIVER ROUTINES ORIGINALLY WITTEN BY JENS HOEFKENS               *
*     UPDATED BY K. MAKINO                                                    *
*                                                                             *
*     USE OF PGPLOT KINDLY ACKNOWLEDGED, SOURCES ARE AVAILABLE AT             *
*     http://www.astro.caltech.edu/~tjp/pgplot/ (SUBJECT TO COPYRIGHT)        *
*                                                                             *
*******************************************************************************
*
      BLOCKDATA PGPDAT
*     ****************
*
*     THIS BLOCKDATA INITIALIZES PGPLOT RELATED COMMON VARIABLES
*
      INTEGER PGPIDS(2,-112:-101)
      COMMON /PGP/ PGPIDS
      DATA PGPIDS / 24*0 /
*
      END
*
      SUBROUTINE PGPBEG(IUNIT)
*     ************************
*
*     THE PGPLOT INTERFACE CURRENTLY SUPPORTS THE UNITS -101 -> -110
*     -20 (IUNIT1=-111) OUTPUTS TO POSTSCRIPT
*     -22 (IUNIT1=-112) OUTPUTS TO LATEX
*
*-----GRAPH FILE OUTPUT -----------------------------------------------------
      INTEGER ICNT
      CHARACTER FBASE*20
      COMMON /GRFILE/ ICNT,FBASE
*----------------------------------------------------------------------------
*
      CHARACTER NUMBER*3,PGPDT*8,PGPDSC*80
      LOGICAL HAVEX,HAVEPX,HAVEPS,HVELTX,HAVEWV,FIRST
      INTEGER PGOPEN,MAXDEV,TLEN,DLEN,INTERA
      INTEGER PGPIDS(2,-112:-101)
      COMMON /PGP/ PGPIDS
*
      SAVE HAVEX,HAVEPX,HAVEPS,HVELTX,HAVEWV,FIRST
*
      DATA FIRST  / .TRUE. /
*
      IF(FIRST) THEN
         HAVEX  = .FALSE.
         HAVEPX = .FALSE.
         HAVEPS = .FALSE.
         HVELTX = .FALSE.
         HAVEWV = .FALSE.
*PGP     CALL PGQNDT(MAXDEV)                                             *PGP
*PGP     DO 100 I=1,MAXDEV                                               *PGP
*PGP        CALL PGQDT(I,PGPDT,TLEN,PGPDSC,DLEN,INTERA)                  *PGP
*PGP        IF(PGPDT.EQ.'/XWINDOW') THEN                                 *PGP
*PGP           HAVEX = .TRUE.                                            *PGP
*PGP        ELSEIF(PGPDT.EQ.'/XSERVE') THEN                              *PGP
*PGP           HAVEPX = .TRUE.                                           *PGP
*PGP        ELSEIF(PGPDT.EQ.'/WV') THEN                                  *PGP
*PGP           HAVEWV = .TRUE.                                           *PGP
*PGP        ELSEIF(PGPDT.EQ.'/CPS') THEN                                 *PGP
*PGP           HAVEPS = .TRUE.                                           *PGP
*PGP        ELSEIF(PGPDT.EQ.'/LATEX') THEN                               *PGP
*PGP           HVELTX = .TRUE.                                           *PGP
*PGP        ENDIF                                                        *PGP
 100     CONTINUE
         FIRST = .FALSE.
      ENDIF
*
      IUNIT1 = IUNIT
      IF(IUNIT.EQ.-20) IUNIT1 = -111
      IF(IUNIT.EQ.-22) IUNIT1 = -112
*
      IF(PGPIDS(2,IUNIT1).EQ.0) THEN
         IF((IUNIT1.LE.-101).AND.(IUNIT1.GE.-110)) THEN
            PGPIDS(2,IUNIT1) = 1
         ELSEIF(IUNIT1.EQ.-112) THEN
            PGPIDS(2,IUNIT1) = 3
         ELSEIF(IUNIT1.EQ.-111) THEN
            PGPIDS(2,IUNIT1) = 4
         ELSE
            PRINT*,'@@@ ERROR IN PGPBEG, PGPLOT INITIALIZATION'
         ENDIF
      ENDIF
*
      IF(PGPIDS(1,IUNIT1).EQ.0) THEN
         IF(PGPIDS(2,IUNIT1).EQ.1) THEN
            IF(HAVEWV) THEN
*PGP           PGPIDS(1,IUNIT1) = PGOPEN('/WV')                          *PGP
            ELSEIF(HAVEPX) THEN
*PGP           PGPIDS(1,IUNIT1) = PGOPEN('/XSERVE')                      *PGP
            ELSEIF(HAVEX) THEN
*PGP           PGPIDS(1,IUNIT1) = PGOPEN('/XWINDOW')                     *PGP
            ELSE
*PGP           PGPIDS(1,IUNIT1) = PGOPEN('/NULL')                        *PGP
               PRINT*,'$$$ WARNING IN PGPBEG, NO PGPLOT GRAPHICS '//
     *                'DEVICE FOUND. AUTOMATICALLY SELECTED /NULL'
            ENDIF
         ELSEIF(PGPIDS(2,IUNIT1).EQ.3.OR.PGPIDS(2,IUNIT1).EQ.4) THEN
            WRITE(NUMBER,'(I3.3)') MOD(ICNT,1000)
            IL = ILAST(FBASE,1,20)
            IF(PGPIDS(2,IUNIT1).EQ.3) THEN
*PGP           PGPIDS(1,IUNIT1) =                                        *PGP
*PGP *            PGOPEN(FBASE(1:IL)//NUMBER//'.tex/LATEX')              *PGP
            ELSEIF(PGPIDS(2,IUNIT1).EQ.4) THEN
*PGP           PGPIDS(1,IUNIT1) = PGOPEN(FBASE(1:IL)//NUMBER//'.ps/CPS') *PGP
            ENDIF
            ICNT = ICNT+1
         ENDIF
*PGP     CALL PGSCR( 0,1.0,1.0,1.0)                                      *PGP
*PGP     CALL PGSCR( 1,0.0,0.0,0.0)                                      *PGP
*PGP     CALL PGSCR(10,1.0,1.0,1.0)                                      *PGP
         IF(PGPIDS(1,IUNIT1).LE.0) THEN
            PRINT*,'@@@ ERROR IN PGPBEG, PGPLOT INITIALIZATION'
            RETURN
         ENDIF
*PGP     CALL PGASK(.FALSE.)                                             *PGP
      ELSE
*PGP     CALL PGSLCT(PGPIDS(1,IUNIT1))                                   *PGP
      ENDIF
*
*PGP  CALL PGPAGE                                                        *PGP
*PGP  CALL PGSVP (0.0, 1.0, 0.0, 1.0)                                    *PGP
*PGP  CALL PGWNAD(0.0, 1.5, 0.0, 1.0)                                    *PGP
*PGP  CALL PGSWIN(0.0, 1.0, 0.0, 1.0)                                    *PGP
*PGP  CALL PGSCF(1)                                                      *PGP
*PGP  CALL PGSCI(1)                                                      *PGP
*PGP  CALL PGSLW(4)                                                      *PGP
*PGP  CALL PGSCH(0.92)                                                   *PGP
*     PGSFS: THE FILL STYLE FOR PGPOLY.  1: SOLID (DEFAULT),  2: OUTLINE
*PGP  CALL PGSFS(1)                                                      *PGP
*
      CALL GENBEG
*
      RETURN
      END
*
      SUBROUTINE PGPMOV(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      REAL X,Y
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      X=REAL(XX)
      Y=REAL(YY)
*PGP  CALL PGMOVE(X,Y)                                                   *PGP
*
      CALL GENMOV(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE PGPDRA(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      REAL X,Y
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      X=REAL(XX)
      Y=REAL(YY)
*PGP  CALL PGDRAW(X,Y)                                                   *PGP
*
      CALL GENDRA(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE PGPDOT(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      REAL X,Y
      LOGICAL LGRDOT
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      X=REAL(XX)
      Y=REAL(YY)
*PGP  CALL PGMOVE(X,Y)                                                   *PGP
*PGP  IF(LGRDOT(XXX,YYY,ZZZ)) CALL PGPT1(X,Y,-1)                         *PGP
*
      CALL GENDOT(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE PGPTRI(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      REAL TX(3),TY(3)
      LOGICAL LGRTRI
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      TX(1) = REAL(XX)
      TY(1) = REAL(YY)
      IF(LGRTRI(XXX,YYY,ZZZ,XP(1),XP(2),XP(3),XC(1),XC(2),XC(3))) THEN
         CALL PROJEC(XP(1),XP(2),XP(3),PXX,PYY)
         TX(2) = REAL(PXX)
         TY(2) = REAL(PYY)
         CALL PROJEC(XC(1),XC(2),XC(3),CXX,CYY)
         TX(3) = REAL(CXX)
         TY(3) = REAL(CYY)
         IF(LWR.EQ.1) THEN
*PGP        CALL PGSFS(2)                                                *PGP
*PGP        CALL PGPOLY(3,TX,TY)                                         *PGP
*PGP        CALL PGSFS(1)                                                *PGP
         ELSE
*PGP        CALL PGPOLY(3,TX,TY)                                         *PGP
         ENDIF
      ENDIF
*PGP  CALL PGMOVE(TX(1),TY(1))                                           *PGP
*
      CALL GENTRI(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE PGPPLY(IA,IPST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION OACC(4),CLR(4),PC(7)
      LOGICAL LGRPLY,LCR
*
      IF(IPST.EQ.-1) GOTO 90
      IF(.NOT.LGRPLY(IA)) GOTO 90
*
      CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
      IF(LCR) THEN
         ND = 6
         DO 10 K=1,4
         OACC(K) = ACC(K)
 10      CONTINUE
      ELSE
         ND = 3
      ENDIF
*
*     CURVE
*     *****
*
      IF(IOBJ.EQ.1) THEN
         CALL PLYEV(IA,ND,0.D0,NP1)
         CALL PGPMOV(PNT(1,0),PNT(2,0),PNT(3,0))
         DO 20 I=1,NP1
         DO 11 J=4,ND
 11      CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
         IF(LCR) CALL PGPCOL(CLR)
         CALL PGPDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 20      CONTINUE
*
*     SURFACE
*     *******
*
      ELSEIF(IOBJ.EQ.2) THEN
*
*      * VARYING OR SOLID COLOR TRIANGLE FILL PAINTING *
         IF(LWR.EQ.0) THEN
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            DO 40 I2=1,NP2
            DO 21 I=0,NP1
            DO 21 J=1,ND
 21         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 40 I=1,NP1
            DO 22 J=1,ND
 22         PC(J) = (PNTO(J,I-1)+PNTO(J,I)+PNT(J,I-1)+PNT(J,I))*0.25D0
*
            CALL PGPMOV(PNT(1,I),PNT(2,I),PNT(3,I))
            CALL PGPMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            DO 31 J=4,ND
 31         CLR(J-3) = (PNT(J,I)+PNTO(J,I)+PC(J))/3.D0
            IF(LCR) CALL PGPCOL(CLR)
            CALL PGPTRI(PC(1),PC(2),PC(3))
*
            DO 32 J=4,ND
 32         CLR(J-3) = (PNTO(J,I)+PC(J)+PNTO(J,I-1))/3.D0
            IF(LCR) CALL PGPCOL(CLR)
            CALL PGPTRI(PNTO(1,I-1),PNTO(2,I-1),PNTO(3,I-1))
*
            DO 33 J=4,ND
 33         CLR(J-3) = (PC(J)+PNTO(J,I-1)+PNT(J,I-1))/3.D0
            IF(LCR) CALL PGPCOL(CLR)
            CALL PGPTRI(PNT(1,I-1),PNT(2,I-1),PNT(3,I-1))
*
            CALL PGPMOV(PC(1),PC(2),PC(3))
            DO 34 J=4,ND
 34         CLR(J-3) = (PNT(J,I-1)+PC(J)+PNT(J,I))/3.D0
            IF(LCR) CALL PGPCOL(CLR)
            CALL PGPTRI(PNT(1,I),PNT(2,I),PNT(3,I))
 40         CONTINUE
*
*      * WIRE FRAME QUADRILATERALS *
         ELSE
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            CALL PGPMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 50 I=1,NP1
            DO 41 J=4,ND
 41         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL PGPCOL(CLR)
            CALL PGPDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 50         CONTINUE
*
            DO 80 I2=1,NP2
            DO 51 I=0,NP1
            DO 51 J=1,ND
 51         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 60 I=0,NP1
            DO 52 J=4,ND
 52         CLR(J-3) = (PNTO(J,I)+PNT(J,I))*0.5D0
            IF(LCR) CALL PGPCOL(CLR)
            CALL PGPMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            CALL PGPDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 60         CONTINUE
*
            CALL PGPMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 70 I=1,NP1
            DO 61 J=4,ND
 61         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL PGPCOL(CLR)
            CALL PGPDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 70         CONTINUE
 80         CONTINUE
         ENDIF
      ENDIF
*
      IF(LCR) CALL PGPCOL(OACC)
*
 90   CALL GENPLY(IA)
*
      RETURN
      END
*
      SUBROUTINE PGPCHA(STR,LST)
*     **************************
*
      CHARACTER STR*1024
      REAL X,Y
*
*PGP  CALL PGQPOS(X,Y)                                                   *PGP
*PGP  CALL PGTEXT(X,Y,STR(1:LST))                                        *PGP
*
*     LATEX FILE OUTPUT (IUNIT=-22) ADVANCES THE POSITION X,Y
*     DESPITE OF THE POSITION RESET BY THE NEXT PGMOVE CALL
*
*PGP  CALL PGMOVE(X,Y)                                                   *PGP
*
      CALL GENCHA
*
      RETURN
      END
*
      SUBROUTINE PGPCOL(CLR)
*     **********************
*
      DOUBLE PRECISION CLR(4)
*
*PGP  CALL PGSCR(16,REAL(CLR(1)),REAL(CLR(2)),REAL(CLR(3)))              *PGP
*PGP  CALL PGSCI(16)                                                     *PGP
*
      CALL GENCOL(CLR)
*
      RETURN
      END
*
      SUBROUTINE PGPWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      IW = MAX(NINT(2.D0*AJ+2.D0),1)
*PGP  CALL PGSLW(IW)                                                     *PGP
*
      CALL GENWID(AJ)
*
      RETURN
      END
*
      SUBROUTINE PGPEND(IUNIT,NARI)
*     *****************************
*
      INTEGER PGPIDS(2,-112:-101)
      COMMON /PGP/ PGPIDS
      CHARACTER STR*100
*
      IF(NARI.GT.0) THEN
*PGP     CALL PGSCI(1)                                                   *PGP
*PGP     CALL PGSLW(4)                                                   *PGP
         WRITE(STR,'(A,I0)') '$$$ WARNING, GRPOLY ARITHMETIC FAILURES: '
     *             ,NARI
         L = ILAST(STR,1,100)
*PGP     CALL PGTEXT(0.007,0.01,STR(1:L))                                *PGP
      ENDIF
*
      IUNIT1 = IUNIT
      IF(IUNIT.EQ.-20) IUNIT1 = -111
      IF(IUNIT.EQ.-22) IUNIT1 = -112
*
      IF(PGPIDS(2,IUNIT1).GT.2) THEN
*PGP     CALL PGCLOS                                                     *PGP
         PGPIDS(1,IUNIT1) = 0
      ENDIF
*
      CALL GENEND
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     GRWIN DRIVER ROUTINES ORIGINALLY WRITTEN BY ALEXEY POKLONSKIY           *
*     UPDATED BY K. MAKINO                                                    *
*                                                                             *
*     USE OF GRWIN 0.99.9b KINDLY ACKNOWLEDGED, SOURCES ARE AVAILABLE AT      *
*     http://spdg1.sci.shizuoka.ac.jp/grwinlib/english                        *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE GRWBEG(IUNIT)
*     ************************
*
      INTEGER IGW(-120:-111)
      LOGICAL FIRST
      SAVE FIRST,IGW
      DATA FIRST / .TRUE. /
      DATA IGW / 10*0 /
*
      IF(FIRST) THEN
*GRW     CALL GWINIT(IRTN)                                               *GRW
         FIRST = .FALSE.
         IF(IUNIT.EQ.-111) THEN
*
*        IF THE FIRST CALL TO GRWBEG IS WITH UNIT -111, BOTH THE #2 AND #1
*        WINDOWS ARE OPENED TO ACHIEVE PROPER TILING OF THESE TWO WINDOWS.
*
            NW = 2
            IRTN = 0
*GRW        CALL GWOPENX(IRTN,NW,3000,2000,0,19,-1,'COSY-'//             *GRW
*GRW *           CHAR(48+MOD(NW,10))//'   ')                             *GRW
            IF(IRTN.EQ.0) PRINT*,'$$$ WARNING IN GRWBEG, OPEN FAILURE'
*GRW        IGW(-112) = IRTN                                             *GRW
*GRW        CALL GWARRANGE(IRTN,2)                                       *GRW
*
            NW = 1
            IRTN = 0
*GRW        CALL GWOPENX(IRTN,NW,3000,2000,0,19,-1,'COSY-'//             *GRW
*GRW *           CHAR(48+MOD(NW,10))//'   ')                             *GRW
            IF(IRTN.EQ.0) PRINT*,'$$$ WARNING IN GRWBEG, OPEN FAILURE'
*GRW        IGW(-111) = IRTN                                             *GRW
*GRW        CALL GWARRANGE(IRTN,2)                                       *GRW
         ENDIF
      ENDIF
*
*     IF WINDOW NOT OPEN YET, OPEN IT CORRESPONDING TO COSY UNIT NUMBER
      IF(IGW(IUNIT).EQ.0) THEN
         NW = -IUNIT-110
         IRTN = 0
*GRW     CALL GWOPENX(IRTN,NW,3000,2000,0,19,-1,'COSY-'//                *GRW
*GRW *        CHAR(48+MOD(NW,10))//'   ')                                *GRW
         IF(IRTN.EQ.0) PRINT*,'$$$ WARNING IN GRWBEG, OPEN FAILURE'
*GRW     IGW(IUNIT) = IRTN                                               *GRW
      ELSE
         NW = IGW(IUNIT)
      ENDIF
*
*GRW  CALL GWSELECT(IRTN,NW)                                             *GRW
*GRW  CALL GWCLEAR(IRTN,-1)                                              *GRW
*
*     SET INITIAL PEN, BRUSH AND TEXT PROPERTIES
*GRW  CALL GWSETPEN(IRTN,0,1,5,4)                                        *GRW
*GRW  CALL GWSETBRS(IRTN,0,1,-1)                                         *GRW
*GRW  CALL GWSETTXT(IRTN, 0.0, 5.0, 1,0,-1,' ')                          *GRW
*
*     SETS WORLD COORDINATES, AND THUS ASPECT RATIO OF WINDOW, AND POSITION
*     - IN A GRZOOM VIEW, EXTRA CONTENTS AT RIGHT MIGHT BE SEEN WHEN
*       CHANGING THE WINDOW SIZE?
*GRW  CALL GWVPORT(IRTN, 0.0, 0.0, 1.5, 1.0)                             *GRW
*GRW  CALL GWINDOW(IRTN, 0.0, 0.0, 1.0, 1.0)                             *GRW
*GRW  CALL GWMOVE2(IRTN, 0.0, 0.0)                                       *GRW
*
*     ARRANGE WINDOWS (2 VERTICALLY)
*GRW  CALL GWARRANGE(IRTN,2)                                             *GRW
*
      CALL GENBEG
*
      RETURN
      END
*
      SUBROUTINE GRWMOV(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      REAL X,Y
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      X=REAL(XX)
      Y=REAL(YY)
*GRW  CALL GWMOVE2(IRTN,X,Y)                                             *GRW
*
      CALL GENMOV(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE GRWDRA(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      REAL X,Y
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      X=REAL(XX)
      Y=REAL(YY)
*GRW  CALL GWLINE2(IRTN,X,Y)                                            *GRW
*
      CALL GENDRA(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE GRWDOT(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      REAL X,Y,W,H,XH,YH
      SAVE IC,W,H,XH,YH
      DATA IC / 0 /
      LOGICAL LGRDOT
*
      IF(IC.EQ.0) THEN
         IC = 1
*GRW     CALL GWGETTXT(IRTN,W,H,XH,YH,'.')                               *GRW
      ENDIF
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      X=REAL(XX)
      Y=REAL(YY)
      IF(LGRDOT(XXX,YYY,ZZZ)) THEN
*        OUTPUTS A PERIOD SUCH THAT THE CENTER IS AT (X,Y)
*GRW     CALL GWPUTTXT(IRTN,X-W*0.55,Y-H*0.05,'.')                       *GRW
      ENDIF
*GRW  CALL GWMOVE2(IRTN,X,Y)                                             *GRW
*
      CALL GENDOT(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE GRWTRI(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      REAL TXY(2,4)
      LOGICAL LGRTRI
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      TXY(1,1)=REAL(XX)
      TXY(2,1)=REAL(YY)
      IF(LGRTRI(XXX,YYY,ZZZ,XP(1),XP(2),XP(3),XC(1),XC(2),XC(3))) THEN
         CALL PROJEC(XP(1),XP(2),XP(3),PXX,PYY)
         TXY(1,2)=REAL(PXX)
         TXY(2,2)=REAL(PYY)
         CALL PROJEC(XC(1),XC(2),XC(3),CXX,CYY)
         TXY(1,3)=REAL(CXX)
         TXY(2,3)=REAL(CYY)
         IF(LWR.EQ.1) THEN
            TXY(1,4)=REAL(XX)
            TXY(2,4)=REAL(YY)
*GRW        CALL GWPOLYLIN(IRTN,TXY,4)                                   *GRW
         ELSE
*GRW        CALL GWPOLYGON(IRTN,TXY,3,0)                                 *GRW
         ENDIF
      ENDIF
*GRW  CALL GWMOVE2(IRTN,TXY(1,1),TXY(2,1))                               *GRW
*
      CALL GENTRI(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE GRWPLY(IA,IPST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION OACC(4),CLR(4),PC(7)
      LOGICAL LGRPLY,LCR
*
      IF(IPST.EQ.-1) GOTO 90
      IF(.NOT.LGRPLY(IA)) GOTO 90
*
      CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
      IF(LCR) THEN
         ND = 6
         DO 10 K=1,4
         OACC(K) = ACC(K)
 10      CONTINUE
      ELSE
         ND = 3
      ENDIF
*
*     CURVE
*     *****
*
      IF(IOBJ.EQ.1) THEN
         CALL PLYEV(IA,ND,0.D0,NP1)
         CALL GRWMOV(PNT(1,0),PNT(2,0),PNT(3,0))
         DO 20 I=1,NP1
         DO 11 J=4,ND
 11      CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
         IF(LCR) CALL GRWCOL(CLR)
         CALL GRWDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 20      CONTINUE
*
*     SURFACE
*     *******
*
      ELSEIF(IOBJ.EQ.2) THEN
*
*      * VARYING OR SOLID COLOR TRIANGLE FILL PAINTING *
         IF(LWR.EQ.0) THEN
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            DO 40 I2=1,NP2
            DO 21 I=0,NP1
            DO 21 J=1,ND
 21         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 40 I=1,NP1
            DO 22 J=1,ND
 22         PC(J) = (PNTO(J,I-1)+PNTO(J,I)+PNT(J,I-1)+PNT(J,I))*0.25D0
*
            CALL GRWMOV(PNT(1,I),PNT(2,I),PNT(3,I))
            CALL GRWMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            DO 31 J=4,ND
 31         CLR(J-3) = (PNT(J,I)+PNTO(J,I)+PC(J))/3.D0
            IF(LCR) CALL GRWCOL(CLR)
            CALL GRWTRI(PC(1),PC(2),PC(3))
*
            DO 32 J=4,ND
 32         CLR(J-3) = (PNTO(J,I)+PC(J)+PNTO(J,I-1))/3.D0
            IF(LCR) CALL GRWCOL(CLR)
            CALL GRWTRI(PNTO(1,I-1),PNTO(2,I-1),PNTO(3,I-1))
*
            DO 33 J=4,ND
 33         CLR(J-3) = (PC(J)+PNTO(J,I-1)+PNT(J,I-1))/3.D0
            IF(LCR) CALL GRWCOL(CLR)
            CALL GRWTRI(PNT(1,I-1),PNT(2,I-1),PNT(3,I-1))
*
            CALL GRWMOV(PC(1),PC(2),PC(3))
            DO 34 J=4,ND
 34         CLR(J-3) = (PNT(J,I-1)+PC(J)+PNT(J,I))/3.D0
            IF(LCR) CALL GRWCOL(CLR)
            CALL GRWTRI(PNT(1,I),PNT(2,I),PNT(3,I))
 40         CONTINUE
*
*      * WIRE FRAME QUADRILATERALS *
         ELSE
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            CALL GRWMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 50 I=1,NP1
            DO 41 J=4,ND
 41         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL GRWCOL(CLR)
            CALL GRWDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 50         CONTINUE
*
            DO 80 I2=1,NP2
            DO 51 I=0,NP1
            DO 51 J=1,ND
 51         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 60 I=0,NP1
            DO 52 J=4,ND
 52         CLR(J-3) = (PNTO(J,I)+PNT(J,I))*0.5D0
            IF(LCR) CALL GRWCOL(CLR)
            CALL GRWMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            CALL GRWDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 60         CONTINUE
*
            CALL GRWMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 70 I=1,NP1
            DO 61 J=4,ND
 61         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL GRWCOL(CLR)
            CALL GRWDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 70         CONTINUE
 80         CONTINUE
         ENDIF
      ENDIF
*
      IF(LCR) CALL GRWCOL(OACC)
*
 90   CALL GENPLY(IA)
*
      RETURN
      END
*
      SUBROUTINE GRWCHA(STR,LST)
*     **************************
*
      CHARACTER STR*1024
      REAL X,Y
*
*GRW  CALL GWGETPOS(IRTN,X,Y)                                            *GRW
*GRW  CALL GWPUTTXT(IRTN,X,Y,STR(1:LST))                                 *GRW
*
      CALL GENCHA
*
      RETURN
      END
*
      SUBROUTINE GRWCOL(CLR)
*     **********************
*
      DOUBLE PRECISION CLR(4)
      INTEGER ICLR(3)
*
      DO 10 I=1,3
      ICLR(I) = INT(CLR(I)*255.D0)
 10   CONTINUE
*GRW  K = KRGB(ICLR(1),ICLR(2),ICLR(3))                                  *GRW
*GRW  CALL GWSETPEN(IRTN,K,1,-1,-1)                                      *GRW
*GRW  CALL GWSETBRS(IRTN,K,-1,-1)                                        *GRW
*GRW  CALL GWCOLOR(IRTN,K,7)                                             *GRW
*
      CALL GENCOL(CLR)
*
      RETURN
      END
*
      SUBROUTINE GRWWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*GRW  CALL GWSETPEN(IRTN,-1,1,MAX(NINT(5.D0*AJ),0),-1)                   *GRW
*
      CALL GENWID(AJ)
*
      RETURN
      END
*
      SUBROUTINE GRWEND(NARI)
*     ***********************
*
      CHARACTER STR*100
      DOUBLE PRECISION CLR(4)
*
      IF(NARI.GT.0) THEN
         DO 10 I=1,3
         CLR(I) = 0.D0
 10      CONTINUE
         CALL GRWCOL(CLR)
*
         WRITE(STR,'(A,I0)') '$$$ WARNING, GRPOLY ARITHMETIC FAILURES: '
     *             ,NARI
         L = ILAST(STR,1,100)
*GRW     CALL GWPUTTXT(IRTN,0.007,0.01,STR(1:L))                         *GRW
      ENDIF
*
      CALL GENEND
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     AQUATERM DRIVER ROUTINES                                                *
*     CONTRIBUTED BY ALEXANDER WITTIG (2006,2012)                             *
*     UPDATED BY K. MAKINO                                                    *
*                                                                             *
*     USE OF AQUATERM KINDLY ACKNOWLEDGED, SOURCES AVAILABLE AT               *
*     http://aquaterm.sf.net/                                                 *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE AQTBEG(IUNIT)
*     ************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(AQTW=1080.D0,AQTH=720.D0)
      CHARACTER TITLE*3
      INTEGER IUNIT
      LOGICAL FIRST
      DATA FIRST / .TRUE. /
      SAVE FIRST
*
      IF(FIRST) THEN
*AQT     CALL AQTINIT                                                    *AQT
         FIRST = .FALSE.
      ENDIF
*AQT  CALL AQTOPENPLOT(ABS(IUNIT))                                       *AQT
*AQT  CALL AQTSETPLOTSIZE(REAL(AQTW),REAL(AQTH))                         *AQT
      WRITE(TITLE,'(I3)') ABS(IUNIT)
*AQT  CALL AQTSETPLOTTITLE('COSY PLOT #'//TITLE)                         *AQT
*AQT  CALL AQTSETCOLOR(0.0,0.0,0.0)                                      *AQT
*AQT  CALL AQTSETBACKGROUNDCOLOR(1.0,1.0,1.0)                            *AQT
*AQT  CALL AQTSETLINECAPSTYLE(1)                                         *AQT
*AQT  CALL AQTSETFONTNAME('Helvetica')                                   *AQT
*AQT  CALL AQTSETFONTSIZE(21.0)                                          *AQT
*AQT  CALL AQTSETLINEWIDTH(1.3)                                          *AQT
*
      CALL GENBEG
*
      RETURN
      END
*
      SUBROUTINE AQTMOV(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(AQTW=1080.D0,AQTH=720.D0)
      REAL X,Y
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      X = REAL(XX*AQTW)
      Y = REAL(YY*AQTH)
*AQT  CALL AQTMOVETO(X,Y)                                                *AQT
*
      CALL GENMOV(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE AQTDRA(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      PARAMETER(AQTW=1080.D0,AQTH=720.D0)
      REAL X,Y
*
*     THIS IS A BUG/FEATURE OF AQUATERM: IT FORGETS THE CURRENT POSITION
*     AFTER CHANGING COLOR, WIDTH, ETC. FIX BY ADDING AN EXPLICIT MOVE
*     BEFORE EACH DRAW UNLESS IT WAS IMMEDIATELY PRECEEDED BY ANOTHER DRAW.
      IF(LCM.NE.1) CALL AQTMOV(XC(1),XC(2),XC(3))
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      X = REAL(XX*AQTW)
      Y = REAL(YY*AQTH)
*AQT  CALL AQTADDLINETO(X,Y)                                             *AQT
*
      CALL GENDRA(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE AQTDOT(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(AQTW=1080.D0,AQTH=720.D0)
      REAL X,Y
      LOGICAL LGRDOT
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      X = REAL(XX*AQTW)
      Y = REAL(YY*AQTH)
      IF(LGRDOT(XXX,YYY,ZZZ)) THEN
*        THERE'S NO CIRCLE, ELLIPSE, OR CURVE COMMAND IN AQUATERM.
*        DRAWING A SHORT LINE WITH ROUND LINE CAPS PRODUCES A FILLED DOT.
*AQT     CALL AQTMOVETO(X-0.25,Y)                                        *AQT
*AQT     CALL AQTADDLINETO(X+0.25,Y)                                     *AQT
      ENDIF
*AQT  CALL AQTMOVETO(X,Y)                                                *AQT
*
      CALL GENDOT(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE AQTTRI(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      PARAMETER(AQTW=1080.D0,AQTH=720.D0)
      REAL TX(4),TY(4)
      LOGICAL LGRTRI
*
      CALL PROJEC(XXX,YYY,ZZZ,XX,YY)
      TX(1) = REAL(XX*AQTW)
      TY(1) = REAL(YY*AQTH)
      IF(LGRTRI(XXX,YYY,ZZZ,XP(1),XP(2),XP(3),XC(1),XC(2),XC(3))) THEN
         CALL PROJEC(XP(1),XP(2),XP(3),PXX,PYY)
         TX(2) = REAL(PXX*AQTW)
         TY(2) = REAL(PYY*AQTH)
         CALL PROJEC(XC(1),XC(2),XC(3),CXX,CYY)
         TX(3) = REAL(CXX*AQTW)
         TY(3) = REAL(CYY*AQTH)
         IF(LWR.EQ.1) THEN
            TX(4) = TX(1)
            TY(4) = TY(1)
*AQT        CALL AQTADDPOLYLINE(TX,TY,4)                                 *AQT
         ELSE
*AQT        CALL AQTADDPOLYGON(TX,TY,3)                                  *AQT
         ENDIF
      ENDIF
*AQT  CALL AQTMOVETO(TX(1),TY(1))                                        *AQT
*
      CALL GENTRI(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE AQTPLY(IA,IPST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION OACC(4),CLR(4),PC(7)
      LOGICAL LGRPLY,LCR
*
      IF(IPST.EQ.-1) GOTO 90
      IF(.NOT.LGRPLY(IA)) GOTO 90
*
      CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
      IF(LCR) THEN
*        FOR RGB, ND=6. IF ALPHA FOR RGBA IS AVAILABLE, SET ND=7.
         ND = 6
         DO 10 K=1,4
         OACC(K) = ACC(K)
 10      CONTINUE
      ELSE
         ND = 3
      ENDIF
*
*     CURVE
*     *****
*
      IF(IOBJ.EQ.1) THEN
         CALL PLYEV(IA,ND,0.D0,NP1)
         CALL AQTMOV(PNT(1,0),PNT(2,0),PNT(3,0))
         DO 20 I=1,NP1
         DO 11 J=4,ND
 11      CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
         IF(LCR) CALL AQTCOL(CLR)
         CALL AQTDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 20      CONTINUE
*
*     SURFACE
*     *******
*
      ELSEIF(IOBJ.EQ.2) THEN
*
*      * VARYING OR SOLID COLOR TRIANGLE FILL PAINTING *
         IF(LWR.EQ.0) THEN
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            DO 40 I2=1,NP2
            DO 21 I=0,NP1
            DO 21 J=1,ND
 21         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 40 I=1,NP1
            DO 22 J=1,ND
 22         PC(J) = (PNTO(J,I-1)+PNTO(J,I)+PNT(J,I-1)+PNT(J,I))*0.25D0
*
            CALL AQTMOV(PNT(1,I),PNT(2,I),PNT(3,I))
            CALL AQTMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            DO 31 J=4,ND
 31         CLR(J-3) = (PNT(J,I)+PNTO(J,I)+PC(J))/3.D0
            IF(LCR) CALL AQTCOL(CLR)
            CALL AQTTRI(PC(1),PC(2),PC(3))
*
            DO 32 J=4,ND
 32         CLR(J-3) = (PNTO(J,I)+PC(J)+PNTO(J,I-1))/3.D0
            IF(LCR) CALL AQTCOL(CLR)
            CALL AQTTRI(PNTO(1,I-1),PNTO(2,I-1),PNTO(3,I-1))
*
            DO 33 J=4,ND
 33         CLR(J-3) = (PC(J)+PNTO(J,I-1)+PNT(J,I-1))/3.D0
            IF(LCR) CALL AQTCOL(CLR)
            CALL AQTTRI(PNT(1,I-1),PNT(2,I-1),PNT(3,I-1))
*
            CALL AQTMOV(PC(1),PC(2),PC(3))
            DO 34 J=4,ND
 34         CLR(J-3) = (PNT(J,I-1)+PC(J)+PNT(J,I))/3.D0
            IF(LCR) CALL AQTCOL(CLR)
            CALL AQTTRI(PNT(1,I),PNT(2,I),PNT(3,I))
 40         CONTINUE
*
*      * WIRE FRAME QUADRILATERALS *
         ELSE
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            CALL AQTMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 50 I=1,NP1
            DO 41 J=4,ND
 41         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL AQTCOL(CLR)
            CALL AQTDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 50         CONTINUE
*
            DO 80 I2=1,NP2
            DO 51 I=0,NP1
            DO 51 J=1,ND
 51         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 60 I=0,NP1
            DO 52 J=4,ND
 52         CLR(J-3) = (PNTO(J,I)+PNT(J,I))*0.5D0
            IF(LCR) CALL AQTCOL(CLR)
            CALL AQTMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            CALL AQTDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 60         CONTINUE
*
            CALL AQTMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 70 I=1,NP1
            DO 61 J=4,ND
 61         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL AQTCOL(CLR)
            CALL AQTDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 70         CONTINUE
 80         CONTINUE
         ENDIF
      ENDIF
*
      IF(LCR) CALL AQTCOL(OACC)
*
 90   CALL GENPLY(IA)
*
      RETURN
      END
*
      SUBROUTINE AQTCHA(STR,LST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      PARAMETER(AQTW=1080.D0,AQTH=720.D0)
      REAL X,Y
      CHARACTER STR*1024
*
      CALL PROJEC(XC(1),XC(2),XC(3),XX,YY)
      X = REAL(XX*AQTW)
      Y = REAL(YY*AQTH)
*     THE LAST PARAMETER IS THE VERTICAL TEXT ALIGNMENT
*     0 - MIDDLE, 4 - BASELINE, 8 - BOTTOM, 16 - TOP
*AQT  CALL AQTADDLABEL(STR(1:LST),X,Y,0,8)                               *AQT
*
      CALL GENCHA
*
      RETURN
      END
*
      SUBROUTINE AQTCOL(CLR)
*     **********************
*
      DOUBLE PRECISION CLR(4)
*
*AQT  CALL AQTSETCOLOR(REAL(CLR(1)),REAL(CLR(2)),REAL(CLR(3)))           *AQT
*
*     NOTE: IF ALPHA FOR RGBA IS AVAILABLE,
*           REPLACE THE ABOVE AQTSETCOLOR CALL BY AQTSETCOLORRGBA?
*           ALSO, SET ND=7 IN AQTPLY.
*
      CALL GENCOL(CLR)
*
      RETURN
      END
*
      SUBROUTINE AQTWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*AQT  CALL AQTSETLINEWIDTH(1.3*REAL(AJ))                                 *AQT
*
      CALL GENWID(AJ)
*
      RETURN
      END
*
      SUBROUTINE AQTEND(NARI)
*     ***********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER(AQTW=1080.D0,AQTH=720.D0)
      REAL X,Y
      CHARACTER STR*100
*
      IF(NARI.GT.0) THEN
*AQT     CALL AQTSETCOLOR(0.0,0.0,0.0)                                   *AQT
         WRITE(STR,'(A,I0)') '$$$ WARNING, GRPOLY ARITHMETIC FAILURES: '
     *             ,NARI
         L = ILAST(STR,1,100)
         XX = 0.007D0
         YY = 0.010D0
         X = REAL(XX*AQTW)
         Y = REAL(YY*AQTH)
*AQT     CALL AQTADDLABEL(STR(1:L),X,Y,0,8)                              *AQT
      ENDIF
*
*AQT  CALL AQTRENDERPLOT                                                 *AQT
*
      CALL GENEND
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     POSTSCRIPT DRIVER ROUTINES                                              *
*                                                                             *
*******************************************************************************
*     EDIT DESCRIPTORS Iw, Fw.d: IF w=0, THE PROCESSOR SELECTS THE FIELD WIDTH.
*
      SUBROUTINE POSESC(STR,LST,STRE,LSTE)
*     ************************************
*
*     THIS SUBROUTINE ESCAPES A STRING TO BE INCLUDED IN POSTSCRIPT OR PDF
*     FROM STR (LENGTH LST) TO STRE (LENGTH LSTE)
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
      CHARACTER STR*(*),STRE*(*)
*
      J = 0
      DO 10 I=1,LST
      IF(STR(I:I).EQ.SCBSL.OR.STR(I:I).EQ.'('.OR.STR(I:I).EQ.')') THEN
         J = J+1
         STRE(J:J) = SCBSL
      ENDIF
      J = J+1
      STRE(J:J) = STR(I:I)
      IF(J+2.GT.LEN(STRE)) GOTO 20
 10   CONTINUE
 20   LSTE = J
*
      RETURN
      END
*
      SUBROUTINE POSBEG
*     *****************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH FILE OUTPUT -----------------------------------------------------
      INTEGER ICNT
      CHARACTER FBASE*20
      COMMON /GRFILE/ ICNT,FBASE
*----------------------------------------------------------------------------
*
      CHARACTER NUMBER*3
*
      WRITE(NUMBER,'(I3.3)') MOD(ICNT,1000)
      IL = ILAST(FBASE,1,20)
      OPEN(4,FILE=FBASE(1:IL)//NUMBER//'.ps',STATUS='UNKNOWN')
      ICNT = ICNT+1
*
      WRITE(4,'(A)') '%!PS-Adobe-3.0'
      WRITE(4,'(A)') '%%Creator: COSY INFINITY Graphics'
      WRITE(4,'(A)') '%%Title: '//FBASE(1:IL)//NUMBER//'.ps'
      WRITE(4,'(A)') '%%BoundingBox: 50 21 550 771'
      WRITE(4,'(A)') '%%Orientation: Landscape'
      WRITE(4,'(A)') '%%Pages: 1'
      WRITE(4,'(A)') '%%Page: 1 1'
      WRITE(4,'(A)') '/sn {stroke newpath} bind def'
      WRITE(4,'(A)') '/w {setlinewidth} def'
      WRITE(4,'(A)') '/m {moveto} bind def'
      WRITE(4,'(A)') '/l {lineto} bind def'
      WRITE(4,'(A)') '/ms {newpath moveto gsave 1 1.5 scale show'
     *   //' grestore} bind def'
      WRITE(4,'(A)') '/t {moveto lineto lineto closepath gsave fill'
     *   //' grestore stroke newpath} bind def'
      WRITE(4,'(A)') '/tw {moveto lineto lineto closepath stroke'
     *   //' newpath} bind def'
      WRITE(4,'(A)') '/c {setrgbcolor} bind def % - color picture'
      WRITE(4,'(A)') '%/c {add add 3 div setgray} bind def'
     *   //' % - grayscale picture except for shfill GRPOLY objects'
      WRITE(4,'(A)') '/Helvetica findfont 0.02 scalefont setfont'
      WRITE(4,'(A)') '500 750 scale 1.100 0.0280 translate 90 rotate'
      WRITE(4,'(A)') '1 setlinejoin 1 setlinecap'
      WRITE(4,'(A)') '0.0015 setlinewidth 0 0 1 1 rectclip'
*
      CALL GENBEG
*
      RETURN
      END
*
      SUBROUTINE POSMOV(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      IF(LCM.EQ.1) WRITE(4,'(A)') 'sn'
*
      CALL GENMOV(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE POSDRA(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      IF(LCM.NE.1) THEN
         CALL PROJEC(XC(1),XC(2),XC(3),CX,CY)
         WRITE(4,'(2(F0.9,X),A)') CX,CY,'m'
      ENDIF
      CALL PROJEC(XXX,YYY,ZZZ,X,Y)
      WRITE(4,'(2(F0.9,X),A)') X,Y,'l'
*
      CALL GENDRA(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE POSDOT(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      LOGICAL LGRDOT
      SAVE SIZE,SX,SY
*
      IF(LCM.EQ.1) WRITE(4,'(A)') 'sn'
*
      IF(LGRDOT(XXX,YYY,ZZZ)) THEN
         IF(LDOTW.EQ.1) THEN
            LDOTW = -1
            SIZE = 0.0052D0*ACW
            SX = -0.00091D0*ACW
            SY = -0.00266D0*ACW
         ENDIF
*
         IF(LDOTC.EQ.0.OR.LDOTW.EQ.-1) THEN
            LDOTW = 0
            LDOTC = 1
            WRITE(4,'(A,F0.9,A)') '/Helvetica findfont ',SIZE,
     *              ' scalefont setfont'
         ENDIF
*
         CALL PROJEC(XXX,YYY,ZZZ,X,Y)
         WRITE(4,'(A,2(F0.9,X),A)') '<B7> ',X+SX,Y+SY,'ms'
      ENDIF
*
      CALL GENDOT(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE POSTRI(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      LOGICAL LGRTRI
*
      IF(LCM.EQ.1) WRITE(4,'(A)') 'sn'
*
      IF(LGRTRI(XXX,YYY,ZZZ,XP(1),XP(2),XP(3),XC(1),XC(2),XC(3))) THEN
         CALL PROJEC(XXX,YYY,ZZZ,X,Y)
         CALL PROJEC(XP(1),XP(2),XP(3),PX,PY)
         CALL PROJEC(XC(1),XC(2),XC(3),CX,CY)
         IF(LWR.EQ.1) THEN
            WRITE(4,'(6(F0.9,X),A)') PX,PY,X,Y,CX,CY,'tw'
         ELSE
            WRITE(4,'(6(F0.9,X),A)') PX,PY,X,Y,CX,CY,'t'
         ENDIF
      ENDIF
*
      CALL GENTRI(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE POSPLY(IA,IPST)
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
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      PARAMETER(KAP=4)
      DOUBLE PRECISION OACC(4),CLR(4)
      INTEGER IX(4),JX(2,2)
      LOGICAL LGRPLY,LCR
*
      IF(LCM.EQ.1) WRITE(4,'(A)') 'sn'
*
      IPLYC = IPLY+1
      IF(IPST.EQ.-1) THEN
         WRITE(4,'(A,I0,5X,A)') '% POLY# ',IPLYC,'arithmetic failure'
         GOTO 90
      ENDIF
      IF(.NOT.LGRPLY(IA)) THEN
         WRITE(4,'(A,I0,5X,A)') '% POLY# ',IPLYC,'out of bound'
         GOTO 90
      ENDIF
*
      CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
      CALL PLYCMT(IA,IOBJ,NP1,NP2,LCR,-10)
*
      IF(LCR) THEN
         ND = 6
         DO 10 K=1,4
         OACC(K) = ACC(K)
 10      CONTINUE
      ELSE
         ND = 3
      ENDIF
*
*     CURVE
*     *****
*
      IF(IOBJ.EQ.1) THEN
         CALL PLYEV(IA,ND,0.D0,NP1)
*
*      * VARYING COLOR CURVE *
         IF(LCR) THEN
            CALL CURVDT0(ND,NP1,N)
*
            MLEN = (KAP+2)*(N+2)
            CALL FOXALL(IX,4,MLEN)
            JX(1,1) = NBEG(IX(1))
            JX(2,1) = NBEG(IX(2))
            JX(1,2) = NBEG(IX(3))
            JX(2,2) = NBEG(IX(4))
*
            CALL CURVDT(N,NCV,JX,KAP)
*
*           THIS MEMORY SIZE ERROR DOES NOT HAPPEN. CHECKING TO BE SURE.
            IF(NCV+1.GT.MLEN) CALL FOXERV('POSPLY')
*
            WRITE(4,96) 2
            DO 20 I=0,NCV
            DO 11 J=4,ND
 11         CLR(J-3) = PNT(J,NC(JX(1,1)+I))
            DO 20 LR=1,2
            WRITE(4,95) (CC(JX(J,LR)+I),J=1,2),(CLR(J),J=1,3)
 20         CONTINUE
            WRITE(4,'(A)') ']>>shfill'
*
            CALL FOXDAL(IX,4)
*
*      * SOLID COLOR CURVE *
         ELSE
            CALL POSMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 30 I=1,NP1
            CALL POSDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 30         CONTINUE
         ENDIF
*
*     SURFACE
*     *******
*
      ELSEIF(IOBJ.EQ.2) THEN
*
*      * VARYING OR SOLID COLOR LATTICE FILL PAINTING *
         IF(LWR.EQ.0) THEN
            WRITE(4,96) NP1+1
            DO 40 I2=0,NP2
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
            CALL PROJPNT(NP1)
*
            DO 40 I1=0,NP1
            IF(LCR) THEN
               WRITE(4,95) (PNT(J,I1),J=1,2),(PNT(J,I1),J=4,6)
            ELSE
               WRITE(4,95) (PNT(J,I1),J=1,2),(ACC(J),J=1,3)
            ENDIF
 40         CONTINUE
            WRITE(4,'(A)') ']>>shfill'
*
*      * WIRE FRAME QUADRILATERALS *
         ELSE
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            CALL POSMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 50 I=1,NP1
            DO 41 J=4,ND
 41         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL POSCOL(CLR)
            CALL POSDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 50         CONTINUE
*
            DO 80 I2=1,NP2
            DO 51 I=0,NP1
            DO 51 J=1,ND
 51         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 60 I=0,NP1
            DO 52 J=4,ND
 52         CLR(J-3) = (PNTO(J,I)+PNT(J,I))*0.5D0
            IF(LCR) CALL POSCOL(CLR)
            CALL POSMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            CALL POSDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 60         CONTINUE
*
            CALL POSMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 70 I=1,NP1
            DO 61 J=4,ND
 61         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL POSCOL(CLR)
            CALL POSDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 70         CONTINUE
 80         CONTINUE
         ENDIF
      ENDIF
*
      WRITE(4,'(A)') 'sn'
      IF(LCR) CALL POSCOL(OACC)
 90   CALL GENPLY(IA)
*
      RETURN
 95   FORMAT(2(F0.9,1X),3(F0.6,1X))
 96   FORMAT('<</ShadingType 5 /ColorSpace /DeviceRGB /VerticesPerRow '
     *       ,I0,' /DataSource [')
 196  FORMAT('<</ShadingType 5 /ColorSpace /DeviceGray /VerticesPerRow '
     *       ,I0,' /DataSource [')
*
*     TO TURN shfill FOR GRAYSCALE, USE 196 FORMAT INSTEAD OF 96 FORMAT,
*     AND GIVE ONLY ONE VALUE FOR THE COLOR PART.
*
      END
*
      SUBROUTINE POSPLY0(IA,IPST)
*     ***************************
*
*     THIS SUBROUTINE IS A VERSION OF POSPLY USING OTHER PRIMITIVE POS ROUTINES,
*     WITHOUT USING POSTSCRIPT'S shfill FOR SMOOTH COLOR TRANSITION.
*     THIS CAN BE USED AS A TEMPLATE FOR OTHER GRAPHICS SYSTEMS
*     WHEN SMOOTH COLOR TRANSITION IS NOT AVAILABLE.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION OACC(4),CLR(4),PC(7)
      LOGICAL LGRPLY,LCR
*
      IF(LCM.EQ.1) WRITE(4,'(A)') 'sn'
*
      IPLYC = IPLY+1
      IF(IPST.EQ.-1) THEN
         WRITE(4,'(A,I0,5X,A)') '% POLY# ',IPLYC,'arithmetic failure'
         GOTO 90
      ENDIF
      IF(.NOT.LGRPLY(IA)) THEN
         WRITE(4,'(A,I0,5X,A)') '% POLY# ',IPLYC,'out of bound'
         GOTO 90
      ENDIF
*
      CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
      CALL PLYCMT(IA,IOBJ,NP1,NP2,LCR,-10)
*
      IF(LCR) THEN
         ND = 6
         DO 10 K=1,4
         OACC(K) = ACC(K)
 10      CONTINUE
      ELSE
         ND = 3
      ENDIF
*
*     CURVE
*     *****
*
      IF(IOBJ.EQ.1) THEN
         CALL PLYEV(IA,ND,0.D0,NP1)
         CALL POSMOV(PNT(1,0),PNT(2,0),PNT(3,0))
         DO 20 I=1,NP1
         DO 11 J=4,ND
 11      CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
         IF(LCR) CALL POSCOL(CLR)
         CALL POSDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 20      CONTINUE
*
*     SURFACE
*     *******
*
      ELSEIF(IOBJ.EQ.2) THEN
*
*      * VARYING OR SOLID COLOR TRIANGLE FILL PAINTING *
         IF(LWR.EQ.0) THEN
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            DO 40 I2=1,NP2
            DO 21 I=0,NP1
            DO 21 J=1,ND
 21         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 40 I=1,NP1
            DO 22 J=1,ND
 22         PC(J) = (PNTO(J,I-1)+PNTO(J,I)+PNT(J,I-1)+PNT(J,I))*0.25D0
*
            CALL POSMOV(PNT(1,I),PNT(2,I),PNT(3,I))
            CALL POSMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            DO 31 J=4,ND
 31         CLR(J-3) = (PNT(J,I)+PNTO(J,I)+PC(J))/3.D0
            IF(LCR) CALL POSCOL(CLR)
            CALL POSTRI(PC(1),PC(2),PC(3))
*
            DO 32 J=4,ND
 32         CLR(J-3) = (PNTO(J,I)+PC(J)+PNTO(J,I-1))/3.D0
            IF(LCR) CALL POSCOL(CLR)
            CALL POSTRI(PNTO(1,I-1),PNTO(2,I-1),PNTO(3,I-1))
*
            DO 33 J=4,ND
 33         CLR(J-3) = (PC(J)+PNTO(J,I-1)+PNT(J,I-1))/3.D0
            IF(LCR) CALL POSCOL(CLR)
            CALL POSTRI(PNT(1,I-1),PNT(2,I-1),PNT(3,I-1))
*
            CALL POSMOV(PC(1),PC(2),PC(3))
            DO 34 J=4,ND
 34         CLR(J-3) = (PNT(J,I-1)+PC(J)+PNT(J,I))/3.D0
            IF(LCR) CALL POSCOL(CLR)
            CALL POSTRI(PNT(1,I),PNT(2,I),PNT(3,I))
 40         CONTINUE
*
*      * WIRE FRAME QUADRILATERALS *
         ELSE
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            CALL POSMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 50 I=1,NP1
            DO 41 J=4,ND
 41         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL POSCOL(CLR)
            CALL POSDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 50         CONTINUE
*
            DO 80 I2=1,NP2
            DO 51 I=0,NP1
            DO 51 J=1,ND
 51         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 60 I=0,NP1
            DO 52 J=4,ND
 52         CLR(J-3) = (PNTO(J,I)+PNT(J,I))*0.5D0
            IF(LCR) CALL POSCOL(CLR)
            CALL POSMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            CALL POSDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 60         CONTINUE
*
            CALL POSMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 70 I=1,NP1
            DO 61 J=4,ND
 61         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL POSCOL(CLR)
            CALL POSDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 70         CONTINUE
 80         CONTINUE
         ENDIF
      ENDIF
*
      WRITE(4,'(A)') 'sn'
      IF(LCR) CALL POSCOL(OACC)
*
 90   CALL GENPLY(IA)
*
      RETURN
      END
*
      SUBROUTINE POSCHA(STR,LST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STR*1024,STRE*1024
*
      IF(LCM.EQ.1) WRITE(4,'(A)') 'sn'
      IF(LDOTC.EQ.1)
     *   WRITE(4,'(A)') '/Helvetica findfont 0.02 scalefont setfont'
*
      CALL PROJEC(XC(1),XC(2),XC(3),X,Y)
      CALL POSESC(STR,LST,STRE,LSTE)
      WRITE(4,'(A,2(F0.9,X),A)') '('//STRE(1:LSTE)//') ',X,Y,'ms'
*
      CALL GENCHA
*
      RETURN
      END
*
      SUBROUTINE POSCOL(CLR)
*     **********************
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CLR(4)
*
      IF(LCM.EQ.1) WRITE(4,'(A)') 'sn'
*
      WRITE(4,'(3(F0.6,X),A)') (CLR(J),J=1,3),'c'
*
      CALL GENCOL(CLR)
*
      RETURN
      END
*
      SUBROUTINE POSWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      IF(LCM.EQ.1) WRITE(4,'(A)') 'sn'
*
      WRITE(4,'(F0.9,A)') 0.0015D0*MAX(0.5D0,AJ),' w'
*
      CALL GENWID(AJ)
*
      RETURN
      END
*
      SUBROUTINE POSEND(NARI)
*     ***********************
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STR*100
*
      IF(NARI.EQ.0) THEN
         IF(LCM.EQ.1) WRITE(4,'(A)') 'stroke'
      ELSE
         IF(LCM.EQ.1) WRITE(4,'(A)') 'sn'
         WRITE(4,'(3(F0.6,X),A)') 0.,0.,0.,'c'
         WRITE(STR,'(A,I0,A,2(F0.9,X),A)') '($$$ WARNING, GRPOLY '
     *      //'ARITHMETIC FAILURES: ',NARI,') ',0.007,0.01,'ms'
         WRITE(4,'(A)') STR(1:ILAST(STR,1,100))
      ENDIF
*
      WRITE(4,'(A)') 'showpage'
      CALL ENDCMT(NARI,-10)
      WRITE(4,'(A)') '%%EOF'
      CLOSE(4)
*
      CALL GENEND
*
      RETURN
      END
*
      SUBROUTINE PLYCMT(IA,IOBJ,NP1,NP2,LCR,IUNIT)
*     ********************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      PARAMETER(ERTOL=0.005D0)
      DOUBLE PRECISION XE(3)
      CHARACTER STA(3)*100,STRE*21
      DATA STRE / '% eps_xyz, eps_color:' /
      LOGICAL LCR,LAUTO
*
      IPLYC = IPLY+1
      CALL PLYPOS(IA,XE)
      WRITE(STA(1),'(A,I0,5X,3(A,E16.8),A)') '% POLY# ',IPLYC,
     *          'end point (',XE(1),',',XE(2),',',XE(3),' )'
*
      IF(IOBJ.EQ.1) THEN
         WRITE(STA(2),'(A)') '% curve,'
      ELSEIF(IOBJ.EQ.2) THEN
         WRITE(STA(2),'(A)') '% surface,'
      ENDIF
      IL = ILAST(STA(2),1,100)
*
      IF(LCR) THEN
         WRITE(STA(2)(IL+1:),'(A)') ' color polynomial,'
      ELSE
         WRITE(STA(2)(IL+1:),'(A)') ' solid color,'
      ENDIF
      IL = ILAST(STA(2),1,100)
*
      IF(IOBJ.EQ.2) THEN
         IF(LWR.EQ.0) THEN
            WRITE(STA(2)(IL+1:),'(A)') ' fill painting,'
         ELSE
            WRITE(STA(2)(IL+1:),'(A)') ' wire frame,'
         ENDIF
      ENDIF
      IL = ILAST(STA(2),1,100)
*
      WRITE(STA(2)(IL+1:),'(A,I0)') ' mesh ',NP1
      IL = ILAST(STA(2),1,100)
*
      IF(IOBJ.EQ.2) THEN
         WRITE(STA(2)(IL+1:),'(A,I0)') 'x',NP2
         IL = ILAST(STA(2),1,100)
      ENDIF
*
      NCP = NC(IA+1)
      INP = ABS(NCP)-IOBJ*1000000
      LAUTO = NCP.LT.0.OR.INP.EQ.0
*
      IF(LAUTO) THEN
         WRITE(STA(2)(IL+1:),'(A)') ' (automatic)'
      ELSE
         WRITE(STA(2)(IL+1:),'(A)') ' (user specified)'
      ENDIF
*
      IF(LAUTO) THEN
         STRE = '% eps_xyz, eps_color:'
         ERP = ERTOL*1.D2
         IF(EPSXYZ.EQ.0.D0.AND.EPSCLR.EQ.0.D0) THEN
            WRITE(STA(3),'(A,2(F4.1,A))') STRE,ERP,'%,',ERP,'%'
         ELSEIF(EPSXYZ.EQ.0.D0) THEN
            WRITE(STA(3),'(A,F4.1,A,E13.5)') STRE,ERP,'%,',EPSCLR
         ELSEIF(EPSCLR.EQ.0.D0) THEN
            WRITE(STA(3),'(A,E13.5,A,F4.1,A)') STRE,EPSXYZ,',',ERP,'%'
         ELSE
            WRITE(STA(3),'(2(A,E13.5))') STRE,EPSXYZ,',',EPSCLR
         ENDIF
      ENDIF
*
*   * POS
      IF(IUNIT.EQ.-10) THEN
         WRITE(4,'(A)') STA(1)(1:ILAST(STA(1),1,100))
         WRITE(4,'(A)') STA(2)(1:ILAST(STA(2),1,100))
         IF(LAUTO) WRITE(4,'(A)') STA(3)(1:ILAST(STA(3),1,100))
*
*   * PDF
      ELSEIF(IUNIT.EQ.-12) THEN
         CALL PDFW(STA(1),0,IE)
         CALL PDFW(STA(2),0,IE)
         IF(LAUTO) CALL PDFW(STA(3),0,IE)
*
*   * SVG
      ELSEIF(IUNIT.EQ.-13) THEN
         WRITE(4,'(A)') '<!--'
         WRITE(4,'(A)') STA(1)(3:ILAST(STA(1),1,100))
         WRITE(4,'(A)') STA(2)(3:ILAST(STA(2),1,100))
         IF(LAUTO) WRITE(4,'(A)') STA(3)(3:ILAST(STA(3),1,100))
         WRITE(4,'(A)') '-->'
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE ENDCMT(NARI,IUNIT)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STA(13)*100,STRX(3)*1
      DATA STRX / 'x','y','z' /
*
      WRITE(STA(1),'(A)') '% 3D bounds of the COSY GR object:'
      WRITE(STA(5),'(A)') '% 2D projection of the 3D bounds:'
      WRITE(STA(8),'(A,3(2X,F0.5))') '% Projection angles (deg):',
     *             (PRANG(J),J=1,3)
      WRITE(STA(9),'(A,I0)') '% Total number of GR polynomial objects: '
     *            ,NPLY+NARI
      IF(IZOOM.EQ.0) THEN
         NLINE = 10
         WRITE(STA(10),'(A)') '% Zoom bounds: none'
      ELSE
         NLINE = 13
         WRITE(STA(10),'(A)') '% Zoom bounds:'
      ENDIF
*
      DO 10 J=1,3
      WRITE(STA( 1+J),'(A,A,A,2(E16.8,X,A))')'%  ',STRX(J),' [',
     *                BDG(J,1),',',BDG(J,2),']'
      WRITE(STA(10+J),'(A,A,A,2(E16.8,X,A))')'%  ',STRX(J),' [',
     *                BD3(J,1),',',BD3(J,2),']'
 10   CONTINUE
*
      DO 20 J=1,2
      WRITE(STA( 5+J),'(A,A,A,2(E16.8,X,A))')'%  ',STRX(J),' [',
     *                BDG2(J,1),',',BDG2(J,2),']'
 20   CONTINUE
*
*   * POS
      IF(IUNIT.EQ.-10) THEN
         DO 30 I=1,NLINE
         WRITE(4,'(A)') STA(I)(1:ILAST(STA(I),1,100))
 30      CONTINUE
*
*   * PDF
      ELSEIF(IUNIT.EQ.-12) THEN
         DO 40 I=1,NLINE
         CALL PDFW(STA(I),0,IE)
 40      CONTINUE
*
*   * SVG
      ELSEIF(IUNIT.EQ.-13) THEN
         WRITE(4,'(A)') '<!--'
         DO 50 I=1,NLINE
         WRITE(4,'(A)') STA(I)(3:ILAST(STA(I),1,100))
 50      CONTINUE
         WRITE(4,'(A)') '-->'
      ENDIF
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     PDF DRIVER ROUTINES                                                     *
*     ORIGINALLY WRITTEN BY ALEXANDER WITTIG (2012)                           *
*     MODIFIED BY K. MAKINO (2013)                                            *
*     http://www.adobe.com/devnet/acrobat/pdfs/PDF32000_2008.pdf              *
*                                                                             *
*******************************************************************************
*     EDIT DESCRIPTORS Iw, Fw.d: IF w=0, THE PROCESSOR SELECTS THE FIELD WIDTH.
*
      SUBROUTINE PDFW(STR,INO,IE)
*     ***************************
*
*     THIS SUBROUTINE WRITES A STRING TERMINATED BY A LINE FEED CHARACTER
*     TO UNIT 4 AND KEEPS TRACK OF THE NUMBER OF CHARACTERS
*     (ASSUMED TO BE BYTES) WRITTEN.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----MEMORY MANAGEMENT -----------------------------------------------------
      PARAMETER(LMEM=140000000,LVAR=10000000,LDIM=1000)
      INTEGER NTYP(LVAR),NBEG(LVAR),NEND(LVAR),NMAX(LVAR),
     *        NC(LMEM),NDIM(LDIM)
      DOUBLE PRECISION CC(LMEM)
      COMMON        NTYP,NBEG,NEND,NMAX, CC,NC, NDIM,IDIM, IVAR,IMEM
      COMMON /TYID/ NRE,NST,NLO,NCM,NVE,NDA,NCD,NGR
*----------------------------------------------------------------------------
*-----PDF DRIVER ------------------------------------------------------------
      INTEGER NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
      COMMON /PDFCOM/ NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
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
      CHARACTER STR*(*),STRE*1201
*
      IF(INO.EQ.0) THEN
         L = ILAST(STR,1,LEN(STR))
         STRE(1:L+1) = STR(1:L)//SCLF
      ELSE
*        THIS MEMORY SIZE ERROR DOES NOT HAPPEN. CHECKING TO BE SURE.
         IF(NBEG(IXREF)+NOBJ.GT.NMAX(IXREF)) CALL FOXERV('PDFW')
         NC(NBEG(IXREF)+NOBJ) = NPOS
         NOBJ = NOBJ+1
         IF(STR.EQ.'obj') THEN
            WRITE(STRE,'(I0,A)') NOBJ,' 0 obj'
            L = ILAST(STRE,1,LEN(STRE))
            STRE(L+1:L+1) = SCLF
         ELSEIF(STR.EQ.'xref') THEN
            L = 4
            STRE(1:L+1) = STR(1:L)//SCLF
         ELSE
            PRINT*,'@@@ ERROR IN PDFW, WRONG KEY WORD'
            CALL FOXSTP(1)
         ENDIF
      ENDIF
*
      WRITE(4,'(A,$)',IOSTAT=IE,ERR=90) STRE(1:L+1)
      NPOS = NPOS+L+1
*
 90   RETURN
      END
*
      SUBROUTINE PDFBEG(INA)
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
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRAPH FILE OUTPUT -----------------------------------------------------
      INTEGER ICNT
      CHARACTER FBASE*20
      COMMON /GRFILE/ ICNT,FBASE
*----------------------------------------------------------------------------
*-----PDF DRIVER ------------------------------------------------------------
      INTEGER NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
      COMMON /PDFCOM/ NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
*----------------------------------------------------------------------------
*
      CHARACTER NUMBER*3,NUMIER*2,FNAME*27,STR*1000
*
*     OPENING A FILE
*
      WRITE(NUMBER,'(I3.3)') MOD(ICNT,1000)
      IL = ILAST(FBASE,1,20)
      FNAME = FBASE(1:IL)//NUMBER//'.pdf'
      ICNT = ICNT+1
*
      IER = 0
 10   OPEN(4,FILE=FNAME,STATUS='UNKNOWN')
      NPOS = 0
      CALL PDFW('%PDF-1.6',0,IE)
*
*     WHEN ACROBAT HAS THE FILE OPEN, IT LOCKS THE FILE AS "READONLY".
*     CAUSING INTEL "forrtl severe (47)" ERROR.
*     THE FOLLOWING OPERATON IS TO AVOID THIS PROBLEM.
      IF(IE.GT.0) THEN
         IER = IER+1
         IF(IER.EQ.1) PRINT*,'$$$ WARNING IN PDFBEG, FILE '//
     *             FNAME(1:ILAST(FNAME,1,27))//' LOCKED.'
         IF(IER.GE.100) THEN
            PRINT*,'$$$ ERROR IN PDFBEG,   FILES THROUGH '//
     *             FNAME(1:ILAST(FNAME,1,27))//' LOCKED.'
            CALL FOXDEB
         ENDIF
         NUMIER = '  '
         WRITE(NUMIER,'(I0)') IER
         FNAME = FBASE(1:IL)//NUMBER//'_'//
     *           NUMIER(1:ILAST(NUMIER,1,2))//'.pdf'
         GOTO 10
      ENDIF
      IF(IER.GT.0) PRINT*,'    INSTEAD WRITING TO FILE '//FNAME
*
*     INITIALIZATION
*
      CALL FOXALL(IXREF,1,NPLY+7)
      NOBJ = 0
      NSH = 0
      IF(NPLY.GE.1) CALL PDFSH(INA)
      ISH = 0
*
*     THE PDF MAIN GRAPHICS CONTENTS - THE (NSH+1)-TH PDF OBJECT
*
      CALL PDFW('obj',1,IE)
      WRITE(STR,'(A,I0,A)') '<</Length ',NSH+2,' 0 R>> stream'
      CALL PDFW(STR,0,IE)
      ISTREAM = NPOS
      CALL PDFW('750 0 0 500 0 0 cm 1 0 0 1 0.028 0.1 cm',0,IE)
      CALL PDFW('1 J 1 j 0.0015 w 0 0 1 1 re W n',0,IE)
      CALL PDFW('/F1 0.03 Tf 66.66 Tz',0,IE)
*
      CALL GENBEG
*
      RETURN
      END
*
      SUBROUTINE PDFMOV(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      IF(LCM.EQ.1) CALL PDFW('S',0,IE)
*
      CALL GENMOV(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE PDFDRA(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STR*100
*
      IF(LCM.NE.1) THEN
         CALL PROJEC(XC(1),XC(2),XC(3),CX,CY)
         WRITE(STR,'(2(F0.9,X),A)') CX,CY,'m'
         CALL PDFW(STR,0,IE)
      ENDIF
      CALL PROJEC(XXX,YYY,ZZZ,X,Y)
      WRITE(STR,'(2(F0.9,X),A)') X,Y,'l'
      CALL PDFW(STR,0,IE)
*
      CALL GENDRA(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE PDFDOT(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STR*100
      LOGICAL LGRDOT
      SAVE SIZE,SX,SY
*
      IF(LCM.EQ.1) CALL PDFW('S',0,IE)
*
      IF(LGRDOT(XXX,YYY,ZZZ)) THEN
         IF(LDOTW.EQ.1) THEN
            LDOTW = -1
            SIZE = 0.008D0*ACW
            SX = -0.00095D0*ACW
            SY = -0.00280D0*ACW+0.00041D0
         ENDIF
*
         IF(LDOTC.EQ.0.OR.LDOTW.EQ.-1) THEN
            LDOTW = 0
            LDOTC = 1
            WRITE(STR,'(A,F0.9,A)') '/F1 ',SIZE,' Tf'
            CALL PDFW(STR,0,IE)
         ENDIF
*
         CALL PROJEC(XXX,YYY,ZZZ,X,Y)
         WRITE(STR,'(A,2(F0.9,X),A)') 'BT ',X+SX,Y+SY,'Td <B7> Tj ET'
         CALL PDFW(STR,0,IE)
      ENDIF
*
      CALL GENDOT(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE PDFTRI(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STR*100
      LOGICAL LGRTRI
*
      IF(LCM.EQ.1) CALL PDFW('S',0,IE)
*
      IF(LGRTRI(XXX,YYY,ZZZ,XP(1),XP(2),XP(3),XC(1),XC(2),XC(3))) THEN
         CALL PROJEC(XXX,YYY,ZZZ,X,Y)
         CALL PROJEC(XP(1),XP(2),XP(3),PX,PY)
         CALL PROJEC(XC(1),XC(2),XC(3),CX,CY)
         IF(LWR.EQ.1) THEN
            WRITE(STR,'(3(2(F0.9,X),A))')CX,CY,'m ',X,Y,'l ',PX,PY,'l s'
         ELSE
            WRITE(STR,'(3(2(F0.9,X),A))')CX,CY,'m ',X,Y,'l ',PX,PY,'l b'
         ENDIF
         CALL PDFW(STR,0,IE)
      ENDIF
*
      CALL GENTRI(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE PDFSH(INA)
*     *********************
*
*     THIS SUBROUTINE OUTPUTS THE PDF TYPE 5 SHADING DICTIONARIES FOR
*     GRPOLY GRADIENT COLOR CURVES/SURFACES OF THE GRAPHICS OBJECT INA.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----PDF DRIVER ------------------------------------------------------------
      INTEGER NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
      COMMON /PDFCOM/ NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*-----PDF GRPOLY RELATED DATA -----------------------------------------------
      DOUBLE PRECISION PLYS(2,0:1),PLYDEC(2,2)
      INTEGER IPLYS,IPLYE
      COMMON /GRPDAT/ PLYS,PLYDEC,IPLYS,IPLYE
*----------------------------------------------------------------------------
*
*     LBX = 2^31-1 FOR '/BitsPerCoordinate 32' WITH ONE BIT LESS
*     2^31-1 IS THE LARGEST POSITIVE INTEGER IN STANDARD FORTRAN
*     LBC = 2^16-1 FOR '/BitsPerComponent 16'
      PARAMETER(LBX=2147483647,LBC=65535,KAP=4)
      INTEGER IX(4),JX(2,2),ICLR(3),IXY(2)
      CHARACTER STR*1000
      LOGICAL LGRPLY,LCR
*
*     SCALING FOR THE PDF 5 SHADING DATA TO KEEP COORDINATE VALUES IN [0,1]^2.
*     IPLYS -1: THE GRAPHICS IS OUTSIDE THE FRAME BD2F(,).
*            1: THE GRAPHICS IS, AT LEAST PARTIALLY, INSIDE THE FRAME.
*               SCALE THE J-DIRECTION AS X*PLYS(J,1)+PLYS(J,0)
*
      IPLYS = 0
      DO 10 J=1,2
      PLYS(J,0) = 0.D0
      PLYS(J,1) = 1.D0
 10   CONTINUE
*
      IF(IZOOM.EQ.1) THEN
         DO 20 J=1,2
         IF(BD2F(J,1).GT.BDG2(J,2).OR.BD2F(J,2).LT.BDG2(J,1)) THEN
            IPLYS = -1
            RETURN
         ELSEIF(BD2F(J,1).GT.BDG2(J,1).OR.BD2F(J,2).LT.BDG2(J,2)) THEN
            DB = BDG2(J,2)-BDG2(J,1)
*           INFLATE THE SIZE BY +-4%. SEE BD2F COMPUTATION IN GRPRE.
            PLYS(J,1) = (BD2F(J,2)-BD2F(J,1))/(DB*1.08D0)
            PLYS(J,0) = (BD2F(J,1)-(BDG2(J,1)-0.04D0*DB))/(DB*1.08D0)
         ENDIF
 20      CONTINUE
      ENDIF
*
      DO 30 J=1,2
      PLYDEC(J,1) = -PLYS(J,0)/PLYS(J,1)
      PLYDEC(J,2) = (2.0D0-PLYS(J,0))/PLYS(J,1)
 30   CONTINUE
*
*     SWEEPING THROUGH THE GRAPHICS OBJECT INA
*     ****************************************
*
      IA = NBEG(INA)+1
      NLEN = 0
 100  IA = IA+NLEN
      IF(IA.GT.IPLYE) RETURN
      NCI = NC(IA)/100000000
      NLEN = MOD(NC(IA),100000000)
      GOTO(200,100,110,120,130,130), NCI-4
      GOTO 100
*
*   * COLOR, WIDTH, OTHER POLY RELATED PROPERTIES
*
 110  CALL GENCOL(CC(IA))
      GOTO 100
 120  CALL GENWID(CC(IA))
      GOTO 100
 130  CALL GENSTL(NCI,NLEN,IA)
      GOTO 100
*
*   * POLY
*
 200  IPST = NINT(CC(IA))
      IF(IPST.EQ.-1) GOTO 100
      IF(.NOT.LGRPLY(IA)) GOTO 100
*
      CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
      IF(LCR) THEN
         ND = 6
      ELSE
         ND = 3
         DO 210 J=1,3
 210     ICLR(J) = INT(ACC(J)*LBC)
      ENDIF
*
*     CURVE * VARYING COLOR CURVE *
*     *****
*
      IF(IOBJ.EQ.1.AND.LCR) THEN
         CALL PLYEV(IA,ND,0.D0,NP1)
         CALL CURVDT0(ND,NP1,N)
*
         MLEN = (KAP+2)*(N+2)
         CALL FOXALL(IX,4,MLEN)
         JX(1,1) = NBEG(IX(1))
         JX(2,1) = NBEG(IX(2))
         JX(1,2) = NBEG(IX(3))
         JX(2,2) = NBEG(IX(4))
*
         CALL CURVDT(N,NCV,JX,KAP)
*
*        THIS MEMORY SIZE ERROR DOES NOT HAPPEN. CHECKING TO BE SURE.
         IF(NCV+1.GT.MLEN) CALL FOXERV('POSPLY')
*
         CALL PDFSHB(2,NCV+1)
         DO 240 I=0,NCV
         DO 220 J=4,ND
 220     ICLR(J-3) = INT(PNT(J,NC(JX(1,1)+I))*LBC)
         DO 240 LR=1,2
         DO 230 J=1,2
         SXY = CC(JX(J,LR)+I)*PLYS(J,1)+PLYS(J,0)
 230     IXY(J) = INT(SXY*LBX)
         WRITE(STR,'(2Z8.8,3Z4.4)') (IXY(J),J=1,2),(ICLR(J),J=1,3)
         CALL PDFW(STR,0,IE)
 240     CONTINUE
         CALL PDFW('endstream endobj',0,IE)
*
         CALL FOXDAL(IX,4)
*
*     SURFACE * VARYING OR SOLID COLOR LATTICE FILL PAINTING *
*     *******
*
      ELSEIF(IOBJ.EQ.2.AND.LWR.EQ.0) THEN
         CALL PDFSHB(NP1+1,NP2+1)
         DO 270 I2=0,NP2
         CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
         CALL PROJPNT(NP1)
         DO 270 I1=0,NP1
         IF(LCR) THEN
            DO 250 J=4,ND
 250        ICLR(J-3) = INT(PNT(J,I1)*LBC)
         ENDIF
         DO 260 J=1,2
         SXY = PNT(J,I1)*PLYS(J,1)+PLYS(J,0)
 260     IXY(J) = INT(SXY*LBX)
         WRITE(STR,'(2Z8.8,3Z4.4)') (IXY(J),J=1,2),(ICLR(J),J=1,3)
         CALL PDFW(STR,0,IE)
 270     CONTINUE
         CALL PDFW('endstream endobj',0,IE)
      ENDIF
*
      GOTO 100
*
      END
*
      SUBROUTINE PDFSHB(NVPR,NROWS)
*     *****************************
*
*     THE HEADER OF THE PDF TYPE 5 SHADING DICTIONARY
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----PDF DRIVER ------------------------------------------------------------
      INTEGER NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
      COMMON /PDFCOM/ NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
*----------------------------------------------------------------------------
*-----PDF GRPOLY RELATED DATA -----------------------------------------------
      DOUBLE PRECISION PLYS(2,0:1),PLYDEC(2,2)
      INTEGER IPLYS,IPLYE
      COMMON /GRPDAT/ PLYS,PLYDEC,IPLYS,IPLYE
*----------------------------------------------------------------------------
*
      CHARACTER STR*1000
*
      NSH = NSH+1
*
      CALL PDFW('obj',1,IE)
      WRITE(STR,'(A,I0)') '% /Sh',NSH
      CALL PDFW(STR,0,IE)
      WRITE(STR,'(A,I0)') '<</ShadingType 5 /ColorSpace /DeviceRGB '
     *        //'/VerticesPerRow ',NVPR
      CALL PDFW(STR,0,IE)
      WRITE(STR,'(A,4(F0.10,X),A)') '/Decode [',
     *          (PLYDEC(J,1),PLYDEC(J,2),J=1,2),'0 1 0 1 0 1]'
      CALL PDFW(STR,0,IE)
      CALL PDFW('/Filter /ASCIIHexDecode '
     *        //'/BitsPerCoordinate 32 /BitsPerComponent 16',0,IE)
*     2x(32/4 bytes) + 3x(16/4 bytes) + 1 = 29
      WRITE(STR,'(A,I0,A)') '/Length ',29*NVPR*NROWS-1,' >> stream'
      CALL PDFW(STR,0,IE)
*
*     GRAYSCALE VERSION -- USE ONLY ONE NUMBER FOR THE COLOR PART.
*     REPLACE THE RESPECTIVE LINES BY THE FOLLOWING.
C     WRITE(STR,'(A,I0)') '<</ShadingType 5 /ColorSpace /DeviceGray '
C    *          (PLYDEC(J,1),PLYDEC(J,2),J=1,2),'0 1]'
C*    2x(32/4 bytes) + 1x(16/4 bytes) + 1 = 21
C     WRITE(STR,'(A,I0,A)') '/Length ',21*NVPR*NROWS-1,' >> stream'
*
      RETURN
      END
*
      SUBROUTINE PDFPLY(IA,IPST)
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
*-----PDF DRIVER ------------------------------------------------------------
      INTEGER NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
      COMMON /PDFCOM/ NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
*----------------------------------------------------------------------------
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*-----PDF GRPOLY RELATED DATA -----------------------------------------------
      DOUBLE PRECISION PLYS(2,0:1),PLYDEC(2,2)
      INTEGER IPLYS,IPLYE
      COMMON /GRPDAT/ PLYS,PLYDEC,IPLYS,IPLYE
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION OACC(4),CLR(4)
      CHARACTER STR*1000
      LOGICAL LGRPLY,LCR
*
      IF(LCM.EQ.1) CALL PDFW('S',0,IE)
*
      IPLYC = IPLY+1
      IF(IPLYS.EQ.-1) THEN
         WRITE(STR,'(A,I0,5X,A)') '% POLY# ',IPLYC,'out of bound'
         CALL PDFW(STR,0,IE)
         GOTO 90
      ENDIF
      IF(IPST.EQ.-1) THEN
         WRITE(STR,'(A,I0,5X,A)') '% POLY# ',IPLYC,'arithmetic failure'
         CALL PDFW(STR,0,IE)
         GOTO 90
      ENDIF
      IF(.NOT.LGRPLY(IA)) THEN
         WRITE(STR,'(A,I0,5X,A)') '% POLY# ',IPLYC,'out of bound'
         CALL PDFW(STR,0,IE)
         GOTO 90
      ENDIF
*
      CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
      CALL PLYCMT(IA,IOBJ,NP1,NP2,LCR,-12)
*
      IF(LCR) THEN
         ND = 6
         DO 10 K=1,4
         OACC(K) = ACC(K)
 10      CONTINUE
      ELSE
         ND = 3
      ENDIF
*
*     CURVE
*     *****
*
      IF(IOBJ.EQ.1) THEN
*
*      * VARYING COLOR CURVE *
         IF(LCR) THEN
            ISH = ISH+1
            WRITE(STR,'(A,I0,A)') '/Sh',ISH,' sh'
            CALL PDFW(STR,0,IE)
*
*      * SOLID COLOR CURVE *
         ELSE
            CALL PLYEV(IA,ND,0.D0,NP1)
            CALL PDFMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 30 I=1,NP1
            CALL PDFDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 30         CONTINUE
         ENDIF
*
*     SURFACE
*     *******
*
      ELSEIF(IOBJ.EQ.2) THEN
*
*      * VARYING OR SOLID COLOR TRIANGLE FILL PAINTING *
         IF(LWR.EQ.0) THEN
            ISH = ISH+1
            WRITE(STR,'(A,I0,A)') '/Sh',ISH,' sh'
            CALL PDFW(STR,0,IE)
*
*      * WIRE FRAME QUADRILATERALS *
         ELSE
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            CALL PDFMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 50 I=1,NP1
            DO 41 J=4,ND
 41         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL PDFCOL(CLR)
            CALL PDFDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 50         CONTINUE
*
            DO 80 I2=1,NP2
            DO 51 I=0,NP1
            DO 51 J=1,ND
 51         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 60 I=0,NP1
            DO 52 J=4,ND
 52         CLR(J-3) = (PNTO(J,I)+PNT(J,I))*0.5D0
            IF(LCR) CALL PDFCOL(CLR)
            CALL PDFMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            CALL PDFDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 60         CONTINUE
*
            CALL PDFMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 70 I=1,NP1
            DO 61 J=4,ND
 61         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL PDFCOL(CLR)
            CALL PDFDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 70         CONTINUE
 80         CONTINUE
         ENDIF
      ENDIF
*
      CALL PDFW('S',0,IE)
      IF(LCR) CALL PDFCOL(OACC)
 90   CALL GENPLY(IA)
*
      RETURN
      END
*
      SUBROUTINE PDFCHA(STR,LST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STR*1024,STRE*1200
*
      IF(LCM.EQ.1) CALL PDFW('S',0,IE)
      IF(LDOTC.EQ.1) CALL PDFW('/F1 0.03 Tf',0,IE)
*
      CALL PROJEC(XC(1),XC(2),XC(3),X,Y)
      WRITE(STRE,'(A,2(F0.9,X),A)') 'BT ',X,Y,'Td'
      CALL PDFW(STRE,0,IE)
      CALL POSESC(STR,LST,STRE,LSTE)
      CALL PDFW('('//STRE(1:LSTE)//') Tj ET',0,IE)
*
      CALL GENCHA
*
      RETURN
      END
*
      SUBROUTINE PDFCOL(CLR)
*     **********************
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CLR(4)
      CHARACTER STR*100
*
      IF(LCM.EQ.1) CALL PDFW('S',0,IE)
*
      WRITE(STR,'(2(3(F0.6,X),A))') (CLR(J),J=1,3),'rg ',
     *                              (CLR(J),J=1,3),'RG'
      CALL PDFW(STR,0,IE)
*
      CALL GENCOL(CLR)
*
      RETURN
      END
*
      SUBROUTINE PDFWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STR*100
*
      IF(LCM.EQ.1) CALL PDFW('S',0,IE)
*
      WRITE(STR,'(F0.9,A)') 0.0015D0*MAX(0.5D0,AJ),' w'
      CALL PDFW(STR,0,IE)
*
      CALL GENWID(AJ)
*
      RETURN
      END
*
      SUBROUTINE PDFEND(NARI)
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
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----PDF DRIVER ------------------------------------------------------------
      INTEGER NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
      COMMON /PDFCOM/ NPOS,NOBJ,IXREF,ISTREAM,NSH,ISH
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
      CHARACTER STR*1000,DATE*8,TIME*10
*
      IF(NSH.NE.ISH) THEN
         PRINT*,'@@@ ERROR IN PDFEND, WRONG NUMBER OF GRPOLY SHADING'
         CALL FOXSTP(1)
      ENDIF
*
      IF(LCM.EQ.1) CALL PDFW('S',0,IE)
*
      IF(NARI.GT.0) THEN
         WRITE(STR,'(2(3(F0.6,X),A))') 0.,0.,0.,'rg ',0.,0.,0.,'RG'
         CALL PDFW(STR,0,IE)
         WRITE(STR,'(A,2(F0.9,X),A)') 'BT /F1 0.03 Tf ',0.007,0.01,'Td'
         CALL PDFW(STR,0,IE)
         WRITE(STR,'(A,I0,A)') '($$$ WARNING, GRPOLY ARITHMETIC '
     *           //'FAILURES: ',NARI,') Tj ET'
         CALL PDFW(STR,0,IE)
      ENDIF
*
      LSTREAM = NPOS-ISTREAM-1
      CALL PDFW('endstream endobj',0,IE)
*
*     THE LENGTH OF THE PDF MAIN GRAPHICS CONTENTS - THE (NSH+2)-TH PDF OBJECT
*
      WRITE(STR,'(I0,A)') LSTREAM,' endobj'
      CALL PDFW('obj',1,IE)
      CALL PDFW(STR,0,IE)
*
*     THE PDF PAGE OBJECT OF THE COSY PDF GRAPHICS - THE (NSH+3)-TH PDF OBJECT
*
      CALL PDFW('obj',1,IE)
      WRITE(STR,'(A,I0,A,I0,A)')
     *'<</Type /Page /Parent ',NOBJ+1,' 0 R /Contents ',NSH+1,' 0 R'
      CALL PDFW(STR,0,IE)
      CALL PDFW('/MediaBox [0 0 792 612] /CropBox [21 50 771 550]',0,IE)
      CALL PDFW('/VP [<</Type /Viewport /BBox [21 50 771 550] '
     *        //'/Name (COSY coordinate system)',0,IE)
      CALL PDFW('/Measure <</Type /Measure /Subtype /RL /R '
     *        //'(COSY units)',0,IE)
      WRITE(STR,'(A,F0.10,A)') '/X [<</U () /C ',
     *          (BD2F(1,2)-BD2F(1,1))/750.D0,' /D 100000>>]'
      CALL PDFW(STR,0,IE)
      WRITE(STR,'(A,F0.10,A)') '/Y [<</U () /C ',
     *          (BD2F(2,2)-BD2F(2,1))/500.D0,' /D 100000>>]'
      CALL PDFW(STR,0,IE)
      WRITE(STR,'(A,2(F0.10,A))') '/O [',
     *          21.D0-750.D0/(BD2F(1,2)-BD2F(1,1))*BD2F(1,1),' ',
     *          50.D0-500.D0/(BD2F(2,2)-BD2F(2,1))*BD2F(2,1),']'
      CALL PDFW(STR,0,IE)
      CALL PDFW('/D [<</U () /C 1 /D 100000>>] '
     *        //'/A [<</U () /C 1 /D 100000>>] /CYX 1>> >>]',0,IE)
      CALL PDFW('/Resources <<',0,IE)
      CALL PDFW('/Font <</F1 <</Type /Font /Subtype /Type1 /BaseFont '
     *        //'/Helvetica>> >>',0,IE)
      CALL PDFW('/Shading <<',0,IE)
      DO 10 I=1,NSH
      WRITE(STR,'(A,I0,A,I0,A)') '/Sh',I,' ',I,' 0 R'
      CALL PDFW(STR,0,IE)
 10   CONTINUE
      CALL PDFW('>> >> >> endobj',0,IE)
*
*     THE PDF PAGE TREE - THE (NSH+4)-TH PDF OBJECT
*
      CALL PDFW('obj',1,IE)
      WRITE(STR,'(A,I0,A)') '<</Type /Pages /Kids [',NOBJ-1,
     *          ' 0 R] /Count 1>> endobj'
      CALL PDFW(STR,0,IE)
*
*     THE PDF ROOT, CATALOG DICTIONARY - THE (NSH+5)-TH PDF OBJECT
*
      CALL PDFW('obj',1,IE)
      WRITE(STR,'(A,I0,A)') '<</Type /Catalog /Pages ',NOBJ-1,
     *          ' 0 R>> endobj'
      CALL PDFW(STR,0,IE)
*
*     THE PDF INFORMATION DICTIONARY - THE (NSH+6)-TH PDF OBJECT
*
      CALL DATE_AND_TIME(DATE,TIME)
      CALL PDFW('obj',1,IE)
      CALL PDFW('<</Title (COSY INFINITY PDF Graphics)',0,IE)
      CALL PDFW('/CreationDate (D:'//DATE//TIME(1:6)//')',0,IE)
      CALL PDFW('/Producer (COSY INFINITY Graphics)',0,IE)
      CALL PDFW('/Creator (COSY INFINITY) >> endobj',0,IE)
*
*     COSY COMMENTS
*
      CALL ENDCMT(NARI,-12)
*
*     THE PDF CROSS-REFERENCE TABLE AND THE TRAILER
*
      CALL PDFW('xref',1,IE)
      WRITE(STR,'(A,I0)') '0 ',NOBJ
      CALL PDFW(STR,0,IE)
      CALL PDFW('0000000000 65535 f'//SCCR,0,IE)
      DO 20 I=1,NOBJ-1
      WRITE(STR,'(I10.10,A)') NC(NBEG(IXREF)+I-1),' 00000 n'//SCCR
      CALL PDFW(STR,0,IE)
 20   CONTINUE
      WRITE(STR,'(A,I0,A,I0,A,I0,A)') 'trailer <</Size ',NOBJ,
     *          ' /Root ',NOBJ-2,' 0 R /Info ',NOBJ-1,' 0 R>>'
      CALL PDFW(STR,0,IE)
      CALL PDFW('startxref',0,IE)
      WRITE(STR,'(I0)') NC(NBEG(IXREF)+NOBJ-1)
      CALL PDFW(STR,0,IE)
      CALL PDFW('%%EOF',0,IE)
*
      CALL FOXDAL(IXREF,1)
*
      CLOSE(4)
*
      CALL GENEND
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     SVG DRIVER ROUTINES                                                     *
*     CONTRIBUTED BY ALEXANDER WITTIG (2012)                                  *
*     MODIFIED BY K. MAKINO (2013)                                            *
*     http://www.w3.org/TR/SVG/                                               *
*                                                                             *
*******************************************************************************
*     EDIT DESCRIPTORS Iw, Fw.d: IF w=0, THE PROCESSOR SELECTS THE FIELD WIDTH.
*     EDIT DESCRIPTOR $ SUPPRESSES TRAILING CARRIAGE RETURN. (INTEL)
*
      SUBROUTINE SVGSTL(STL,IL,IW,IWR)
*     ********************************
*
*     THIS SUBROUTINE PREPARES A CSS STYLE STRING (STL) OF LENGTH IL
*     BASED ON THE CURRENT GRAPH CONTEXT.
*     IW=1:  INCLUDE THE CURRENT LINE WIDTH STYLE
*     IWR=1: INCLUDE THE WIRE FRAME STYLE FOR SVGTRI
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----SVG DRIVER ------------------------------------------------------------
      LOGICAL LPLSVG
      COMMON /SVGCOM/ LPLSVG
*----------------------------------------------------------------------------
*
      CHARACTER STL*200
      LOGICAL LSTYLE
*
      STL = ''
      IL = 0
*
*     BYPASS WHEN SVGSTL IS CALLED BY SVGDRA AND SVGTRI THAT ARE CALLED
*     WITHIN SVGPLY FOR SOLID COLOR POLYNOMIALS.
*
      IF(LPLSVG) RETURN
*
      LSTYLE = .FALSE.
      LSTYLE = ((ACC(1)+ACC(2)+ACC(3)+1.D0-ACC(4)).GT.0.D0)
     *         .OR.(IWR.EQ.1).OR.(IW.EQ.1.AND.ACW.NE.1.D0)
      IF(.NOT.LSTYLE) RETURN
*
      IL = 8
      STL(1:IL) = ' style="'
*
      IF((ACC(1)+ACC(2)+ACC(3)).GT.0.D0) THEN
         WRITE(STL(IL+1:),'(A,3Z2.2,A)') 'color:#',
     *        (INT(ACC(J)*255.D0),J=1,3),';'
         IL = ILAST(STL,1,200)
      ENDIF
*
      IF(ACC(4).LT.1.D0) THEN
         WRITE(STL(IL+1:),'(A,F0.3,A)') 'opacity:',ACC(4),';'
         IL = ILAST(STL,1,200)
      ENDIF
*
      IF(IW.EQ.1.AND.ACW.NE.1.D0) THEN
         WRITE(STL(IL+1:),'(A,F0.1,A)') 'stroke-width:',2.D0*ACW,';'
         IL = ILAST(STL,1,200)
      ENDIF
*
      IF(IWR.EQ.1) THEN
         WRITE(STL(IL+1:),'(A)') 'fill: none;'
         IL = ILAST(STL,1,200)
      ENDIF
*
      STL(IL+1:IL+1) = '"'
      IL = IL+1
*
      RETURN
      END
*
      SUBROUTINE SVGESC(STR,LST,STRE,LSTE)
*     ************************************
*
*     THIS SUBROUTINE ESCAPES A STRING TO BE INCLUDED IN SVG OR OTHER XML
*     FROM STR (LENGTH LST) TO STRE (LENGTH LSTE)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CHARACTER STR*(*),STRE*(*)
*
      J = 0
      DO 10 I=1,LST
      IF(STR(I:I).EQ.'<') THEN
         STRE(J+1:J+4) = '&lt;'
         J = J+4
      ELSEIF(STR(I:I).EQ.'&') THEN
         STRE(J+1:J+5) = '&amp;'
         J = J+5
      ELSE
         J = J+1
         STRE(J:J) = STR(I:I)
      ENDIF
      IF(J+5.GT.LEN(STRE)) GOTO 20
 10   CONTINUE
 20   LSTE = J
*
      RETURN
      END
*
      SUBROUTINE SVGPOS(X,Y,XS,YS)
*     ****************************
*
*     THIS SUBROUTINE CONVERTS THE COSY 2D (X,Y) COORDINATES IN [0,1]^2
*     TO THE SVG viewBox COORDINATES (XS,YS).
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER(SVGW=1500.D0,SVGH=1000.D0)
*
      XS = SVGW*X
      YS = SVGH*(1.D0-Y)
*
      RETURN
      END
*
      SUBROUTINE SVGBEG
*     *****************
*
*-----GRAPH FILE OUTPUT -----------------------------------------------------
      INTEGER ICNT
      CHARACTER FBASE*20
      COMMON /GRFILE/ ICNT,FBASE
*----------------------------------------------------------------------------
*
      PARAMETER(SVGW=1500.D0,SVGH=1000.D0)
      CHARACTER NUMBER*3
*
      WRITE(NUMBER,'(I3.3)') MOD(ICNT,1000)
      IL = ILAST(FBASE,1,20)
      OPEN(4,FILE=FBASE(1:IL)//NUMBER//'.svg',STATUS='UNKNOWN')
      ICNT = ICNT+1
*
      WRITE(4,'(A)') '<?xml version="1.0" encoding="UTF-8"?>'
      WRITE(4,'(A)') '<svg xmlns="http://www.w3.org/2000/svg" '//
     *               'version="1.1" style="color:black;"'
      WRITE(4,'(A,2I5,A)') 'width="9in" height="6in" '//
     *               'viewBox="0 0',NINT(SVGW),NINT(SVGH),
     *               '" overflow="hidden">'
      WRITE(4,'(A)') '<!--'
      WRITE(4,'(A)') 'To change the size of the picture, change the '//
     *               'width and height in the line above in this file.'
      WRITE(4,'(A)') '-->'
      WRITE(4,'(A)') '<style type="text/css">'
      WRITE(4,'(A)') 'polyline { stroke: currentcolor; fill: none; '//
     *               'stroke-width: 2; opacity: 1;'
      WRITE(4,'(A)') '           stroke-linecap: round; '//
     *               'stroke-linejoin: round; }'
      WRITE(4,'(A)') 'polygon  { fill: currentcolor; opacity: 1; '//
     *               'stroke: currentcolor; stroke-width: 2; }'
      WRITE(4,'(A)') 'circle   { fill: currentcolor; opacity: 1; '//
     *               'stroke: currentcolor; stroke-width: 0; }'
      WRITE(4,'(A)') 'text     { fill: currentcolor; opacity: 1; '//
     *               'font-family: Helvetica, Arial, sans-serif; '//
     *               'font-size: 30px; }'
      WRITE(4,'(A)') '</style>'
      WRITE(4,'(A)') '<title>COSY INFINITY SVG Graphic (pic'//NUMBER//
     *               '.svg)</title>'
*
      CALL GENBEG
*
      RETURN
      END
*
      SUBROUTINE SVGMOV(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      IF(LCM.EQ.1) WRITE(4,'(A)') '"/>'
*
      CALL GENMOV(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE SVGDRA(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STL*200
*
      IF(LCM.NE.1) THEN
         CALL PROJEC(XC(1),XC(2),XC(3),X,Y)
         CALL SVGPOS(X,Y,XS,YS)
         CALL SVGSTL(STL,IL,1,0)
         WRITE(4,'(2(A,F0.1),$)') '<polyline'//STL(1:IL)//' points="',
     *           XS,' ',YS
      ENDIF
      CALL PROJEC(XXX,YYY,ZZZ,X,Y)
      CALL SVGPOS(X,Y,XS,YS)
      WRITE(4,'(2(A,F0.1),$)') ', ',XS,' ',YS
*
      CALL GENDRA(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE SVGDOT(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STL*200
      LOGICAL LGRDOT
*
      IF(LCM.EQ.1) WRITE(4,'(A)') '"/>'
*
      IF(LGRDOT(XXX,YYY,ZZZ)) THEN
         CALL PROJEC(XXX,YYY,ZZZ,X,Y)
         CALL SVGPOS(X,Y,XS,YS)
         CALL SVGSTL(STL,IL,0,0)
         WRITE(4,'(3(A,F0.1),A)') '<circle cx="',XS,'" cy="',YS,
     *           '" r="',ACW,'"'//STL(1:IL)//'/>'
*
*        CAUTION:
*        IT TAKES VERY LONG TO DISPLAY A HUGE SIZED COSY TRACKING PICTURE,
*        OR DISPLAYING SOFTWARE (WEB BROWSERS, INKSCAPE) MAY BE HUNG.
*        USING A CHARACTER SUCH AS A PERIOD INSTEAD OF A CIRCLE ABOVE
*        DOES NOT HELP, OR MAY BE EVEN WORSE.
*
      ENDIF
*
      CALL GENDOT(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE SVGTRI(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STL*200
      LOGICAL LGRTRI
*
      IF(LCM.EQ.1) WRITE(4,'(A)') '"/>'
*
      IF(LGRTRI(XXX,YYY,ZZZ,XP(1),XP(2),XP(3),XC(1),XC(2),XC(3))) THEN
         CALL PROJEC(XC(1),XC(2),XC(3),CX,CY)
         CALL SVGPOS(CX,CY,CXS,CYS)
         CALL PROJEC(XP(1),XP(2),XP(3),PX,PY)
         CALL SVGPOS(PX,PY,PXS,PYS)
         CALL PROJEC(XXX,YYY,ZZZ,X,Y)
         CALL SVGPOS(X,Y,XS,YS)
         CALL SVGSTL(STL,IL,0,LWR)
         WRITE(4,'(A,6(F0.1,A))') '<polygon'//STL(1:IL)//' points="',
     *           CXS,' ',CYS,',',PXS,' ',PYS,',',XS,' ',YS,'"/>'
      ENDIF
*
      CALL GENTRI(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE SVGPLY(IA,IPST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*-----SVG DRIVER ------------------------------------------------------------
      LOGICAL LPLSVG
      COMMON /SVGCOM/ LPLSVG
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION OACC(4),CLR(4),PC(7)
      CHARACTER STL*200
      LOGICAL LGRPLY,LCR
*
      IF(LCM.EQ.1) WRITE(4,'(A)') '"/>'
      LPLSVG = .FALSE.
*
      IPLYC = IPLY+1
      IF(IPST.EQ.-1) THEN
         WRITE(4,'(A)') '<!--'
         WRITE(4,'(A,I0,5X,A)') 'POLY# ',IPLYC,'arithmetic failure'
         WRITE(4,'(A)') '-->'
         GOTO 90
      ENDIF
      IF(.NOT.LGRPLY(IA)) THEN
         WRITE(4,'(A)') '<!--'
         WRITE(4,'(A,I0,5X,A)') 'POLY# ',IPLYC,'out of bound'
         WRITE(4,'(A)') '-->'
         GOTO 90
      ENDIF
*
      CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
      CALL PLYCMT(IA,IOBJ,NP1,NP2,LCR,-13)
*
      IF(LCR) THEN
         ND = 7
         DO 10 K=1,4
         OACC(K) = ACC(K)
 10      CONTINUE
      ELSE
         ND = 3
         CALL SVGSTL(STL,IL,0,0)
         WRITE(4,'(A)') '<g'//STL(1:IL)//'>'
         LPLSVG = .TRUE.
      ENDIF
*
*     CURVE
*     *****
*
      IF(IOBJ.EQ.1) THEN
         CALL PLYEV(IA,ND,0.D0,NP1)
         CALL SVGMOV(PNT(1,0),PNT(2,0),PNT(3,0))
         DO 20 I=1,NP1
         DO 11 J=4,ND
 11      CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
         IF(LCR) CALL SVGCOL(CLR)
         CALL SVGDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 20      CONTINUE
*
*     SURFACE
*     *******
*
      ELSEIF(IOBJ.EQ.2) THEN
*
*      * VARYING OR SOLID COLOR TRIANGLE FILL PAINTING *
         IF(LWR.EQ.0) THEN
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            DO 40 I2=1,NP2
            DO 21 I=0,NP1
            DO 21 J=1,ND
 21         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 40 I=1,NP1
            DO 22 J=1,ND
 22         PC(J) = (PNTO(J,I-1)+PNTO(J,I)+PNT(J,I-1)+PNT(J,I))*0.25D0
*
            CALL SVGMOV(PNT(1,I),PNT(2,I),PNT(3,I))
            CALL SVGMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            DO 31 J=4,ND
 31         CLR(J-3) = (PNT(J,I)+PNTO(J,I)+PC(J))/3.D0
            IF(LCR) CALL SVGCOL(CLR)
            CALL SVGTRI(PC(1),PC(2),PC(3))
*
            DO 32 J=4,ND
 32         CLR(J-3) = (PNTO(J,I)+PC(J)+PNTO(J,I-1))/3.D0
            IF(LCR) CALL SVGCOL(CLR)
            CALL SVGTRI(PNTO(1,I-1),PNTO(2,I-1),PNTO(3,I-1))
*
            DO 33 J=4,ND
 33         CLR(J-3) = (PC(J)+PNTO(J,I-1)+PNT(J,I-1))/3.D0
            IF(LCR) CALL SVGCOL(CLR)
            CALL SVGTRI(PNT(1,I-1),PNT(2,I-1),PNT(3,I-1))
*
            CALL SVGMOV(PC(1),PC(2),PC(3))
            DO 34 J=4,ND
 34         CLR(J-3) = (PNT(J,I-1)+PC(J)+PNT(J,I))/3.D0
            IF(LCR) CALL SVGCOL(CLR)
            CALL SVGTRI(PNT(1,I),PNT(2,I),PNT(3,I))
 40         CONTINUE
*
*      * WIRE FRAME QUADRILATERALS *
         ELSE
            CALL PLYEV(IA,ND,-1.D0,NP1)
*
            CALL SVGMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 50 I=1,NP1
            DO 41 J=4,ND
 41         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL SVGCOL(CLR)
            CALL SVGDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 50         CONTINUE
*
            DO 80 I2=1,NP2
            DO 51 I=0,NP1
            DO 51 J=1,ND
 51         PNTO(J,I) = PNT(J,I)
            CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
            DO 60 I=0,NP1
            DO 52 J=4,ND
 52         CLR(J-3) = (PNTO(J,I)+PNT(J,I))*0.5D0
            IF(LCR) CALL SVGCOL(CLR)
            CALL SVGMOV(PNTO(1,I),PNTO(2,I),PNTO(3,I))
            CALL SVGDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 60         CONTINUE
*
            CALL SVGMOV(PNT(1,0),PNT(2,0),PNT(3,0))
            DO 70 I=1,NP1
            DO 61 J=4,ND
 61         CLR(J-3) = (PNT(J,I-1)+PNT(J,I))*0.5D0
            IF(LCR) CALL SVGCOL(CLR)
            CALL SVGDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 70         CONTINUE
 80         CONTINUE
         ENDIF
      ENDIF
*
      IF(LCM.EQ.1) WRITE(4,'(A)') '"/>'
      LCM = 0
      IF(LCR) THEN
         CALL SVGCOL(OACC)
      ELSE
         WRITE(4,'(A)') '</g>'
      ENDIF
*
 90   CALL GENPLY(IA)
      LPLSVG = .FALSE.
*
      RETURN
      END
*
      SUBROUTINE SVGCHA(STR,LST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STR*1024,STRE*1024,STL*200
*
      IF(LCM.EQ.1) WRITE(4,'(A)') '"/>'
*
      CALL SVGESC(STR,LST,STRE,LSTE)
      CALL PROJEC(XC(1),XC(2),XC(3),X,Y)
      CALL SVGPOS(X,Y,XS,YS)
      CALL SVGSTL(STL,IL,0,0)
      WRITE(4,'(2(A,F0.1),A)') '<text x="',XS,'" y="',YS,
     *        '"'//STL(1:IL)//'>'//STRE(1:LSTE)//'</text>'
*
      CALL GENCHA
*
      RETURN
      END
*
      SUBROUTINE SVGCOL(CLR)
*     **********************
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CLR(4)
*
      IF(LCM.EQ.1) WRITE(4,'(A)') '"/>'
*
      CALL GENCOL(CLR)
*
      RETURN
      END
*
      SUBROUTINE SVGWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      IF(LCM.EQ.1) WRITE(4,'(A)') '"/>'
*
      CALL GENWID(AJ)
*
      RETURN
      END
*
      SUBROUTINE SVGEND(NARI)
*     ***********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CHARACTER STR*100
*
      IF(LCM.EQ.1) WRITE(4,'(A)') '"/>'
*
      IF(NARI.GT.0) THEN
         X = 0.007D0
         Y = 0.010D0
         CALL SVGPOS(X,Y,XS,YS)
         WRITE(STR,'(2(A,F0.1),A,I0,A)') '<text x="',XS,'" y="',YS,
     *      '">$$$ WARNING, GRPOLY ARITHMETIC FAILURES: ',NARI,'</text>'
         WRITE(4,'(A)') STR(1:ILAST(STR,1,100))
      ENDIF
*
      CALL ENDCMT(NARI,-13)
      WRITE(4,'(A)') '</svg>'
      CLOSE(4)
*
      CALL GENEND
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     STL DRIVER ROUTINES                                                     *
*     CONTRIBUTED BY ALEXANDER WITTIG (2012)                                  *
*     MODIFIED BY K. MAKINO (2013)                                            *
*     https://secure.wikimedia.org/wikipedia/en/wiki/STL_%28file_format%29    *
*                                                                             *
*******************************************************************************
*
      SUBROUTINE STLBEG(IUNIT)
*     ************************
*
*-----STL DRIVER ------------------------------------------------------------
      INTEGER*4 NTRI,IMODE
      COMMON /STLCOM/ NTRI,IMODE
*----------------------------------------------------------------------------
*-----GRAPH FILE OUTPUT -----------------------------------------------------
      INTEGER ICNT
      CHARACTER FBASE*20
      COMMON /GRFILE/ ICNT,FBASE
*----------------------------------------------------------------------------
*
      CHARACTER NUMBER*3
*
      WRITE(NUMBER,'(I3.3)') MOD(ICNT,1000)
      IL = ILAST(FBASE,1,20)
      NTRI = 0
      IF(IUNIT.EQ.-14) THEN
         IMODE = 0
         OPEN(4,FILE=FBASE(1:IL)//NUMBER//'.stl',STATUS='UNKNOWN')
         WRITE(4,'(A)') 'solid COSY GENERATED ASCII STL FILE'
      ELSE
         IMODE = 1
C        SOME FORTRAN COMPILERS SUCH AS GNU CANNOT HANDLE ACCESS='STREAM'.
         OPEN(4,FILE=FBASE(1:IL)//NUMBER//'.stl',STATUS='UNKNOWN',       *STL
     *          ACCESS='STREAM')                                         *STL
         WRITE(4) 'COSY GENERATED BINARY STL FILE                    '//
     *            '                              ',NTRI
      ENDIF
      ICNT = ICNT+1
*
      CALL GENBEG
*
      RETURN
      END
*
      SUBROUTINE STLMOV(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL GENMOV(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE STLDRA(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL GENDRA(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE STLDOT(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL GENDOT(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE STLTRI(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*-----STL DRIVER ------------------------------------------------------------
      INTEGER*4 NTRI,IMODE
      COMMON /STLCOM/ NTRI,IMODE
*----------------------------------------------------------------------------
*
      INTEGER*2 NACT
      DOUBLE PRECISION AFL
      DATA AFL / 1.D0 /
      SAVE AFL
*
      NTRI = NTRI+1
      IF(LCM.NE.3) AFL=1.D0
*
      V1 = XP(1)-XC(1)
      V2 = XP(2)-XC(2)
      V3 = XP(3)-XC(3)
      W1 = XXX-XP(1)
      W2 = YYY-XP(2)
      W3 = ZZZ-XP(3)
      A1 = AFL*(V2*W3-V3*W2)
      A2 = AFL*(V3*W1-V1*W3)
      A3 = AFL*(V1*W2-V2*W1)
      AN = SQRT(A1*A1+A2*A2+A3*A3)
*
      AFL = -AFL
*
      IF(IMODE.EQ.0) THEN
         WRITE(4,'(3X,A,3(X,E16.8))') 'facet normal',A1/AN,A2/AN,A3/AN
         WRITE(4,'(6X,A)') 'outer loop'
         IF(AFL.LT.0.D0) THEN
            WRITE(4,'(9X,A,3(X,E16.8))') 'vertex',XP(1),XP(2),XP(3)
            WRITE(4,'(9X,A,3(X,E16.8))') 'vertex',XC(1),XC(2),XC(3)
         ELSE
            WRITE(4,'(9X,A,3(X,E16.8))') 'vertex',XC(1),XC(2),XC(3)
            WRITE(4,'(9X,A,3(X,E16.8))') 'vertex',XP(1),XP(2),XP(3)
         ENDIF
         WRITE(4,'(9X,A,3(X,E16.8))') 'vertex',XXX,YYY,ZZZ
         WRITE(4,'(6X,A)') 'endloop'
         WRITE(4,'(3X,A)') 'endfacet'
      ELSE
         WRITE(4) REAL(A1/AN),REAL(A2/AN),REAL(A3/AN)
         IF(AFL.LT.0.D0) THEN
            WRITE(4) REAL(XP(1)),REAL(XP(2)),REAL(XP(3))
            WRITE(4) REAL(XC(1)),REAL(XC(2)),REAL(XC(3))
         ELSE
            WRITE(4) REAL(XC(1)),REAL(XC(2)),REAL(XC(3))
            WRITE(4) REAL(XP(1)),REAL(XP(2)),REAL(XP(3))
         ENDIF
         WRITE(4) REAL(XXX),REAL(YYY),REAL(ZZZ)
         NACT = 0
         WRITE(4) NACT
      ENDIF
*
      CALL GENTRI(XXX,YYY,ZZZ)
*
      RETURN
      END
*
      SUBROUTINE STLPLY(IA,IPST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GRPOLY DISCRETIZED POINTS ---------------------------------------------
      PARAMETER(NGP=999)
      DOUBLE PRECISION PNT(7,0:NGP),PNTO(7,0:NGP)
      COMMON /GRPNT/ PNT,PNTO
*----------------------------------------------------------------------------
*
      LOGICAL LCR
*
      IF(IPST.EQ.-1) GOTO 90
*
      CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
      ND = 3
*
*     CURVE
*     *****
*
      IF(IOBJ.EQ.1) THEN
         CALL PLYEV(IA,ND,0.D0,NP1)
         CALL STLMOV(PNT(1,0),PNT(2,0),PNT(3,0))
         DO 10 I=1,NP1
         CALL STLDRA(PNT(1,I),PNT(2,I),PNT(3,I))
 10      CONTINUE
*
*     SURFACE
*     *******
*
      ELSEIF(IOBJ.EQ.2) THEN
         CALL PLYEV(IA,ND,-1.D0,NP1)
*
         DO 30 I2=1,NP2
         DO 20 I=0,NP1
         DO 20 J=1,ND
 20      PNTO(J,I) = PNT(J,I)
         CALL PLYEV(IA,ND,(2.D0*I2)/NP2-1.D0,NP1)
*
         CALL STLMOV(PNTO(1,0),PNTO(2,0),PNTO(3,0))
         CALL STLMOV(PNT(1,0),PNT(2,0),PNT(3,0))
         DO 30 I=1,NP1
         CALL STLTRI(PNTO(1,I),PNTO(2,I),PNTO(3,I))
         CALL STLTRI(PNT(1,I),PNT(2,I),PNT(3,I))
 30      CONTINUE
      ENDIF
*
 90   CALL GENPLY(IA)
*
      RETURN
      END
*
      SUBROUTINE STLCHA(STR,LST)
*     **************************
*
      CHARACTER STR*1024
*
      CALL GENCHA
*
      RETURN
      END
*
      SUBROUTINE STLCOL(CLR)
*     **********************
*
      DOUBLE PRECISION CLR(4)
*
      CALL GENCOL(CLR)
*
      RETURN
      END
*
      SUBROUTINE STLWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL GENWID(AJ)
*
      RETURN
      END
*
      SUBROUTINE STLEND(NARI)
*     ***********************
*
*-----STL DRIVER ------------------------------------------------------------
      INTEGER*4 NTRI,IMODE
      COMMON /STLCOM/ NTRI,IMODE
*----------------------------------------------------------------------------
*-----GRAPH FILE OUTPUT -----------------------------------------------------
      INTEGER ICNT
      CHARACTER FBASE*20
      COMMON /GRFILE/ ICNT,FBASE
*----------------------------------------------------------------------------
*
      IF(IMODE.EQ.0) THEN
         WRITE(4,'(A)') 'endsolid COSY GENERATED ASCII STL FILE'
      ELSE
C        SOME FORTRAN COMPILERS SUCH AS GNU CANNOT HANDLE "POS=" IDENTIFIER.
         WRITE(4,POS=81) NTRI                                            *STL
      ENDIF
      CLOSE(4)
*
      IF(NARI.GT.0) THEN
         WRITE(6,'(A,I0,A,I3.3,A)') ' $$$ WARNING, GRPOLY ARITHMETIC '
     *      //'FAILURES: ',NARI,'  IN STL OUTPUT FILE '//FBASE(1:IL),
     *      MOD(ICNT-1,1000),'.stl'
      ENDIF
*
      CALL GENEND
*
      RETURN
      END
*
*******************************************************************************
*                                                                             *
*     GUI DRIVER ROUTINES                                                     *
*     CONTRIBUTED BY ALEXANDER WITTIG (2012)                                  *
*                                                                             *
*******************************************************************************
*     EDIT DESCRIPTORS Iw, Fw.d: IF w=0, THE PROCESSOR SELECTS THE FIELD WIDTH.
*     EDIT DESCRIPTOR $ SUPPRESSES TRAILING CARRIAGE RETURN. (INTEL)
*
      SUBROUTINE GUIBEG(IUNIT,NELE)
*     *****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*-----GRAPHICS BOUNDS AND PROJECTION PARAMETERS -----------------------------
      DOUBLE PRECISION BD2F(2,2),BD3(3,2),BDG(3,2),BDG2(2,2),
     *       PRANG(3),PHI,THETA,SP,CP,ST,CT
      INTEGER IZOOM
      COMMON /GRCOM/ BD2F,BD3,BDG,BDG2,PRANG,PHI,THETA,SP,CP,ST,CT,IZOOM
*----------------------------------------------------------------------------
*
      IF(LGUI.EQ.1) THEN
*        OUTPUT NUMBERS TO GUI USING HEXADECIMAL VALUES, Z8.8 FOR REAL
         WRITE(6,'(A,I0)')    'E',NELE
         WRITE(6,'(A,6Z8.8)') 'B',REAL(BD3(1,1)),REAL(BD3(1,2)),
     *                            REAL(BD3(2,1)),REAL(BD3(2,2)),
     *                            REAL(BD3(3,1)),REAL(BD3(3,2))
         WRITE(6,'(A,2Z8.8)') 'V',REAL(PHI),REAL(THETA)
         CALL GENBEG
      ELSE
         CALL TTYBEG
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GUIMOV(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      IF(LGUI.EQ.1) THEN
         WRITE(6,'(A,3Z8.8)') 'M',REAL(XXX),REAL(YYY),REAL(ZZZ)
         CALL GENMOV(XXX,YYY,ZZZ)
      ELSE
         CALL TTYMOV(XXX,YYY,ZZZ)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GUIDRA(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      IF(LGUI.EQ.1) THEN
         WRITE(6,'(A,3Z8.8)') 'D',REAL(XXX),REAL(YYY),REAL(ZZZ)
         CALL GENDRA(XXX,YYY,ZZZ)
      ELSE
         CALL TTYDRA(XXX,YYY,ZZZ)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GUIDOT(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      IF(LGUI.EQ.1) THEN
         WRITE(6,'(A,3Z8.8)') 'X',REAL(XXX),REAL(YYY),REAL(ZZZ)
         CALL GENDOT(XXX,YYY,ZZZ)
      ELSE
         CALL TTYDOT(XXX,YYY,ZZZ)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GUITRI(XXX,YYY,ZZZ)
*     ******************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      IF(LGUI.EQ.1) THEN
         WRITE(6,'(A,3Z8.8)') 'T',REAL(XXX),REAL(YYY),REAL(ZZZ)
         CALL GENTRI(XXX,YYY,ZZZ)
      ELSE
         CALL TTYTRI(XXX,YYY,ZZZ)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GUIPLY(IA,IPST)
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
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      INTEGER IPLEN(7),IPRER(7)
      LOGICAL LCR
*
      IF(LGUI.EQ.1) THEN
         IF(IPST.NE.-1) THEN
*           OUTPUT POLYNOMIAL TO GUI (HEX RANGE: Z2.2 = 0,255, Z3.3 = 0,4095)
            CALL PLYPRE(IA,IOBJ,NP1,NP2,LCR,0)
            WRITE(6,'(A,2Z3.3,$)') 'P',NP1,NP2
            CALL PLYLEN(IA,IPLEN,IPRER)
            IPOS = IA+2
            DO 20 I=1,7
            WRITE(6,'(Z3.3,$)') IPLEN(I)
            DO 10 J=IPOS,IPOS+IPLEN(I)-IPRER(I)-1
            IP2 = INT(NC(J)/100)
            IP1 = NC(J)-IP2*100
            WRITE(6,'(2Z2.2,Z8.8,$)') IP1,IP2,REAL(CC(J))
 10         CONTINUE
*
*           TM REMAINDER BOUND
*
            IF(IPRER(I).EQ.1) THEN
               WRITE(6,'(A,Z8.8,$)') 'R000',REAL(CC(J))
            ENDIF
            IPOS = IPOS+IPLEN(I)
 20         CONTINUE
            WRITE(6,'(A)') ''
         ELSE
*
*           TM ARITHMETIC FAILURE
*
            WRITE(6,'(A)') 'A'
         ENDIF
      ELSE
         CALL TTYPLY(IA,IPST)
      ENDIF
*
      CALL GENPLY(IA)
*
      RETURN
      END
*
      SUBROUTINE GUICHA(STR,LST)
*     **************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      CHARACTER STR*1024
*
      IF(LGUI.EQ.1) THEN
         WRITE(6,'(A)') 'S'//STR(1:LST)
         CALL GENCHA
      ELSE
         CALL TTYCHA(STR,LST)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GUICOL(CLR)
*     **********************
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      DOUBLE PRECISION CLR(4)
      INTEGER ICLR(4)
*
*     IAREVS=1 TURNS THE OPACITY OPPOSITE, I.E. RGBA'S A=0 FOR FULL OPACITY.
*
      PARAMETER(IAREVS=0)
*
      IF(LGUI.EQ.1) THEN
         DO 10 J=1,4
         ICLR(J) = INT(CLR(J)*255.D0)
 10      CONTINUE
         IF(IAREVS.EQ.1) ICLR(4) = 255-ICLR(4)
         WRITE(6,'(A,4Z2.2)') 'R',(ICLR(J),J=1,4)
         CALL GENCOL(CLR)
      ELSE
         CALL TTYCOL(CLR)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GUIWID(AJ)
*     *********************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      IF(LGUI.EQ.1) THEN
         WRITE(6,'(A,I0)') 'W',NINT(AJ)
         CALL GENWID(AJ)
      ELSE
         CALL TTYWID(AJ)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GUISTL(ID,NLEN,I)
*     ****************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*-----GRAPH CONTEXT ---------------------------------------------------------
      DOUBLE PRECISION XC(3),XP(3),ACC(4),ACW,EPSXYZ,EPSCOL
      INTEGER LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
      COMMON /GRAPHCOM/ XC,XP,ACC,ACW,EPSXYZ,EPSCOL,
     *       LCM,LWR,LDOTW,LDOTC,NPLY,IPLY
*----------------------------------------------------------------------------
*
      CALL GENSTL(ID,NLEN,I)
*
      IF(LGUI.EQ.1) THEN
*
*        WIRE FRAME DRAWING FOR TRI AND PLY SURFACE (1: WIRE FRAME, 0: FILL)
*
         IF(ID.EQ.9) THEN
            WRITE(6,'(A,I1)') 'Y',LWR
*
*        DISCRETIZATION ACCURACY EPSXYZ, EPSCOL FOR PLY
*
         ELSEIF(ID.EQ.10) THEN
            WRITE(6,'(A,2Z8.8)') 'Q',REAL(EPSXYZ),REAL(EPSCOL)
         ENDIF
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE GUIEND(IUNIT,NARI)
*     *****************************
*
*-----GUI -------------------------------------------------------------------
      PARAMETER(LGV=100,LGW=10)
      INTEGER LGUI,NGVAR(LGW,2),IGVAR(LGW,LGV,2)
      COMMON /GUICOM/ LGUI,NGVAR,IGVAR
*----------------------------------------------------------------------------
*
      IF(LGUI.EQ.1) THEN
         WRITE(6,'(A)') 'F'
         CALL GENEND
      ELSE
         CALL TTYEND(6,NARI)
      ENDIF
*
      RETURN
      END
*
