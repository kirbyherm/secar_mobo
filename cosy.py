#!/usr/bin/env python

import sys, math
import os, shutil
#import commands
import subprocess as commands
import re
from numpy.random import random as rng
from numpy import array, append

qNom = array([-0.39773, 0.217880+0.001472, 0.242643-0.0005+0.000729, -0.24501-0.002549, 0.1112810+0.00111, 0.181721-0.000093+0.00010-0.000096, -0.0301435+0.0001215] )
qNew = array([-0.319612, 0.198440, 0.221305,-0.179820, 0.090087, 0.194583,-0.040515])

# Function that runs cosy given field gradients and outputs resolution at FP3. 
# Output is written in file temp-results
def cosyrun(qs=qNom):
    input_len = len(qs)
    if (len(qNom)-input_len>0):
        for i in range(len(qNom)-input_len):
            qs = append(qs,qNom[i+input_len])
    [q1s, q2s, q3s, q4s, q5s, q6s, q7s] = qs
    rand = rng()
    cosyFilename = 'pygmoCosy'+str(rand)+'.fox'
    while os.path.exists(cosyFilename):
        rand = rng()
        cosyFilename = 'pygmoCosy'+str(rand)+'.fox'
    lisFilename = 'pygmoCosy'+str(rand)+'.lis'
    # creating input file
    f = open(cosyFilename,'w')
    
    f.write('INCLUDE \'COSY\';\n')
    f.write('PROCEDURE RUN ;\n')
    f.write('VARIABLE B1N 1 7;\n')
    f.write('VARIABLE B1S1 1 7; VARIABLE B1S2 1 7;\n')
    f.write('VARIABLE B2S1 1 7; VARIABLE B2S2 1 7;\n')
    f.write('VARIABLE B3S1 1 7; VARIABLE B3S2 1 7;\n')
    f.write('VARIABLE B4S1 1 7; VARIABLE B4S2 1 7;\n')
    f.write('VARIABLE B5S1 1 7; VARIABLE B5S2 1 7;\n')
    f.write('VARIABLE B6S1 1 7; VARIABLE B6S2 1 7;\n')
    f.write('VARIABLE B7S1 1 7; VARIABLE B7S2 1 7;\n')
    f.write('VARIABLE B8S1 1 7; VARIABLE B8S2 1 7;\n')
    f.write('\n')
    f.write('{Wien filter higher order components}\n')
    f.write('VARIABLE NE1 1 4;\n')
    f.write('VARIABLE NM1 1 4;\n')
    f.write('\n')
    f.write('VARIABLE WV 1;\n')
    f.write('VARIABLE WW 1 2;\n')
    f.write('VARIABLE CENTER 1;\n')
    f.write('VARIABLE NMAX 1;\n')
    f.write('\n')
    f.write('VARIABLE XX 1;{X-SIZE}\n')
    f.write('VARIABLE AX 1;\n')
    f.write('VARIABLE YY 1;{Y-SIZE}\n')
    f.write('VARIABLE AY 1;\n')
    f.write('VARIABLE DE 1;\n')
    f.write('\n')
    f.write('{Ray definitions}\n')
    f.write('VARIABLE SRXX 1; \n')
    f.write('VARIABLE SRAX 1;\n')
    f.write('VARIABLE SRYY 1;\n')
    f.write('VARIABLE SRAY 1;\n')
    f.write('VARIABLE SRDE 1;\n')
    f.write('\n')
    f.write('VARIABLE MRESOL_P1 1;\n')
    f.write('VARIABLE MRESOL_P2 1;\n') 
    f.write('VARIABLE MRESOL_P3 1;\n')
    f.write('VARIABLE RESOLU_P2 1;\n')
    f.write('\n')
    f.write('VARIABLE NX 1;\n')
    f.write('VARIABLE NA 1;\n')
    f.write('VARIABLE NY 1;\n')
    f.write('VARIABLE NB 1;\n')
    f.write('VARIABLE NE 1;\n')
    f.write('\n')
    f.write('VARIABLE N1 1;\n')
    f.write('VARIABLE N2 1;\n')
    f.write('vARIABLE N3 1;\n')
    f.write('VARIABLE N4 1;\n')
    f.write('VARIABLE N5 1;\n')
    f.write('VARIABLE READPARA 1;\n')
    f.write('\n')
    f.write('\n')
    f.write('{**************** READ_RAY ************}\n')
    f.write('\n')
    f.write('PROCEDURE READ_RAY;\n')
    f.write('\n')
    f.write('VARIABLE SAVEMAP 10000 8;\n')
    f.write('SM SAVEMAP;\n')
    f.write('UM; CR;\n')
    f.write('\n')
    f.write('{N1=1 {# of rays: 3} -> -1 0 +1}\n')
    f.write('{N1=2 {# of rays: 5} -> -2 -1 0 +1 +2}\n')
    f.write('{N1=3 {# of rays: 7} -> -3 -2 -1 0 +1 +2 +3}\n')
    f.write('{N1=4 {# of rays: 9} -> -4 -3 -2 -1 0 +1 +2 +3 +4}\n')
    f.write('\n')
    f.write('N1:=1; N2:=3; N3:=1;N4:=1; N5:=2;\n')
    f.write('\n')
    f.write('LOOP NX 1 2*N1+1;LOOP NA 1 2*N2+1;LOOP NY 1 2*N3+1;LOOP NB 1 2*N4+1;\n')
    f.write('LOOP NE 1 2*N5+1;\n')
    f.write('\n')
    f.write('SRXX:= XX*(NX-(N1+1))/N1;\n')
    f.write('SRAX:= AX*(NA-(N2+1))/N2;\n')
    f.write('SRYY:= YY*(NY-(N3+1))/N3;\n')
    f.write('SRAY:= AY*(NB-(N4+1))/N4;\n')
    f.write('SRDE:= DE*(NE-(N5+1))/N5;\n')
    f.write('\n')
    f.write('IF (((NA-(N2+1))/N2)^2+((NB-(N4+1))/N4)^2+((NE-(N5+1))/N5)^2)<1.01;\n')
    f.write('SR SRXX SRAX SRYY SRAY 0 SRDE 0 0 1;\n')
    f.write('ENDIF;\n')
    f.write('ENDLOOP;ENDLOOP;ENDLOOP;ENDLOOP;ENDLOOP;\n')
    f.write('\n')
    f.write('AM SAVEMAP;\n')
    f.write('ENDPROCEDURE;\n')
    f.write('\n')
    f.write('{**************** RECOIL_LINE ************}\n')
    f.write('\n')
    f.write('PROCEDURE RECOIL_BL;\n')
    f.write('VARIABLE X_MAP 1;\n')
    f.write('\n')
    f.write('{Section 1}\n')
    f.write('\n')
    f.write('FR 3; {Fringe field flag}\n')
    f.write('\n')
    f.write('DL 0.80+0.000106; {z}                                         {DL1}\n')
    f.write('\n')
    f.write('TA -0.0299 -0.0435; {Pitch Yaw}\n')
    f.write('RA 0.0035; {Roll}\n')
    f.write('SA 0.000290 -0.000002; {x, y}\n')
    f.write('M5 0.250 {0:.6f} 0.004679 0 -0.00318 0 0.055;'.format(q1s))   
    f.write('                               {Q1+Hex}\n')
    f.write('SA -0.000290 0.000002;\n')
    f.write('RA -0.0035;\n')
    f.write('TA 0.0299 0.0435;\n')
    f.write('\n')
    f.write('DL 0.19+0.00105-0.000106; {z}                                 {DL2}\n')
    f.write('\n')
    f.write('TA 0.0496 0.0620; {Pitch, Yaw}\n')
    f.write('RA 0.0252; {Roll}\n')
    f.write('SA 0.0002 0.000068; {x, y}\n')
    f.write('MQ 0.300-0.0021 {0:.6f} 0.068;'.format(q2s))   
    f.write('                               {Q2}\n')
    f.write('SA -0.0002 -0.000068; {-x, -y}\n')
    f.write('RA -0.0252; {-Roll}\n')
    f.write('TA -0.0496 -0.0620; {-Pitch, -Yaw}\n')
    f.write('\n')
    f.write('DL 0.58-0.000125+0.00105;                                     {DL3}\n')
    f.write('\n')
    f.write('TA 0.0094 0.0108; {Pitch, Yaw}\n')
    f.write('RA 0.0058; {Roll}\n')
    f.write('SA -0.000036 0.000112; {x, y}\n')
    f.write('MC 1.25 22.5+0.01145 0.030 B1N B1S1 B1S2 7;                   {B1}\n')  
    f.write('SA 0.000036 -0.000112; {-x, -y}\n') 
    f.write('RA -0.0058; {-Roll}\n')
    f.write('TA -0.0094 -0.0108; {-Pitch, -Yaw}\n')
    f.write('\n')
    f.write('DL 1.00-0.000125-0.00006-0.000073;                            {DL4}\n')
    f.write('\n')
    f.write('TA 0.0044 0.0044; {Pitch, Yaw}\n')
    f.write('RA 0.0012; {Roll}\n')
    f.write('SA -0.000007 0.000147; {x, y}\n')
    f.write('MC 1.25 22.5+0.0055+0.0066 0.030 B1N B2S1 B2S2 7;             {B2}\n')
    f.write('SA 0.000007 0.000147; {-x, -y}\n')
    f.write('RA -0.0012; {-Roll}\n')
    f.write('TA -0.0044 -0.0044; {-Pitch, -Yaw}\n')
    f.write('\n')
    f.write('\n')
    f.write('\n')
    f.write('DL 0.77-0.00006-0.000073;                                     {DL5}\n')
    f.write('\n')
    f.write('DL 0.40;                                                      {DL6}\n')
    f.write('\n')
    f.write('MH 0.26 0.0103064+0.0000001 0.11;                             {HEX1}\n')
    f.write('\n')  
    f.write('DL 0.27+0.00005;                                              {DL7}\n')
    f.write('MQ 0.350-0.0001 {0: .6f} 0.11;'.format(q3s))
    f.write('                               {Q3}\n')
    f.write('DL 0.35+0.00005+0.00165;                                      {DL8}\n')
    f.write('M5 0.350-0.0033 {0: .6f} 0 0 0 0 0.08;'.format(q4s))
    f.write('                               {Q4}\n')
    f.write('DL 0.21+0.00165+0.0017;                                       {DL9}\n')
    f.write('MQ 0.350-0.0034 {0: .6f} 0.06;'.format(q5s))
    f.write('                               {Q5}\n')
    f.write('DL 0.145+0.0017 ;                                             {DL10}\n')  
    f.write('PS 0.01;\n')
    f.write('\n')
    #f.write('WRITE 6 \'ME(1,2) AT FP1 =\' ME(1,2);\n')
    #f.write('WRITE 6 \'ME(1,1) AT FP1 =\' ME(1,1);\n')
    f.write('\n') 
    f.write('{Section 2}\n')
    f.write('\n')
    f.write('DL 0.185 ;                                                    {DL11}\n')
    f.write('\n')
    f.write('DL 0.17-0.00035 ;                                             {DL12}\n')
    f.write('MC 1.25 22.5+0.0321 0.05 B1N B3S1 B3S2 7;                     {B3}\n')
    f.write('DL 0.51-0.00035-0.00192-0.00088;                              {DL13}\n')
    f.write('\n')
    f.write('MC 1.25 22.5+0.0807{-0.012} 0.05 B1N B4S1 B4S2 7;             {B4}\n')
    f.write('DL 0.30-0.00088;                                              {DL14}\n')
    f.write('\n')
    f.write('M5 0.26 0 0.01449-0.003433-0.00071+0.00016 0 0 0 0.12;        {HEX2}\n')
    f.write('DL 0.27;                                                      {DL15}\n')
    f.write('DL 0.27+0.0001;                                               {DL16}\n')
    f.write('MQ 0.34-0.0002 {0: .6f} 0.14;'.format(q6s))
    f.write('                               {Q6}\n')
    f.write('DL 0.20+0.0001-0.00005 ;                                      {DL17}\n')
    f.write('MQ 0.34+0.0001 {0: .6f} 0.13;'.format(q7s))
    f.write('                               {Q7}\n')
    f.write('DL 0.50-0.00005 ;                                             {DL18}\n')
    f.write('\n')
    f.write('{FC 1=dipole, 1/2=entrance/exit, 1/2=magn/elect, a1...a6}\n')
    f.write('\n')
    f.write('FC 1 1 1 -0.16 1.603 -0.0105 0.015 -0.0226 0.0038; {entr,magnet dipole}\n')
    f.write('FC 1 2 1 -0.16 1.603 -0.0105 0.015 -0.0226 0.0038; {exit,magnet dipole}\n')
    f.write('\n')
    f.write('FC 1 1 2 -0.167 1.874 0.246 -0.052 0.0142 0.066;{entr, elect dipole}\n')
    f.write('FC 1 2 2 -0.167 1.874 0.246 -0.052 0.0142 0.066;{exit, elect dipole}\n')
    f.write('\n')
    f.write('CB;\n')
    f.write('WC 7.0 7.0 2.365  0.11 NE1 NM1 4 ;                            {WF1}\n')
    f.write('CB;\n')
    f.write('FD;\n')
    f.write('FR 3;\n')
    f.write('\n')
    f.write('\n')
    f.write('DL 0.50; {THIS IS PART OF WF1!}                               {DL19}\n')
    f.write('\n')
    f.write('M5 0.26 0 -0.01251*(0.09/0.11)^2+0.0000 0 0 0 {0.11}0.09;     {HEX3}\n')
    f.write('DL 0.28;                                                      {DL20}\n')
    f.write('\n')
    f.write('M5 0.26 0 0 (0.040233+0.01-0.01700)-0.00195 0 0 {0.07}0.09;   {OCT1}\n')
    f.write('DL 1.75;                                                      {DL21}\n')
    f.write('PS 0.005;\n')
    f.write('\n')
    f.write('{Section 3}\n')
    f.write('\n')
    f.write('{{\n')
    f.write('DL    0.872;                                                  {DL22}\n')
    f.write('\n')
    f.write('{M5 0.25 -0.15032+0.0017{+0.00325} -0.0001 0 0 0 0.05;        {Q8}}\n')
    f.write('M5 0.25 -0.15032+0.0017+0.0004 0 -0.0001 0 0 0.05;            {Q8_oct}\n')
    f.write('\n')
    f.write('DL 0.395;                                                     {DL23}\n')
    f.write('\n')
    f.write('{M5 0.30 0.23438 0.00015 0 0 0 0.07;                          {Q9}}\n')
    f.write('M5 0.30 0.23438 0 0.00015 0 0 0.07;                           {Q9_oct}\n')
    f.write('\n')
    f.write('DL 0.36;                                                      {DL24}\n')
    f.write('\n')
    f.write('MC 1.25 42.5 0.03+0.006 B1N B5S1 B5S2 7;                      {B5}\n')
    f.write('DL 0.35;                                                      {DL25}\n')
    f.write('MC 1.25 42.5 0.03+0.006 B1N B6S1 B6S2 7;                      {B6}\n')
    f.write('DL 0.83 ;                                                     {DL26}\n')
    f.write('\n')
    f.write('{MQ 0.26 -0.03250 0.09 ;                                      {Q10}}\n')
    f.write('{M5 0.26 -0.03250 0.00005 0 0 0 0.09;                         {Q10}}\n')
    f.write('M5 0.26 -0.03250-0.00050-0.00067 0 0.00005 0 0 0.09;          {Q10_oct}\n')
    f.write('\n')
    f.write('DL 0.65;                                                      {DL27}\n')
    f.write('{MQ 0.34 0.1616-0.00002 0.12;                                 {Q11}}\n')
    f.write('M5 0.34 0.1616-0.00002 0 0.000125 0 0 0.12;                   {Q11_oct}\n')
    f.write('\n')
    f.write('DL 1.0 ;                                                      {DL28}\n')
    f.write('\n')
    f.write('FC 1 1 1 -0.16 1.603 -0.0105 0.015 -0.0226 0.0038; {entr,magnet dipole}\n')
    f.write('FC 1 2 1 -0.16 1.603 -0.0105 0.015 -0.0226 0.0038; {exit,magnet dipole}\n')
    f.write('\n')
    f.write('FC 1 1 2 -0.167 1.874 0.246 -0.052 0.0142 0.066;{entr, elect dipole}\n')
    f.write('FC 1 2 2 -0.167 1.874 0.246 -0.052 0.0142 0.066;{exit, elect dipole}\n')
    f.write('\n')
    f.write('CB;\n')
    f.write('WC 7.0*(1+0.0000)  7.0 2.365  0.11 NE1 NM1 4 ;                {WF2}\n')
    f.write('CB;\n')
    f.write('FD;\n')
    f.write('FR 3;\n')
    f.write('\n')
    f.write('DL 4.60;                                                      {DL29}\n')
    f.write('PS 0.035;   {FP3}\n')
    f.write('}}\n')
    f.write('{{{\n')
    f.write('\n')
    f.write('{SECTION 4}\n')
    f.write('\n')
    f.write('DL 0.25+0.00115;                                              {DL30}\n')
    f.write('MQ 0.3-0.0023 -0.1820 0.07;                                   {Q12}\n')
    f.write('DL 0.30+0.05+0.00115-0.0004;                                  {DL31}\n') 
    f.write('MQ 0.3+0.0008 0.1910 0.05;                                    {Q13}\n')
    f.write('DL 0.66-0.0004;                                               {DL32}\n')
    f.write('\n')
    f.write('MC 1.25 55. 0.03 B1N B7S1 B7S2 4;                             {B7}\n')
    f.write('DL 0.68;                                                      {DL33}\n')
    f.write('MC 1.25 55. 0.03 B1N B8S1 B8S2 4;                             {B8}\n')
    f.write('\n')
    f.write('DL 0.86+0.00025;                                              {DL34}\n')
    f.write('MQ 0.3-0.0005 0.1290 0.05;                                    {Q14}\n')
    f.write('DL 0.45+0.00025-0.0006;                                       {DL35}\n')
    f.write('MQ 0.3+0.0012 -0.138 0.05;                                    {Q15}\n')
    f.write('DL 1.21-0.0006;                                               {DL36}\n')
    f.write('\n')
    f.write('PS 0.020;\n')
    f.write('DL 1.1;                                                       {DL37}\n')
    f.write('\n')
    f.write('PS 0.015;                                                     {Detector}\n')
    f.write('{MASS_SEPARATION MRES;}\n')
    f.write('\n')
    #f.write('MRESOL_P3 := ABS(ME(1,7))/(2*XX*ME(1,1));\n')
    #f.write('WRITE 6 \'Mass Res.Power in Det plane=\' MRESOL_P3;\n')
    f.write('\n')
    #f.write('PM 12 ; {print map to file fort.12}\n')
    f.write('\n')
    #f.write('WRITE 6 \'ME(1,1),ME(1,2),ME(1,6)=\' ME(1,1) ME(1,2) ME(1,6);\n')
    #f.write('WRITE 6 \'M11*M22=\' ME(1,1)*ME(2,2);\n')
    f.write('\n')
    f.write('{}\n')
    f.write('DL 0.4;                                                       {DL38}\n')
    f.write('PS 0.02;           {FP3}\n')
    f.write('\n')
    f.write('DL 0.50;                                                      {DL39}\n')
    f.write('{}\n')
    f.write('}}}\n')
    f.write('ENDPROCEDURE;\n')
    f.write('\n')
    
    f.write('{############################################################################}\n')
    f.write('{########################### DEFINITION AND COMMAND #########################}\n')
    f.write('{############################################################################}\n')
    
    
    f.write('{OV 2 3 3 ;}\n') 
    f.write('{OV 3 3 3 ;}\n')
    f.write('OV 4 3 2 ;\n')
    f.write('\n')
    f.write('WSET 1.5E-3;\n')
    f.write('{Optimized Ey- and By-fields, including higher orders}\n')
    f.write('{in parentheses are the limits for MRES2 = 750}\n')
    f.write('\n')
    f.write('NE1(1) := 0;\n')
    f.write('NE1(2) := -0.01 {y=0mm};\n')
    f.write('{NE1(2) := -0.10{y=35mm};}\n')
    f.write('NE1(3) := 0;\n')
    f.write('NE1(4) := 11. {y=0mm};\n')
    f.write('{NE1(4) := 12. {y=35mm};}\n')
    f.write('\n')
    f.write('NM1(1) := 0;\n')
    f.write('NM1(2) := 0.033{y=0};\n')
    f.write('{NM1(2) := 0.063{y=35};}\n')
    f.write('NM1(3) := 0;\n')
    f.write('NM1(4) := -3.8{y=0};\n')
    f.write('{NM1(4) := -4{y=35};}\n')
    f.write('\n')
    f.write('XX:=0.00075{-0.00075}{+0.0005};\n')
    f.write('AX:=0.025{-0.025};\n')
    f.write('YY:=0.00075{-0.00075}{+0.0005};\n')
    f.write('AY:=0.025{-0.025};\n')
    f.write('DE:= 0.031{-0.031};\n')
    f.write('\n')
    f.write('NMAX:=190;\n')
    f.write('\n')
    f.write('B1S1(1):={0.19}{0.189}{map}{0.189}0.189;\n')
    f.write('B1S1(2):={0.0025}{0.0036}{map}{0.0025}0.0115;\n')
    f.write('B1S1(3) := {0.154}{0.198}{map}{0.252}0.2438;\n')
    f.write('B1S1(4):= {0.78}{0.69}{map}{1.1}0.9504;\n')
    f.write('\n')
    f.write('B1S2(1):={0.15}{0.150}{map}{0.149}0.1504;\n')
    f.write('B1S2(2):={-0.019}{-0.0022}{map}{-0.0332}-0.0278;\n')
    f.write('B1S2(3) := {0.147}{0.209}{map}{0.240}0.1095;\n')
    f.write('B1S2(4):= {0.10}{-1.29}{map}{1.65}0.8070;\n')
    f.write('\n')
    f.write('B2S1(1):={0.115}{0.115}{map}{0.115}{0.115}0.115;\n')
    f.write('B2S1(2):={0.0125}{0.0244}{map}{0.0154}{0.0017}0.0083;\n')
    f.write('B2S1(3) := {0.198}{0.23}{map}{0.28}{0.186}0.180;\n')
    f.write('B2S1(4):= {-40.77}{-42.52}{map}{-40.74}{-40.04}-40.15;\n')
    f.write('\n')
    f.write('B2S2(1):={0.115}{0.115}{map}{0.114}0.115;\n')
    f.write('B2S2(2):={-0.2448}{-0.2449}{map}{-0.2539}-0.2499;\n')
    f.write('B2S2(3) := {1.411}{1.40}{map}{1.526}1.430;\n')
    f.write('B2S2(4):= {37.47}{37.57}{map}{38.96}38.34;\n')
    f.write('\n')
    f.write('{Below enter B3 measured results of final magnet}\n')
    f.write('B3S1(1):={map}0.190;\n')
    f.write('B3S1(2):={map}1.0541;\n')
    f.write('B3S1(3) :={map}-9.475;\n')
    f.write('B3S1(4) :={map}5.54;\n')
    f.write('\n')
    f.write('B3S2(1):={map}0.115;\n')
    f.write('B3S2(2):={map}-0.0499;\n')
    f.write('B3S2(3) :={map}34.09;\n')
    f.write('B3S2(4) :={map}52.87;\n')
    f.write('\n')
    f.write('{\n')
    f.write('B3S1(1):={0.19}{0.185}0.190;\n')
    f.write('B3S1(2):= {1.07}{1.1025}1.054;\n')
    f.write('B3S1(3) := {-9.10}{-7.839}-9.21;\n')
    f.write('B3S1(4) := {0.}{-2.373}5.51;\n')
    f.write('\n')
    f.write('B3S2(1):={0.115}{0.110}0.1141;\n')
    f.write('B3S2(2):={-0.3315}{0.07}{-0.20}{0.0410}-0.159;\n')
    f.write('B3S2(3) := {32.70}{35.747}33.97;\n')
    f.write('B3S2(4) := {89.56}{-90}{-30}{-57.}56.77;\n')
    f.write('}\n')
    f.write('\n')
    f.write('{Start: Report email 8/26/2016, with presence of B3 iron}\n')
    f.write('{{{}}}\n')
    f.write('B4S1(1):={spec}{0.19}{0.191}{0.193}{map}0.191;\n')
    f.write('B4S1(2):={spec}{-0.339}{-0.3394}{-0.2990}{map}-0.3390;\n')
    f.write('B4S1(3) := {spec}{-5.51}{-5.536}{-5.623}{map}-5.464;\n')
    f.write('B4S1(4) := {spec}{-0.84}{-0.48}{-3.19}{map}0.70;\n')
    f.write('\n')
    f.write('B4S2(1):={0.19}{0.190}{0.190}{map}0.190;\n')
    f.write('B4S2(2):={-0.030}{-0.0265}{-0.0279}{map}-0.0287;\n')
    f.write('B4S2(3) := {-0.364}{-0.347}{-0.376}{map}-0.332;\n')
    f.write('B4S2(4) := {-0.15}{0.28}{0.24}{map}0.08;\n')
    f.write('{{{}}}\n')
    f.write('{End: Report email 8/26/2016, with presence of B3 iron}\n')
    f.write('{Start: Final B4 stand alone, Report of 9/14/2016}\n')
    f.write('{{{\n')
    f.write('B4S1(1):=0.190;\n')
    f.write('B4S1(2):= -0.3409;\n')
    f.write('B4S1(3) := -5.540;\n')
    f.write('B4S1(4) := -0.44;\n')
    f.write('\n')
    f.write('B4S2(1):= 0.190;\n')
    f.write('B4S2(2):= -0.0283;\n')
    f.write('B4S2(3) := -0.363;\n')
    f.write('B4S2(4) := 0.13;\n')
    f.write('}}}\n')
    f.write('{End: Final B4 stand alone, Report of 9.14.2016}\n')
    f.write('\n')
    f.write('{{{Opt. Settings of p. 384-385}\n')
    f.write('B5S1(1):=0.189;\n')
    f.write('B5S1(2):= 0.696;\n')
    f.write('B5S1(3):= -0.953;\n')
    f.write('B5S1(4):= -53.;\n')
    f.write('\n')
    f.write('B5S2(1):=-0.172;\n')
    f.write('B5S2(2):= -5.928;\n')
    f.write('B5S2(3):= 0;\n')
    f.write('B5S2(4):= 0;\n')
    f.write('\n')
    f.write('B6S1(1):=0.197;\n')
    f.write('B6S1(2):= -1.66;\n')
    f.write('B6S1(3):= 0.0;\n')
    f.write('B6S1(4):=  0.0;\n')
    f.write('\n')
    f.write('B6S2(1):=0.200;\n')
    f.write('B6S2(2):= -4.00;\n')
    f.write('B6S2(3):= 69.;\n')
    f.write('B6S2(4):= 0.0;\n')
    f.write('}}\n')
    f.write('\n')
    f.write('{DF Design settings}\n')
    f.write('B5S1(1):={0.189}{0.190}0.189;\n')
    f.write('B5S1(2):= {0.696}{0.6994}0.712;\n')
    f.write('B5S1(3):= {-0.953}{-1.220}-0.825;\n')
    f.write('B5S1(4):= {-53.}{-50.499}-53.36;\n')
    f.write('\n')
    f.write('B5S2(1):={-0.172}{-0.172}{-0.179}{-0.181}-0.180;\n')
    f.write('B5S2(2):= {-5.928}{-5.6}{-5.6722}{-5.5503}-5.549;\n')
    f.write('B5S2(3):= {-26.5}{0}{4.978}{3.632}3.288;\n')
    f.write('B5S2(4):= {-940.}{0}{76.733}{25.066}28.82;\n')
    f.write('\n')
    f.write('B6S1(1):={0.197}{0.198}0.198;\n')
    f.write('B6S1(2):= {-1.66}{0.0324}0.023;\n')
    f.write('B6S1(3):= {-50.0}{0}{-0.238}-0.227;\n')
    f.write('B6S1(4):=  {0.0}{-3.543}1.02;\n')
    f.write('\n')
    f.write('B6S2(1):={0.200}{0.211}0.203;\n')
    f.write('B6S2(2):= {-4.00}{-4.5041}-4.13;\n')
    f.write('B6S2(3):= {69.}{67.528}67.40;\n')
    f.write('B6S2(4):= {0.0}{103.325}52.95;\n')
    f.write('\n')
    f.write('B7S1(1):= {0}{map}0; \n')
    f.write('B7S1(2):= {0.0036}{map}0.0115;\n')
    f.write('B7S1(3):= {0.032}{map}0.066;\n')
    f.write('B7S1(4):= {1.98}{map}3.44;\n')
    f.write('\n')
    f.write('B7S2(1):= {0}{map}0;\n')
    f.write('B7S2(2):= {0.003}{map}0.0228;\n')
    f.write('B7S2(3):= {0.032}{map}0.036;\n')
    f.write('B7S2(4):= {1.98}{map}-1.02;\n')
    f.write('\n')
    f.write('B8S1(1):= {0}{map}0;\n')
    f.write('B8S1(2):= {0}{map}0.0115;\n')
    f.write('B8S1(3):= {0}{map}0.066;\n')
    f.write('B8S1(4):= {0}{map}3.44;\n')
    f.write('\n')
    f.write('B8S2(1):= {0}{map}0;\n')
    f.write('B8S2(2):= {0}{map}0.0228;\n')
    f.write('B8S2(3):= {0}{map}0.036;\n')
    f.write('B8S2(4):= {0}{map}-1.02;\n')
    f.write('\n')
    f.write('RP 206 66*PARA(1) 21*PARA(2) ;\n')
    f.write('\n')
    f.write('UM;\n')
    f.write('CR {Clears all rays}; UM;\n')
    f.write('\n')
    f.write('{N1=1 {# of rays: 3} -> -1 0 +1} \n')                   
    f.write('{N1=2 {# of rays: 5} -> -2 -1 0 +1 +2} \n')             
    f.write('{N1=3 {# of rays: 7} -> -3 -2 -1 0 +1 +2 +3} \n')       
    f.write('{N1=4 {# of rays: 9} -> -4 -3 -2 -1 0 +1 +2 +3 +4} \n') 
    f.write('\n')
    f.write('N1:=1;N2:=3;N3:=1;N4:=1;N5:=2;\n')
    f.write('\n')
    f.write('LOOP NX 1 2*N1+1;LOOP NA 1 2*N2+1;LOOP NY 1 2*N3+1;LOOP NB 1 2*N4+1;\n')
    f.write('LOOP NE 1 2*N5+1;\n')
    f.write('\n')
    f.write('SRXX:= XX*(NX-(N1+1))/N1;\n')
    f.write('SRAX:= AX*(NA-(N2+1))/N2;\n')
    f.write('SRYY:= YY*(NY-(N3+1))/N3;\n')
    f.write('SRAY:= AY*(NB-(N4+1))/N4;\n')
    f.write('SRDE:= DE*(NE-(N5+1))/N5;\n')
    f.write('\n')
    f.write('IF (((NA-(N2+1))/N2)^2+((NB-(N4+1))/N4)^2+((NE-(N5+1))/N5)^2)<1.01;\n')
    f.write('SR SRXX SRAX SRYY SRAY 0 SRDE 0 0 1;\n')
    f.write('ENDIF;\n')
    f.write('ENDLOOP;ENDLOOP;ENDLOOP;ENDLOOP;ENDLOOP;\n')
    f.write('\n')
    f.write('SR 0 0 0 0 0 0.035 0 0 2;\n')
    f.write('SR 0 0 0 0 0 -0.035 0 0 2;\n')
    f.write('{{{\n')
    f.write('{Beam A = 15, Target(4He) A=4, dE=1/15=0.066667, dM=1/19=0.052632}\n')
    f.write('SR 0.0     0.0  0.0 0.0 0.0  0.066667      -0.052632 0.0 1 ;{Beam}\n')
    f.write('SR 0.0     0.01  0.0 0.0 0.0  0.066667      -0.052632 0.0 1 ;{Beam}\n')
    f.write('SR 0.0     -0.01  0.0 0.0 0.0  0.066667      -0.052632 0.0 1 ;{Beam}\n')
    f.write('\n')
    f.write('{Beam A = 20, Target(4He) A=4, dE=1/20=0.05, dM=1/24=0.041667}\n')
    f.write('SR 0.0     0.0  0.0 0.0 0.0  0.05      -0.041667 0.0 1 ;{Beam}\n')
    f.write('SR 0.0     0.01  0.0 0.0 0.0  0.05      -0.041667 0.0 1 ;{Beam}\n')
    f.write('SR 0.0     -0.01  0.0 0.0 0.0  0.05      -0.041667 0.0 1 ;{Beam}\n')
    f.write('SR 0.0     0.0  0.0 0.0 0.0  0.0     0 1/15 4 ;\n')
    f.write('SR 0.0     0.0  0.0 0.0 0.0  0.0     0 1/33 5 ;\n')
    f.write('SR 0.0     0.0  0.0 0.0 0.0  0.002857      -0.002849 0.0 2 ;{Beam}\n')
    f.write('\n')
    f.write('{Beam A = 350, Target(H) A=1, dE=1/350=0.002857, dM=1/351=0.002849}\n')
    f.write('SR 0.0     0.0  0.0 0.0 0.0  0.001667     -0.001664 0.0 3 ;{Beam}\n')
    f.write('\n')
    f.write('{Beam A = 600, Target(H) A=1, dE=1/600=0.001667, dM=1/601=0.001664}\n')
    f.write('SR 0.0     0.0  0.0 0.0 0.0  0.001111     -0.00111 0.0 4 ;{Beam}\n')
    f.write('\n')
    f.write('{Beam A = 900, Target(H) A=1, dE=1/900=0.00125, dM=1/901=0.001248}\n')
    f.write('\n')
    f.write('SR 0.0     0  0.001 0.0 0.0  0.0     0.0 0.0 3 ;{red}\n')
    f.write('SR 0.0     0  -0.001 0.0 0.0  0.0     0.0 0.0 3 ;{red}\n')
    f.write('SR 0.001    0  0.00000 0.0 0.0  0.0     0.0 0.0 3 ;{red}\n')
    f.write('SR -0.001    0  0.00000 0.0 0.0  0.0     0.0 0.0 3 ;{red}\n')
    f.write('\n')
    f.write('SR 0.0     0.0  0.0 0.0 0.0  0.0005      -0.00049975 0.0 1 ;{Beam}\n')
    f.write('{Beam A = 2000, Target(H) A=1, dE=1/2000=0.0005, dM=1/2001=0.00049975}\n')
    f.write('}}}\n')
    f.write('\n')
    f.write('BP;\n')
    f.write('RECOIL_BL;\n')
    f.write('EP;\n')
    f.write('\n')
    #f.write('PG -1 -2;\n')
    #f.write('PP -10 0. 0.;\n')
    #f.write('PP -10 0. 90.;\n')
    f.write('\n')
    f.write('READ_RAY;\n')
    f.write('WW(1):= VMAX(RAY(1));\n')
    #f.write('WRITE 6 \'VMAX_SIZE)=\' WW(1);\n')
    f.write('WW(2):= VMIN(RAY(1));\n')
    #f.write('WRITE 6 \'VMIN_SIZE)=\' WW(2);\n')
    f.write('WV:=(WW(1)-WW(2))/2;\n')
    #f.write('WRITE 6 \'HALF Image size\' WV;\n')
    #f.write('\n')
    #f.write('CENTER := (WW(1)+WW(2))/2;\n')
    #f.write('WRITE 6 \'Image Center\' CENTER;\n')
    #f.write('{\n')
    #f.write('WRITE 6 \'RAY1\' RAY(1);\n')
    #f.write('}\n')
    #f.write('MRESOL_P1 := ABS(ME(1,7))/(2*XX*ME(1,1)); \n')
    #f.write('WRITE 6 \'Mass Res.Power at FP2 or FP3 =\' MRESOL_P1;\n')
    #f.write('\n')
    #f.write('WRITE 6 \'Resol.by max.ray =\' ABS(ME(1,7))/(2*WV);\n')
    #f.write('WRITE 6 \'ME(1,1),ME(1,2),ME(1,6),ME(1,7)=\' ME(1,1) ME(1,2) ME(1,6) ME(1,7);\n')
    #f.write('WRITE 6 \'M11*M22=\' ME(1,1)*ME(2,2);\n')
    #f.write('OPENF 99 \'temp-results\' \'NEW\';\n')
#    f.write('OPENF 99 \'{0}\' \'NEW\';\n'.format(tempresFilename))
    f.write('WRITE 6 -ABS(ME(1,7))/(2*WV);  \n')
    f.write('WRITE 6 ABS(ME(1,6));  \n') # min
    f.write('WRITE 6 ABS(ME(2,6));  \n') # min
    f.write('WRITE 6 ABS(ME(1,2));  \n') # min
    f.write('\n')
    f.write('ENDPROCEDURE ;\n')
    f.write('RUN ;\n')
    f.write('END ;\n') 
    f.close()
    
    #Removing files from older runs
#    failure, output = commands.getstatusoutput(cmd)
    
    #Run file
    cmd = 'cosy ' + cosyFilename
#    failure, output = commands.getstatusoutput(cmd)
    output = commands.run(['cosy',cosyFilename], capture_output=True)
    stripped = output.stdout.strip()
#    print(stripped.split())
    resol = (stripped.split())
    for i in range(len(resol)):
        resol[i] = float(resol[i])
    print(resol)            
    commands.run(['rm','-f',cosyFilename])
    commands.run(['rm','-f',lisFilename])

#    try:
#        f = open(tempresFilename,'r')  # Opening temp file with results
#        resol = f.readline()
#        f.close()
#    except (OSError, IOError) as e:
#        print('q1 = %.6f q2 = %.6f q3 = %.6f q4 = %.6f q5 = %.6f q6 = %.6f q7 = %.6f'  %(q1s, q2s, q3s, q4s, q5s, q6s, q7s) )   
#        resol = 0
#    f = open('results.txt','a')  # Writing results file: magnet field, resolution
#    f.write( '{0:d} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.1f}\n' .format(b, q1s, q2s, q3s, q4s, q5s, q6s, q7s, float(resol)) )
#    f.close()

#    print('Finished bee %d'%b)
    return (resol)
    

if __name__ == "__main__":
    print(cosyrun(qNew))



