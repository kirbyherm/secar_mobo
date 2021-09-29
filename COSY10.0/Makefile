# Makefile for COSY on Linux with ifort of Intel Fortran
#
# Definition of FC, FFLAGS, LINK and LIBS depends on each local system.
# We appreciate receiving information on your local makefile.
#
#
# Default Makefile Notes:  
# If compiling with PGPlot, the corresponding LIBS line needs to be uncommented
#
# FFlags are primarily made for verification and optimization  
#  -m32 is required on 64 bit systems
#  -fp-model strict is for verified numerics
#  -f77rtl is probably not required
#        (and will likely be deprecated in future versions of ifort)

FC = ifort 
FFLAGS= -m64 -f77rtl -fp-model strict 
LINK = ifort
LIBS = 
#LIBS = -L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11

OBJ = foxy.o dafox.o foxfit.o foxgraf.o
OBJSH = foxy2.o dafox.o foxfit.o foxgraf.o

all: cosy cosysh

cosy: $(OBJ)
	$(LINK) -o cosy $(OBJ) $(LIBS)

cosysh: $(OBJSH)
	$(LINK) -o cosysh $(OBJSH) $(LIBS)

clean:
	rm -f cosy cosysh $(OBJ) $(OBJSH)
