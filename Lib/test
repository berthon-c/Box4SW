# ----------------------------------------
# makefile pour la compilation de Box4SW
# ----------------------------------------
#
# --- SUFFIXES ---
.SUFFIXES: .out .o .f .h .f90

CF       = gfortran

FFLAGS   =  -O2 \
        -Wall -Wcharacter-truncation -Wunderflow \
        -ffpe-trap=invalid,zero,overflow,underflow \
        -fbounds-check

#FFLAGS   =  -O2

LLFLAGS= $(FFLAGS)
#PVMLIBS = $(HOME)/pvm3/lib/ALPHA
GEOM= ../GEOM/
#FFLAGS    = 
# --- OBJETS ---Flux.o 
MODULES0=Var_Types.o

MODULES1=Var_Types.o readmesh.o IO.o Flux_Numeriques.o Topographie.o \
         Source.o Conditions_limites.o

MODULES_ALL = Var_Types.o readmesh.o IO.o Initialisation.o SegmtsGeom.o \
        Flux_Numeriques.o MvGeomparameters.o Physique.o Topographie.o \
	Source.o Conditions_limites.o

RM = rm
RCP = rcp
EXE = ../
REP = ../SRC/
#
#.f90.o:   ;   $(F90) $(FFLAGS) -c $*.f90
#

# --- EXECUTABLE ---
PG_VF:  mes  $(MODULES_ALL) Un_Pas_Temps.o Main.o 
	$(CF)  $(LLFLAGS)  -o  $(EXE)sw_box $(MODULES_ALL) \
	Un_Pas_Temps.o Main.o 

mes:
	@echo "  ......................................................."
	@echo "  Vous avez sollicite la compilation du code Box4SW "
	@echo "  ......................................................."

clean:  
	@echo "  Netoyage en cours *.o *.mod"
	$(RM) -f *.o *.mod  *~
	$(RM) ../sw_box

#  --------------------------------------------------------------------------
#  dependances de compilation  (a compiler) : (quand ceux ci sont compiles)
#  --------------------------------------------------------------------------
Var_Types.o : $(REP)Var_Types.f90
	$(CF)  $(FFLAGS)  -c $(REP)Var_Types.f90

Physique.o : $(REP)Physique.f90 $(MODULES0)
	$(CF) $(FFLAGS) -c $(REP)Physique.f90

Topographie.o : $(REP)Topographie.f90 $(MODULES0)
	$(CF) $(FFLAGS) -c $(REP)Topographie.f90

Conditions_limites.o : $(REP)Conditions_limites.f90 $(MODULES0) Physique.o
	$(CF) $(FFLAGS) -c $(REP)Conditions_limites.f90

Source.o : $(REP)Source.f90 $(MODULES0) Physique.o
	$(CF) $(FFLAGS) -c $(REP)Source.f90

IO.o : $(REP)IO.f90 $(MODULES0) Physique.o
	$(CF) $(FFLAGS) -c $(REP)IO.f90

Initialisation.o : $(REP)Initialisation.f90 $(MODULES0) Physique.o
	$(CF) $(FFLAGS) -c $(REP)Initialisation.f90

Flux_Numeriques.o : $(REP)Flux_Numeriques.f90 $(MODULES0) Physique.o
	$(CF)  $(FFLAGS) -c $(REP)Flux_Numeriques.f90

Main.o : $(REP)Main.f90 $(MODULES_ALL) Un_Pas_Temps.o 
	$(CF)  $(FFLAGS) -c $(REP)Main.f90 

Un_Pas_Temps.o : $(REP)Un_Pas_Temps.f90 $(MODULES1) 
	$(CF)  $(FFLAGS) -c $(REP)Un_Pas_Temps.f90

readmesh.o : $(GEOM)readmesh.f90 $(MODULES0) 
	$(CF) $(FFLAGS) -c $(GEOM)readmesh.f90

SegmtsGeom.o : $(GEOM)SegmtsGeom.f90 $(MODULES0) 
	$(CF) $(FFLAGS) -c $(GEOM)SegmtsGeom.f90

Geomparameters.o : $(GEOM)Geomparameters.f90 $(MODULES0) 
	$(CF) $(FFLAGS) -c $(GEOM)Geomparameters.f90

MvGeomparameters.o : $(GEOM)MvGeomparameters.f90 $(MODULES0) 
	$(CF) $(FFLAGS) -c $(GEOM)MvGeomparameters.f90



