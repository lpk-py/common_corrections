#
IASQ=emtvelin.o
TAUL=libtau.o libsun.o
#
#  path of executables
BIN = 
#
#  resampling of model 
REM = remodl.o
remodlv: $(REM) $(TAUL) $(IASQ)
	gfortran -std=legacy -o $(BIN)remodlv $(REM) $(IASQ) $(TAUL)
#
#  set branches  program
SBR = setbrn.o
setbrn: $(SBR) $(TAUL) 
	gfortran -std=legacy -o $(BIN)setbrn $(SBR) $(TAUL)
#
# travel time branches program
TBR =  ttimes.o       
ttimes: $(TBR) $(TAUL)
	gfortran -std=legacy -o $(BIN)ttimes $(TBR) $(TAUL)
#
