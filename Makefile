SRC=./src
BIN=./bin
OBJ=./obj
TST=./tst

EXECUTABLE=extract_pod_sirovich
EXECUTABLE2=project_pod_sirovich

SRCS:= declarations.f90 input.f90 output.f90 pierce.f90
OBJS:= $(SRCS:%.f90=%.o)

LINK=-c
CF=gfortran

ifeq ($(CF), gfortran)
FFLAGS=-I $(OBJ) -J $(OBJ) -std=f2008 -cpp -O2 -openmp
endif
ifeq ($(CF), nagfor)
FFLAGS=-I $(OBJ) -mdir $(OBJ) -f2008 -C=all -info -colour -fpp  -g -gline
endif

LDFLAGS=-llapack -openmp

build: $(OBJS) sirovich.o
	$(CF) $(addprefix $(OBJ)/,$^)  $(LDFLAGS) -o $(BIN)/$(EXECUTABLE)

projection: $(OBJS) projection.o
	$(CF) $(addprefix $(OBJ)/,$^)  $(LDFLAGS) -o $(BIN)/$(EXECUTABLE2)

$(OBJS):
	$(CF) $(FFLAGS) $(LINK) $< -o $(OBJ)/$@

sirovich.o:
	$(CF) $(FFLAGS) $(LINK) $< -o $(OBJ)/$@

projection.o:
	$(CF) $(FFLAGS) $(LINK) $< -o $(OBJ)/$@

depend: .depend

.depend: $(addprefix $(SRC)/,$(SRCS)) $(SRC)/sirovich.f90 $(SRC)/projection.f90
	rm -f ./.depend
	gfortran -I $(OBJ) -J $(OBJ) -cpp -M $^  >> ./.depend

clean:
	rm -f $(BIN)/$(EXECUTABLE) $(OBJ)/*.o $(OBJ)/*.mod $(TST)/*.bin $(TST)/*.asc ./.depend

try: build
	make && cd $(TST) && tar -xvf test_fields.tar.gz && sed 's|CURRENT_DIR|'`pwd`'|' input.dat | .$(BIN)/$(EXECUTABLE) && rm field*.bin

retry: build
	make && cd $(TST) && tar -xvf test_fields.tar.gz && sed 's|CURRENT_DIR|'`pwd`'|' input2.dat | .$(BIN)/$(EXECUTABLE) && rm field*.bin

proj: build
	make projection && cd $(TST) && tar -xvf test_fields.tar.gz && sed 's|CURRENT_DIR|'`pwd`'|' input.dat | .$(BIN)/$(EXECUTABLE2) && rm field*.bin && gnuplot plot.gnuplot && rm proj-*comparison-F.asc


include .depend
