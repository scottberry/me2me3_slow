CC = gcc
IDIR = ./include
ODIR = ./obj
LDIR = ./lib

CFLAGS = -O2 -Wall
IFLAGS = -I$(IDIR) -I/usr/local/include

LIBS = -lm -lgsl -lgslcblas -lscottsmatrices	
LFLAGS = -L$(LDIR)
STATLIB = $(LDIR)/libscottsmatrices.a

_DEPS = definitions.h scottsmatrices.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = random.o modifications.o nonprocessive.o gillespie.o results.o parse.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJPROCMETH = random.o modifications.o processive_methylation.o gillespie.o results.o parse.o
OBJPROCMETH = $(patsubst %,$(ODIR)/%,$(_OBJPROCMETH))

_OBJPROCDEMETH = random.o modifications.o processive_demethylation.o gillespie.o results.o parse.o
OBJPROCDEMETH = $(patsubst %,$(ODIR)/%,$(_OBJPROCDEMETH))

all: $(STATLIB) $(OBJ) me2me3 Dynamic HistoneTurnover ProcMeth ProcDemeth Silac Tests ConstTimeInterpolate

# make libscottsmatrices object file
$(LDIR)/scottsmatrices.o: $(LDIR)/scottsmatrices.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(IFLAGS)

# compile libscottsmatrices.a
$(STATLIB): $(LDIR)/scottsmatrices.o
	ar -rcs $@ $^

# make object files for .c files in working directory
$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(IFLAGS)

# make executables
me2me3: $(ODIR)/Main.o $(STATLIB) $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

Dynamic: $(ODIR)/Dynamic.o $(STATLIB) $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

HistoneTurnover: $(ODIR)/HistoneTurnover.o $(STATLIB) $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

ProcMeth: $(ODIR)/Main.o $(STATLIB) $(OBJPROCMETH)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

ProcDemeth: $(ODIR)/Main.o $(STATLIB) $(OBJPROCDEMETH)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

Silac: $(ODIR)/Silac.o $(STATLIB) $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

Tests: $(ODIR)/FunctionTests.o $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

ConstTimeInterpolate: $(ODIR)/ConstTimeInterpolate.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

.PHONY: clean empty nores

clean:
	rm -f $(LDIR)/*.o $(ODIR)/*.o *~ $(IDIR)/*~ 

empty:
	rm -f $(LDIR)/*.o $(ODIR)/*.o *~ $(IDIR)/*~ $(STATLIB) *.pdf *.rds me2me3 Tests ConstTimeInterpolate Silac ProcMeth ProcDemeth HistoneTurnover Dynamic

nores:
	rm -f *.txt *.out
