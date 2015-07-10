CC =gcc
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

_OBJ = random.o modifications.o gillespie.o results.o parse.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJP = random.o modificationsProcK27me2.o gillespie.o results.o parse.o
OBJP = $(patsubst %,$(ODIR)/%,$(_OBJP))

all: $(STATLIB) $(OBJ) me2me3 Dynamic HistoneTurnover TranscriptionInhibit ProcK27me2 Silac Tests ConstTimeInterpolate

# make libscottsmatrices object file
$(LDIR)/scottsmatrices.o: $(LDIR)/scottsmatrices.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(IFLAGS)

# compile libscottsmatrices.a
$(STATLIB): $(LDIR)/scottsmatrices.o
	ar -rcs $@ $^

# make object files for .c files in working directory
$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(IFLAGS)

me2me3: $(ODIR)/Main.o $(STATLIB) $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

Dynamic: $(ODIR)/Dynamic.o $(STATLIB) $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

HistoneTurnover: $(ODIR)/HistoneTurnover.o $(STATLIB) $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

TranscriptionInhibit: $(ODIR)/TranscriptionInhibit.o $(STATLIB) $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

ProcK27me2: $(ODIR)/Main.o $(STATLIB) $(OBJP)
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
	rm -f $(LDIR)/*.o $(ODIR)/*.o *~ $(IDIR)/*~ $(STATLIB) *.pdf *.rds me2me3 Tests ConstTimeInterpolate Silac ProcK27me2 TranscriptionInhibit HistoneTurnover

nores:
	rm -f *.txt *.out
