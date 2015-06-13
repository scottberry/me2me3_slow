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

_OBJ = random.o modifications.o gillespie.o results.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: $(STATLIB) $(OBJ) me2me3 Tests ConstTimeInterpolate

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

Tests: $(ODIR)/FunctionTests.o $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

ConstTimeInterpolate: $(ODIR)/ConstTimeInterpolate.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS) $(LFLAGS) $(IFLAGS)

.PHONY: clean empty nores

clean:
	rm -f $(LDIR)/*.o $(ODIR)/*.o *~ $(IDIR)/*~ 

empty:
	rm -f $(LDIR)/*.o $(ODIR)/*.o *~ $(IDIR)/*~ $(STATLIB) me2me3 Tests ConstTimeInterpolate

nores:
	rm -f *.txt
