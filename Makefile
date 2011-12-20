.KEEP_STATE:
SHELL=/bin/sh

#
#	define version of c compiler and linker
#
CC=	g++
LINK=	g++
#
#	define compiler flags
#
CFLAGS	= -Wall -Wextra	-Wshadow -Wpointer-arith -Wcast-qual \
	-Wcast-align -Wwrite-strings -fshort-enums -fno-common \
	-g -O2 -DHAVE_INLINE
LDLIBS	= -lgsl -lcblas -latlas -lm

DEPS	= zb_def.hh wz_def.hh gwz_def.hh s_def.hh coulomb_def.hh \
	dielectric_def.hh well_def.hh pseudopotential_def.hh \
	hamiltonian.hh matrix_term.hh

ODIR	= obj
_OBJ	= main.o
OBJ	= $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN	= eivals.bin

GENERATED = $(OBJ) $(BIN) s_def.hh zb_def.hh wz_def.hh gwz_def.hh \
	coulomb_def.hh dielectric_def.hh well_def.hh \
	pseudopotential_def.hh sc_def.hh

.PHONY	:	all
all	:	$(BIN)

run	:	$(BIN)
	./$^

$(BIN)	:	$(OBJ)
	$(LINK) -o $@ $^ $(CFLAGS) $(LDLIBS)

.PHONY	:	objs
objs	:	$(OBJ)

$(ODIR)/%.o : %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

%_def.hh :	%.txt math2matrix.pl
	./math2matrix.pl < $< > $@

.PHONY	:	gdb
gdb	:	$(BIN)
	gdb ./$(BIN)

.PHONY	:	valgrind
valgrind :	$(BIN)
	valgrind --trace-children=yes --leak-check=full \
		--show-reachable=yes ./$(BIN)

.PHONY	:	clean
clean	:
	-rm -f $(OBJ)

.PHONY	:	distclean
distclean :
	-rm -f $(GENERATED)
