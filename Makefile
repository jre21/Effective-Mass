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
	-fno-var-tracking -g -O3 -DHAVE_INLINE
LDLIBS	= -lgsl -lcblas -latlas -lm

DEFNS	= exp_zb_def.hh exp_wz_def.hh exp_gwz_def.hh exp_coulomb_def.hh \
	exp_dielectric_def.hh exp_well_def.hh exp_lcz_def.hh \
	exp_overlap_def.hh \
#	gauss_zb_def.hh gauss_wz_def.hh gauss_gwz_def.hh \
#	gauss_coulomb_def.hh gauss_dielectric_def.hh gauss_well_def.hh \
#	gauss_lcz_def.hh gauss_hgh_def.hh gauss_overlap_def.hh


DEPS	= hamiltonian.hh matrix_term.hh enums.hh defs.hh

ODIR	= obj
_OBJ	= main.o hamiltonian.o matrix_term.o \
	matrix_term_exp.o enums.o # matrix_term_gauss.o
OBJ	= $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN	= evals.bin

GENERATED = $(OBJ) $(BIN)

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

$(ODIR)/matrix_term_%.o : matrix_term_%.cc matrix_term.hh defs.hh $(DEFNS)
	$(CC) -c -o $@ $< $(CFLAGS)

# regenerate all the _def.hh files
# For some reason, the above command errors out rather than generating
# them when they don't exist.
.PHONY	:	defns
defns	:	$(DEFNS)

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
