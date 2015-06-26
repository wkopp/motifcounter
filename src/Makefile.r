SHELL=/bin/bash

SRC=motifdistribution.c background.c sequence.c \
    forground.c matrix.c minmaxscore.c scorefunctions.c \
    score1d.c score2d.c compoundpoisson.c \
    simulate.c countdist.c overlap.c posteriorcount.c
OBJDIR=obj

OBJS=$(patsubst %.c,$(OBJDIR)/%.o,$(SRC))
PROG=./mdist
VAL=
INCLUDES=-I../lapack-3.5.0/lapacke/include/
#INCLUDES=#-I../shared
LIBPATH=../lapack-3.5.0/lib
LIBS=-lm -llapack -lblas -lgfortran 

ACTDIR=motifdistribution

$(OBJDIR)/%.o : %.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(PROG):$(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) -L$(LIBPATH) $(LIBS)

all:$(SRC) $(PROG)

#TMOT=x1.pwm x2.pwm x3.pwm x4.pwm x5.pwm

#T1=$(patsubst data/%.pwm,data/%.t1,$(TMOT))

#data/%.t1 : data/%.pwm
#	test1d motif=$< PVAL=0.1 ORDER=0


#alltest: $(TMOT)
#	echo $<

clean:
	rm -f $(OBJDIR)/*.o $(PROG)

