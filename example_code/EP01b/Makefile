
CC    =  gcc
FLAG  = -O3 -Wall 

UPPE_LOC     = /Users/kolesik/Vyvoj/w04gUPPEcore
UPPE_LOC_INC = ${UPPE_LOC}/include
UPPE_LOC_LIB = ${UPPE_LOC}/lib
UPPE_LOC_SRC = ${UPPE_LOC}/src


INCL   = -I ${UPPE_LOC_INC} -I .
LIBS   = -L ${UPPE_LOC_LIB} -lfftw3  -lm

###############################################################


all: clean g b c clean

clean:
	rm -f *.o *~ core* *.out  DBG*.dat TEST*.dat x*.out

cleanall: 
	rm -f *.o *~ core* *.out  *.dat *.log *.aux TEST*.dat  lib*


%.o : %.cc %.h
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.cc 
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.c %.h
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.c 
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

b:   YeeDispersion_1.c
	${CC} ${FLAG} -o b.out YeeDispersion_1.c ${INCL} ${LIBS} 

c:   YeeDispersion_2.c
	${CC} ${FLAG} -o c.out YeeDispersion_2.c ${INCL} ${LIBS} 








