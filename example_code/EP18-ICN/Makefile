
CC    =  g++
FLAG  = -O3 -Wall -static

UPPE_LOC     = /home/data1/Programy/UPPEcore
UPPE_LOC_INC = ${UPPE_LOC}/include
UPPE_LOC_LIB = ${UPPE_LOC}/lib
UPPE_LOC_SRC = ${UPPE_LOC}/src


INCL   = -I ${UPPE_LOC_INC} -I .
LIBS   = -L ${UPPE_LOC_LIB} -lgsl -lgslcblas  -lm

###############################################################

OBJS   =  bpm_icn_initialization.o bpm_icn_step.o bpm_observer.o 

all:    cleanall a clean

clean:
	rm -f *.o *~ core*  DBG*.dat BPM_* a.out X*

cleanall: 
	rm -f *.o *~ core* *.out  DBG*.dat TEST*.dat x*.out lib*


%.o : %.cc %.h
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.cc 
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.c %.h
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.c 
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 


a:   ${OBJS} main.cc
	${CC} ${FLAG} -o a.out ${OBJS} main.cc ${INCL} ${LIBS} 











