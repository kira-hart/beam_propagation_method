
CC    =  g++
FLAG  = -O3 -Wall -static

UPPE_LOC     = /home/data1/Programy/UPPEcore
UPPE_LOC_INC = ${UPPE_LOC}/include
UPPE_LOC_LIB = ${UPPE_LOC}/lib
UPPE_LOC_SRC = ${UPPE_LOC}/src


INCL   = -I ${UPPE_LOC_INC} -I .
LIBS   = -L ${UPPE_LOC_LIB} -lgsl -lgslcblas  -lm

###############################################################

all:    cleanall a clean

clean:
	rm -f *.o *~ core*  DBG*.dat TEST*.dat a.out X*

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


a:   cmp.cc
	${CC} ${FLAG} -o a.out cmp.cc ${INCL} ${LIBS} 

b:   cmp1.cc
	${CC} ${FLAG} -o a.out cmp1.cc ${INCL} ${LIBS} 








