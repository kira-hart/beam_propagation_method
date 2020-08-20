
CC    =  gcc
FLAG  = -O3 -Wall 

UPPE_LOC     = /Users/kolesik/Vyvoj/w04gUPPEcore
#UPPE_LOC     = /home/data1/Programy/UPPEcore
UPPE_LOC_INC = ${UPPE_LOC}/include
UPPE_LOC_LIB = ${UPPE_LOC}/lib
UPPE_LOC_SRC = ${UPPE_LOC}/src


INCL   = -I ${UPPE_LOC_INC} -I .
LIBS   = -L ${UPPE_LOC_LIB} -lfftw3  -lm

###############################################################


all: clean uppe cleanall

clean:
	rm -f *.o *~ core* *.out  DBG*.dat TEST*.dat x*.out

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


a:   cmp.c
	${CC} ${FLAG} -o a.out cmp.c ${INCL} ${LIBS} 

b:   cmpomp.c
	${CC} ${FLAG} -fopenmp -o b.out cmpomp.c ${INCL} ${LIBS} 

c:   Solution/instructors.c
	${CC} ${FLAG} -fopenmp -o c.out Solution/instructors.c ${INCL} ${LIBS} 







