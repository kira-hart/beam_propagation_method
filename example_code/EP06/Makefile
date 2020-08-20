
CC    =  g++
MPICC = /usr/local/bin/mpic++
FLAG  = -O3  -Wno-deprecated

#WRKDIR_LOC     = /Users/kolesik/Programy/ngUPPE
WRKDIR_LOC     = /Users/kolesik/Vyvoj/w04gUPPEcore
WRKDIR_LOC_INC = ${WRKDIR_LOC}/include
WRKDIR_LOC_LIB = ${WRKDIR_LOC}/lib
WRKDIR_LOC_SRC = ${WRKDIR_LOC}/src

INCL   = -I ${WRKDIR_LOC_INC} -I .
LIBS   =  -L ${WRKDIR_LOC_LIB}    -lguppecore -lfftw3 -lgsl -lgslcblas   

###############################################################


OBJS = 


all: cleanall guppe t clean

clean:
	rm -f *.o *~ core*  DBG*.dat TEST*.dat x*.out

cleanall: 
	rm -f *.o *~ core* *.out  DBG*.dat TEST*.dat x*.out lib* RRR*


%.o : %.cc %.hg
	${MPICC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.cc 
	${MPICC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.c
	${MPICC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 


objs : ${OBJS}


t:
	${CC} ${FLAG}  -o t.out myDHT.cc ${INCL} ${LIBS}

a:
	${CC} ${FLAG}  -o a.out DHT-BPM.cc ${INCL} ${LIBS}



