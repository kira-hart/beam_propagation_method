CC    =  g++
FLAG  = -O3 -Wall

UPPE_LOC     = ~/Vyvoj/w04gUPPEcore
UPPE_LOC_INC = ${UPPE_LOC}/include
UPPE_LOC_LIB = ${UPPE_LOC}/lib
UPPE_LOC_SRC = ${UPPE_LOC}/src


INCL   = -I ${UPPE_LOC_INC} -I .
LIBS   = -L ${UPPE_LOC_LIB} 

###############################################################


all: clean a b clean

clean:
	rm -f *.o *~ core*   TEST*.dat x*.out

cleanall: 
	rm -f *.o *~ core* *.out  TEST*.dat x*.out lib*


%.o : %.cc %.h
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.cc 
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.c %.h
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.c 
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

a:   cn-bpm-radial-Kerr.cc
	${CC} ${FLAG} -o a.out cn-bpm-radial-Kerr.cc ${INCL} ${LIBS} 

b:   show_sf_distance.cc
	g++ -o b.out show_sf_distance.cc

