CC    =  g++
FLAG  = -O3 -Wall 

UPPE_LOC     = /Users/kolesik/Vyvoj/w04gUPPEcore
UPPE_LOC_INC = ${UPPE_LOC}/include
UPPE_LOC_LIB = ${UPPE_LOC}/lib
UPPE_LOC_SRC = ${UPPE_LOC}/src


INCL   = -I ${UPPE_LOC_INC} -I .
LIBS   = -L ${UPPE_LOC_LIB} 

###############################################################


all: clean a b c d cleanall

clean:
	rm -f *.o *~ core* *.out  *.dat TEST*.dat x*.out

cleanall: 
	rm -f *.o *~ core* *.out  *.dat TEST*.dat x*.out lib*


%.o : %.cc %.h
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.cc 
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.c %.h
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.c 
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

a:   cn-bpm-radial.cc
	${CC} ${FLAG} -o a.out cn-bpm-radial.cc ${INCL} ${LIBS} 

b:   cn-bpm-radial-B.cc
	${CC} ${FLAG} -o b.out cn-bpm-radial-B.cc ${INCL} ${LIBS} 

c:   cn-bpm-radial-C.cc
	${CC} ${FLAG} -o c.out cn-bpm-radial-C.cc ${INCL} ${LIBS} 

d:   cn-bpm-radial-D.cc
	${CC} ${FLAG} -o d.out cn-bpm-radial-D.cc ${INCL} ${LIBS} 





