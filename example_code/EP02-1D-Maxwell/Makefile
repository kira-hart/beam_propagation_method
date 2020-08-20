
CC    =  gcc
FLAG  = -O3 -Wall 

###############################################################


all: clean a cleanall

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


a:   OneDMaxwell-WithBoundary.c
	${CC} ${FLAG} -o a.out OneDMaxwell-WithBoundary.c ${INCL} ${LIBS}










