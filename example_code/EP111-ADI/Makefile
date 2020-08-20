CC    =  g++
FLAG  = -O3  -Wno-deprecated

INCL= -I .
LIBS=

###############################################################

OBJS = tdm_utils.o observer.o adi_solver.o

all: clean adi cleanall

clean:
	rm -f *.o *~ core* *.out  

cleanall: 
	rm -f *.o *~ core*  map_*.dat sol_*.dat 


%.o : %.cc %.h
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

%.o : %.cc 
	${CC} ${FLAG} -c $< -o $@ ${ILOC} ${INCL} 

adi:    test_ADI_BPM.cc ${OBJS}
	${CC} ${FLAG} -o adi.out ${OBJS} test_ADI_BPM.cc  ${INCL} ${LIBS}

d:      test_ADI_BPM.cc ${OBJS}
	${CC} ${FLAG} -o ndi.out ${OBJS} dispersion_ADI_BPM.cc  ${INCL} ${LIBS}








