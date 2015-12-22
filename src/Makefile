### BayeScan makefile

# set variables
CXX = g++
CXX_FLAGS = -lpthread

UNAME := $(shell uname -s)
ifeq ($(UNAME), Darwin)
	CXX = g++-4.8 
else
	CXX_FLAGS += -static
endif

# build commands
bayescan_2.1: start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o 
	$(CXX) -fopenmp $(CXX_FLAGS) -o bayescan_2.1 start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o
	cp bayescan_2.1 ../inst/bin/bayescan_2.1

start.o: start.cpp errors.cpp anyoption.h global_defs.h
	$(CXX) -fopenmp -c start.cpp errors.cpp 

beta.o: beta.cpp global_defs.h
	$(CXX) -fopenmp -c beta.cpp 
      
dirichlet.o: dirichlet.cpp global_defs.h
	$(CXX) -fopenmp -c dirichlet.cpp 

RJupdates.o: RJupdates.cpp global_defs.h
	$(CXX) -fopenmp -c RJupdates.cpp 

MHupdates.o: MHupdates.cpp global_defs.h
	$(CXX) -fopenmp -c MHupdates.cpp 

likelihood.o: likelihood.cpp global_defs.h
	$(CXX) -fopenmp -c likelihood.cpp 

read_write.o: read_write.cpp errors.cpp global_defs.h
	$(CXX) -fopenmp -c read_write.cpp errors.cpp 

anyoption.o: anyoption.cpp anyoption.h 
	$(CXX) -fopenmp -c anyoption.cpp 

clean: 
	rm *.o bayescan_2.1
