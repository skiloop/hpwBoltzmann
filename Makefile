CXX = g++
CXXFLAGS = -g -D_DEBUG -Wall #-DMATLAB_SIMULATION 
#CXXFLAGS = -O3 -Wall #-DMATLAB_SIMULATION 
LFLAGS = -g -O3
SRC=src
MATLAB_ROOT=/opt/MATLAB/R2011a
MATLAB_LIB=$(MATLAB_ROOT)/bin/glnxa64
MATLAB_INC=$(MATLAB_ROOT)/extern/include
#INC= -I$(MATLAB_INC)
#LIB= -L$(MATLAB_LIB) -leng -lmx -Wl,-rpath-link $(MATLAB_LIB)
COBJS = datastruct.o commondata.o
OBJS = main.o $(COBJS) bounddata.o\
connectingdef.o matlabsim.o \
initials.o\
pmlbound.o\
updatedensity.o\
updatefields.o\
fdtdloop.o FileInterp.o interp.o


.PHONY: all clean

all:hpw
hpw: $(OBJS) $(SRC)/*.cpp
	$(CXX) $(LFLAGS) -o $@  $(OBJS) $(LIB)
$(OBJS):%.o:$(SRC)/%.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INC)
clean:
	-rm -f *.o hpw
