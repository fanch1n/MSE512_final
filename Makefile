.SUFFIXES: .cpp .o 
.cpp.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INC) -c $<

.cc.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INC) -c $<

# ----- Make Macros -----
OPTFLAGS = -O3
CXX	=	g++

INC = -I/usr/include
LIB = -lgsl -lgslcblas 

CXXFLAGS = -g -std=c++11 $(OPTFLAGS) 

CUR_DIR = `pwd`
TARGETS = run
OBJECTS =	run.o 
# ----- Make Rules -----
all:	$(TARGETS) 

run: $(OBJECTS)  	
	$(CXX) $(OBJECTS) $(CXXFLAGS) -o run $(INC) $(LIB)

clean:
	rm -f $(TARGETS) $(OBJECTS) *.o

