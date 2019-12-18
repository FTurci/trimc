################################################################################
### CUSTOMIZE BELOW HERE #######################################################

CXX=g++-9 # define the compiler to use
TARGET=main # define the name of the executable
SOURCES=main.cpp box.cpp particle.cpp system.cpp histogram.cpp interactions.cpp cells.cpp
CXXFLAGS=-Ofast   -std=c++11 -I/Users/francesco/mymc/eigen/ 
LFLAGS=-lm 

################################################################################
### DO NOT EDIT THE FOLLOWING LINES ############################################

# define list of objects
OBJSC=$(SOURCES:.c=.o)
OBJS=$(OBJSC:.cpp=.o)


# DEBUG ?= 1
# ifeq ($(DEBUG), 1)
#     CXXFLAGS += -DDEBUG
# else
#     CXXFLAGS += -DNDEBUG
# endif




%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -MMD $<

# the target is obtained linking all .o files
all: $(SOURCES) $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LFLAGS)  $(OBJS) -o $(TARGET)

purge: clean
	rm -f $(TARGET)

clean:
	rm -f *.o


-include *.d
################################################################################
################################################################################
