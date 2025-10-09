CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -g
#CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -g
TARGET = mc_simulation
SRCS = main.cpp simulation.cpp rng.cpp file_handler.cpp lattice.cpp input_parser.cpp
OBJS = $(SRCS:.cpp=.o)
DOXYGEN = doxygen
DOXYFILE = Doxyfile

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGET) $(OBJS)

doc: $(SRCS) $(DOXYFILE)
	@echo "--- Generating Doxygen documentation... ---"
	$(DOXYGEN) $(DOXYFILE)
