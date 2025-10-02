CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -g
TARGET = mc_simulation
SRCS = main.cpp simulation.cpp rng.cpp file_handler.cpp
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
