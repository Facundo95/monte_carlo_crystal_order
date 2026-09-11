CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -g
#CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -g

# Default settings for Unix/Linux
TARGET = mc_simulation
TEST_TARGET = test_heisenberg
RM = rm -f

# Check for Windows environments (cmd.exe or PowerShell)
ifeq ($(OS),Windows_NT)
    TARGET = mc_simulation.exe
    RM = del /Q
# Check for Unix-like environments on Windows (MinGW, Cygwin)
else
    UNAME_S := $(shell uname -s)
    ifeq ($(findstring CYGWIN,$(UNAME_S)),CYGWIN)
        TARGET = mc_simulation.exe
    endif
    ifeq ($(findstring MINGW,$(UNAME_S)),MINGW)
        TARGET = mc_simulation.exe
    endif
endif

SRCDIR = ../src
SRCS = $(SRCDIR)/main.cpp $(SRCDIR)/simulation.cpp $(SRCDIR)/rng.cpp $(SRCDIR)/file_handler.cpp $(SRCDIR)/lattice.cpp \
    $(SRCDIR)/lattice_output.cpp $(SRCDIR)/input_parser.cpp $(SRCDIR)/mc_helpers.cpp \
    $(SRCDIR)/neighbor_sums.cpp $(SRCDIR)/sweep_list.cpp $(SRCDIR)/beg_hamiltonian.cpp \
    $(SRCDIR)/ising_hamiltonian.cpp $(SRCDIR)/heisenberg_hamiltonian.cpp
    SRCS += $(SRCDIR)/logger.cpp
OBJS = $(SRCS:.cpp=.o)
TEST_SRCS = ../tests/test_heisenberg.cpp $(filter-out $(SRCDIR)/main.cpp,$(SRCS))
TEST_OBJS = $(TEST_SRCS:.cpp=.o)
DOXYGEN = doxygen
DOXYFILE = Doxyfile

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

test-heisenberg: $(TEST_TARGET)
	./$(TEST_TARGET)

$(TEST_TARGET): $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(SRCDIR) -c -o $@ $<

clean:
	-$(RM) $(TARGET) $(TEST_TARGET) $(OBJS) $(TEST_OBJS)

doc: $(SRCS) $(DOXYFILE)
	@echo "--- Generating Doxygen documentation... ---"
	$(DOXYGEN) $(DOXYFILE)
