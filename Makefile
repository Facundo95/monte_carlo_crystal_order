CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -g
#CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -g

# Default settings for Unix/Linux
TARGET = mc_simulation
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
	# Use -$(RM) to ignore errors if files don't exist
	-$(RM) $(TARGET) $(OBJS)

doc: $(SRCS) $(DOXYFILE)
	@echo "--- Generating Doxygen documentation... ---"
	$(DOXYGEN) $(DOXYFILE)
