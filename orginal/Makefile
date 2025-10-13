CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2
SRCDIR = src
INCDIR = inc
OBJDIR = obj

# Source files
SOURCES = main.cpp $(SRCDIR)/graph.cpp $(SRCDIR)/common.cpp $(SRCDIR)/helpers.cpp
OBJECTS = $(SOURCES:.cpp=.o)

TARGET = bk_algorithm

# Default target
all: $(TARGET)

# Create object directory
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Build target
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

# Clean build files
clean:
	rm -f $(OBJECTS) $(TARGET)

# Phony targets
.PHONY: all clean