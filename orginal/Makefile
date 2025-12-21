CXX = clang++
CXXFLAGS = -std=c++17 -Wall -g -O0
SRCDIR = src
INCDIR = inc
OBJDIR = obj

SOURCES = main.cpp \
          $(SRCDIR)/graph.cpp \
          $(SRCDIR)/common.cpp \
          $(SRCDIR)/helpers.cpp

OBJECTS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(SOURCES))

TARGET = bk_algorithm

# Default target
all: $(TARGET)

# Create obj directories
$(OBJDIR):
	mkdir -p $(OBJDIR)/$(SRCDIR)

# Link
$(TARGET): $(OBJDIR) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS)

# Compile
$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

# Clean
clean:
	rm -rf $(OBJDIR) $(TARGET)

.PHONY: all clean
