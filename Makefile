# Makefile

# Compiler
CXX      := g++
CXXFLAGS := -O2 -std=c++17 -fopenmp
INCLUDES := -I/usr/include/eigen3 -I./lis-2.1.10/include
LDFLAGS  := -L./lis-2.1.10/src/.libs
LDLIBS   := -llis -lm

# Targets
TARGET   := challenge1
SRC      := challenge1.cpp

# Default rule: build lis first, then your program
all: lis $(TARGET)

# Rule to build LIS (in-place, no install directory)
lis:
	cd lis-2.1.10 && make clean
	cd lis-2.1.10 && CFLAGS="-fPIC" ./configure
	cd lis-2.1.10 && make

# Rule to build your program
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ $(LDFLAGS) $(LDLIBS) -o $@

# Cleanup
clean:
	rm -f $(TARGET) *.o
	cd lis-2.1.10 && make clean

