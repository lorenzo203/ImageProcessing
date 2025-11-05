# Compiler
CXX      := g++
CXXFLAGS := -O2 -std=c++17 -fopenmp
mkEigenInc ?= /usr/include/eigen3
EIGEN_INC ?= ${mkEigenInc}

# LIS
LIS_DIR  := lis-2.1.10
LIS_LIB  := $(LIS_DIR)/src/.libs/liblis.a

INCLUDES := -I$(EIGEN_INC) -I$(LIS_DIR)/include -Isrc
LDFLAGS  := -L$(LIS_DIR)/src/.libs
LDLIBS   := -llis -lm

# Program
TARGET   := challenge1
SRC      := src/challenge1.cpp

.PHONY: all run clean
all: $(TARGET)


$(TARGET): $(SRC) $(LIS_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< $(LDFLAGS) $(LDLIBS) -o $@

$(LIS_LIB):
	@echo "==> Building LIS (only if missing)"
	@if [ ! -f "$(LIS_DIR)/config.status" ]; then \
        cd $(LIS_DIR) && CFLAGS="-fPIC" ./configure; \
    fi
	$(MAKE) -C $(LIS_DIR)

run: $(TARGET)
	@./$(TARGET)

clean:
	@rm -f $(TARGET) *.o
	@echo "Cleaned"