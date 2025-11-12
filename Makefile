# Compiler
CXX      := g++
CXXFLAGS := -O2 -std=c++17 -fopenmp
mkEigenInc ?= /usr/include/eigen3
EIGEN_INC ?= ${mkEigenInc}

# LIS
LIS_DIR  := lis-2.1.10
LIS_ZIP  := $(LIS_DIR).zip
LIS_URL  := https://www.ssisc.org/lis/dl/lis-2.1.10.zip
LIS_LIB  := $(LIS_DIR)/src/.libs/liblis.a

INCLUDES := -I$(EIGEN_INC) -I$(LIS_DIR)/include -Isrc
LDFLAGS  := -L$(LIS_DIR)/src/.libs
LDLIBS   := -llis -lm

# Program
TARGET   := challenge1
SRC      := src/challenge1.cpp

.PHONY: all run clean distclean
all: $(TARGET)


$(TARGET): $(SRC) $(LIS_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< $(LDFLAGS) $(LDLIBS) -o $@

$(LIS_DIR):
	@if [ ! -f $(LIS_ZIP) ]; then \
		echo "==> Downloading LIS $(LIS_DIR) from official site..."; \
		curl -L -f -o $(LIS_ZIP) $(LIS_URL) || \
		wget --no-check-certificate -O $(LIS_ZIP) $(LIS_URL) || \
		(rm -f $(LIS_ZIP); \
		 echo "ERROR: Could not download LIS. Please download manually from $(LIS_URL)"; \
		 exit 1); \
		echo "==> Download successful"; \
	else \
		echo "==> Found existing $(LIS_ZIP)"; \
	fi
	@echo "==> Extracting LIS..."
	@unzip -q $(LIS_ZIP) || \
	(echo "ERROR: Failed to extract $(LIS_ZIP). File may be corrupted."; \
	 echo "Run 'make distclean' and try again."; \
	 exit 1)
	@echo "==> Extraction complete"
	@echo "==> Using pre-built LIS library"

$(LIS_LIB): $(LIS_DIR)
	@if [ ! -f "$(LIS_LIB)" ]; then \
		echo "==> Configuring and building LIS (this may take a few minutes)..."; \
		cd $(LIS_DIR) && \
		if [ ! -f config.status ]; then \
			CFLAGS="-fPIC" ./configure --quiet > /dev/null 2>&1; \
		fi && \
		$(MAKE) -C src all > /dev/null 2>&1; \
		echo "==> LIS library built successfully"; \
	else \
		echo "==> LIS library already built"; \
	fi

run: $(TARGET)
	@./$(TARGET)

clean:
	@rm -f $(TARGET) *.o
	@echo "Cleaned"

distclean: clean
	@rm -rf $(LIS_DIR) $(LIS_ZIP)
	@echo "Deep cleaned (removed LIS)"