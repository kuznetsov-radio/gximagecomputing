program_base := RenderGRFF
program_name := $(program_base).so

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

CXX ?= g++
CPPFLAGS ?=
CXXFLAGS ?= -std=c++11 -O3 -fPIC
LDFLAGS ?=
LDLIBS ?=

DEFS :=
OMP_CXXFLAGS :=
OMP_LDLIBS :=
ARCH_SUFFIX :=

ifeq ($(UNAME_S),Darwin)
    DEFS += -DMACOS
    OMP_CXXFLAGS += -Xpreprocessor -fopenmp
    OMP_LDLIBS += -lomp
    ARCH_SUFFIX := _$(UNAME_M)

    # Prefer Homebrew libomp if available.
    HOMEBREW_PREFIX ?= $(shell brew --prefix 2>/dev/null)
    ifneq ($(HOMEBREW_PREFIX),)
        CPPFLAGS += -I$(HOMEBREW_PREFIX)/include -I$(HOMEBREW_PREFIX)/opt/libomp/include
        LDFLAGS += -L$(HOMEBREW_PREFIX)/lib -L$(HOMEBREW_PREFIX)/opt/libomp/lib
    endif
else ifeq ($(UNAME_S),Linux)
    DEFS += -DLINUX
    OMP_CXXFLAGS += -fopenmp
    OMP_LDLIBS += -fopenmp
else
    $(warning Unsupported platform $(UNAME_S). Build may fail.)
endif

OBJ := \
    Coulomb.o \
    DEM.o \
    ExtMath.o \
    FF.o \
    GR.o \
    MWtransfer.o \
    Neutrals.o \
    Plasma.o \
    RenderIrregular.o \
    Rmain.o \
    Zeta.o \
    EUVmain.o \
    GXexports.o \
    PyInterface.o

.PHONY: all install-binary clean

all: $(program_name) install-binary

$(program_name): $(OBJ)
	$(CXX) -shared -o $@ $^ $(LDFLAGS) $(OMP_LDLIBS) $(LDLIBS)

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(DEFS) $(OMP_CXXFLAGS) -o $@ $<

install-binary: $(program_name)
	mkdir -p ../binaries
ifeq ($(UNAME_S),Darwin)
	cp $(program_name) ../binaries/$(program_base)$(ARCH_SUFFIX).so
else
	cp $(program_name) ../binaries/$(program_base).so
endif

clean:
	rm -f $(program_name) $(OBJ)
