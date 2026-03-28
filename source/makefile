program_base := RenderGRFF
program_name := $(program_base).so

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

CXX ?= g++
CPPFLAGS ?=
CXXFLAGS ?=
LDFLAGS ?=
LDLIBS ?=

PROJECT_CPPFLAGS :=
PROJECT_CXXFLAGS := -std=c++11 -O3 -fPIC
PROJECT_LDFLAGS :=
PROJECT_LDLIBS :=
LINK_OBJS = $(OBJ)

DEFS :=
OMP_CXXFLAGS :=
OMP_LDLIBS :=
ARCH_SUFFIX :=
LINK_SHARED_FLAGS := -shared
PYTHON_BIN ?= $(shell command -v python 2>/dev/null || command -v python3 2>/dev/null)
PYTHON_PREFIX := $(strip $(shell if [ -n "$(PYTHON_BIN)" ]; then "$(PYTHON_BIN)" -c 'import sys; print(sys.prefix)' 2>/dev/null; fi))
DARWIN_LINK_PREFIX_FLAGS :=
DARWIN_LINK_SUFFIX_FLAGS :=

ifeq ($(UNAME_S),Darwin)
    DEFS += -DMACOS
    OMP_CXXFLAGS += -Xpreprocessor -fopenmp
    OMP_LDLIBS += -lomp
    ARCH_SUFFIX := _$(UNAME_M)
    LINK_SHARED_FLAGS := -bundle -undefined dynamic_lookup
    LINK_OBJS = pyinit.o $(OBJ)

    # Prefer Homebrew libomp if available.
    HOMEBREW_PREFIX ?= $(shell brew --prefix 2>/dev/null)
    ifneq ($(HOMEBREW_PREFIX),)
        PROJECT_CPPFLAGS += -I$(HOMEBREW_PREFIX)/include -I$(HOMEBREW_PREFIX)/opt/libomp/include
        PROJECT_LDFLAGS += -L$(HOMEBREW_PREFIX)/lib -L$(HOMEBREW_PREFIX)/opt/libomp/lib
    endif
    ifneq ($(PYTHON_PREFIX),)
        DARWIN_LINK_PREFIX_FLAGS += -Wl,-rpath,$(PYTHON_PREFIX)/lib -L$(PYTHON_PREFIX)/lib
        DARWIN_LINK_PREFIX_FLAGS += -Wl,-rpath,$(PYTHON_PREFIX)/lib -L$(PYTHON_PREFIX)/lib
    endif
    ifneq ($(HOMEBREW_PREFIX),)
        DARWIN_LINK_PREFIX_FLAGS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
        DARWIN_LINK_SUFFIX_FLAGS += -L$(HOMEBREW_PREFIX)/lib -L$(HOMEBREW_PREFIX)/opt/libomp/lib
    endif
    DARWIN_LINK_SUFFIX_FLAGS += -lomp
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
    EUVmain.o \
    ExtMath.o \
    FF.o \
    GR.o \
    GXexports.o \
    MWtransfer.o \
    Messages.o \
    Neutrals.o \
    Plasma.o \
    PyInterface.o \
    RenderIrregular.o \
    Rmain.o \
    Zeta.o

.PHONY: all install-binary clean

all: $(program_name) install-binary

ifeq ($(UNAME_S),Darwin)
$(program_name): $(LINK_OBJS)
	$(CXX) $(LINK_SHARED_FLAGS) $(DARWIN_LINK_PREFIX_FLAGS) $(LINK_OBJS) -o $@ $(DARWIN_LINK_SUFFIX_FLAGS)
else
$(program_name): $(LINK_OBJS)
	$(CXX) $(LINK_SHARED_FLAGS) -o $@ $(LINK_OBJS) $(LDFLAGS) $(PROJECT_LDFLAGS) $(OMP_LDLIBS) $(LDLIBS)
endif

pyinit.o: ../pyinit.cpp
	$(CXX) -c $(CPPFLAGS) $(PROJECT_CPPFLAGS) $(CXXFLAGS) $(PROJECT_CXXFLAGS) -o $@ $<

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(PROJECT_CPPFLAGS) $(CXXFLAGS) $(PROJECT_CXXFLAGS) $(DEFS) $(OMP_CXXFLAGS) -o $@ $<

install-binary: $(program_name)
	mkdir -p ../binaries
ifeq ($(UNAME_S),Darwin)
	cp $(program_name) ../binaries/$(program_base)$(ARCH_SUFFIX).so
else
	cp $(program_name) ../binaries/$(program_base).so
endif

clean:
	rm -f $(program_name) $(OBJ) pyinit.o
