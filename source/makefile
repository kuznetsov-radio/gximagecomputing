program_name   := RenderGRFF.so
COMPILER     := g++
compile_flags  := -std=c++11 -O3 -fPIC -fopenmp -D LINUX

.PHONY: all clean

all: $(program_name)

$(program_name): Coulomb.o DEM.o ExtMath.o FF.o GR.o MWtransfer.o Neutrals.o Plasma.o RenderIrregular.o Rmain.o Zeta.o EUVmain.o GXexports.o PyInterface.o
	$(COMPILER) -shared -fopenmp -o $(program_name) $^

%.o: %.cpp
	$(COMPILER) -c $(compile_flags) -o $@ $<

clean:
	rm -f $(program_name) *.o