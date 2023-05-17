RenderX			:	Coulomb.o DEM.o ExtMath.o FF.o GR.o MWtransfer.o Neutrals.o Plasma.o RenderIrregular.o Rmain.o Zeta.o
				g++ $^ -shared -o RenderX.so
Coulomb.o		:	Coulomb.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX Coulomb.cpp
DEM.o			:	DEM.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX DEM.cpp
ExtMath.o		:	ExtMath.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX ExtMath.cpp
FF.o			:	FF.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX FF.cpp
GR.o			:	GR.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX GR.cpp
MWtransfer.o		:	MWtransfer.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX MWtransfer.cpp
Neutrals.o		:	Neutrals.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX Neutrals.cpp
Plasma.o		:	Plasma.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX Plasma.cpp
RenderIrregular.o	:	RenderIrregular.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX RenderIrregular.cpp
Rmain.o			:	Rmain.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX Rmain.cpp
Zeta.o			:	Zeta.cpp Coulomb.h DEM.h ExtMath.h FF.h GR.h IDLinterface.h MWtransfer.h Neutrals.h Plasma.h RenderIrregular.h Zeta.h
				g++ -c -std=c++11 -O3 -fPIC -D LINUX Zeta.cpp
