IDIR =../include
CXX=g++
CXXFLAGS= -std=c++11
CPPFLAGS = -I$(IDIR) -I ../eigen-3.4.0/ -I../boost_1_79_0

ODIR=../bin

OBJ = \
MultiProc.o \
VectorMath.o \
OutputHandler.o \
InputHandler.o \
Orbital.o \
Basis.o \
Atom.o \
MatrixElement.o \
Molecule.o \
ElectronIntegrals.o \
RHF.o \
hfcalc.o

hfcalc: $(OBJ) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ 

$(OBJ): %.o : %.cpp 

%.cpp:
	$(CXX) $(CPPFLAGS) -c -o $@ $(CXXFLAGS)

.PHONY: clean check distcheck

install:
	if [ ! -d $(ODIR) ]; then mkdir $(ODIR); fi
	mv *.o $(ODIR)
	mv hfcalc $(ODIR)

clean:
	if [ ! -d $(ODIR) ]; then rm *.o hfcalc; fi
	if [ -d $(ODIR) ]; then rm -r $(ODIR); fi
