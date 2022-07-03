IDIR =../include
CXX=g++
CXXFLAGS=-I$(IDIR) -I ../eigen-3.4.0/ -I../boost_1_79_0

ODIR=../bin

OBJ = \
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
	$(CXX) -o $@ $^ $(CFLAGS)

$(OBJ): %.o : %.cpp

%.cpp:
	$(CXX) $(CXXFLAGS) -c -o $@ 

.PHONY: clean check distcheck

install:
	if [ ! -d $(ODIR) ]; then mkdir $(ODIR); fi
	mv *.o $(ODIR)
	mv hfcalc $(ODIR)

clean:
	if [ -d $(ODIR) ]; then rm -r $(ODIR); fi
	if [ ! -d $(ODIR) ]; then rm *.o hfcalc; fi
