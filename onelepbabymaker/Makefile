ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#######################################
# -- DEFINE ARCH to something sensible!
#######################################
 
#
CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -pg -Wextra -pedantic -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel  -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option 
#EXTRACFLAGS=-Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel  -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option 
LD            = g++
LDFLAGS       = -g -pg
SOFLAGS       = -shared


CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

COREPATH       = ../CORE
TOOLSPATH      = $(COREPATH)/Tools
DORKYPATH      = $(TOOLSPATH)/dorky


CXXFLAGS       += -I$(COREPATH)
CXXFLAGS       += -I$(TOOLSPATH)
CXXFLAGS       += -I$(DORKYPATH)

CORELIB        = $(COREPATH)/CMS3_CORE.so
STOPVARLIB     = ../analysisutils/libAnalysisUtils.so

#runBabyMaker: runBabyMaker.o looper.o EventTree.o LeptonTree.o JetTree.o IsoTracksTree.o TauTree.o GenParticleTree.o PhotonTree.o StopSelections.o LinkDef_out.so
#	$(LD) $(LDFLAGS)  -o runBabyMaker runBabyMaker.o looper.o EventTree.o JetTree.o LeptonTree.o PhotonTree.o IsoTracksTree.o TauTree.o GenParticleTree.o StopSelections.o $(GLIBS) -lGenVector -lEG $(CORELIB) $(STOPVARLIB) LinkDef_out.so -Wl,-rpath,./

runBabyMakerCondor: runBabyMakerCondor.o looper.o EventTree.o LeptonTree.o JetTree.o IsoTracksTree.o TauTree.o GenParticleTree.o PhotonTree.o StopSelections.o LinkDef_out.so
	$(LD) $(LDFLAGS)  -o runBabyMakerCondor runBabyMakerCondor.o looper.o EventTree.o JetTree.o LeptonTree.o PhotonTree.o IsoTracksTree.o TauTree.o GenParticleTree.o StopSelections.o $(GLIBS) -lGenVector -lEG $(CORELIB) $(STOPVARLIB) LinkDef_out.so -Wl,-rpath,./
#runBabyMaker.o: runBabyMaker.C looper.h 
#	$(CXX) $(CXXFLAGS) -c runBabyMaker.C 
runBabyMakerCondor.o: runBabyMakerCondor.C looper.h 
	$(CXX) $(CXXFLAGS) -c runBabyMakerCondor.C 

looper.o: looper.C looper.h EventTree.h TauTree.h IsoTracksTree.h GenParticleTree.h LeptonTree.h PhotonTree.h StopSelections.h
	$(CXX) $(CXXFLAGS) -c looper.C

EventTree.o: EventTree.cc EventTree.h StopSelections.h 
	$(CXX) $(CXXFLAGS) -c EventTree.cc 

LeptonTree.o: LeptonTree.cc LeptonTree.h StopSelections.h
	$(CXX) $(CXXFLAGS) -c LeptonTree.cc 

PhotonTree.o: PhotonTree.cc PhotonTree.h StopSelections.h
	$(CXX) $(CXXFLAGS) -c PhotonTree.cc 

JetTree.o: JetTree.cc JetTree.h StopSelections.h
	$(CXX) $(CXXFLAGS) -c JetTree.cc

TauTree.o: TauTree.cc TauTree.h
	$(CXX) $(CXXFLAGS) -c TauTree.cc

IsoTracksTree.o: IsoTracksTree.cc IsoTracksTree.h
	$(CXX) $(CXXFLAGS) -c IsoTracksTree.cc

GenParticleTree.o: GenParticleTree.cc GenParticleTree.h
	$(CXX) $(CXXFLAGS) -c GenParticleTree.cc

StopSelections.o: StopSelections.cc StopSelections.h
	$(CXX) $(CXXFLAGS) -c StopSelections.cc

LinkDef_out.so: LinkDef.h
	@rootcling -v -f LinkDef.cc -c -p Math/Vector4D.h LinkDef.h
	$(CXX) $(CXXFLAGS) $(SOFLAGS) LinkDef.cc -o $@

.PHONY: clean
clean:  
	rm -v -f \
	runBabyMaker \
	LinkDef.cc \
	LinkDef_rdict.pcm \
	*.o *.d *.so; echo "Done"
