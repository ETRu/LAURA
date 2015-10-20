
################ template makefile ##############


BIN	= Project
OBJ	= draw.o form.o Frame.o Main.o igraphutil.o logdebug.o randomgen.o dinamica.o networkop.o dataanalysis.o

# We don't know what compiler to use to build fltk on this machine - but fltk-config does...
CC	= $(shell fltk-config --cc)
CXX	= $(shell fltk-config --cxx)

# Set the flags for compiler: fltk-config knows the basic settings, then we can add our own...
CFLAGS		= $(shell fltk-config --cflags) -O3 
CXXFLAGS	= $(shell fltk-config --cxxflags) -w -O3 

# We don't know what libraries to link with: fltk-config does...
LINKFLTK = $(shell fltk-config --ldstaticflags)
LINKFLTK_GL = $(shell fltk-config --use-gl --ldstaticflags)
LINKFLTK_IMG = $(shell fltk-config --use-images --ldstaticflags)

LINK_IGRAPH = -I/usr/local/include/igraph -L/usr/local/lib -ligraph

# Possible steps to run after linking...
STRIP		= strip
POSTBUILD	= fltk-config --post # Required on OSX, does nothing on other platforms, so safe to call


# Define what your target application is called
all: $(BIN)

# Define how to build the various object files...

draw.o:  draw.cpp draw.h
	$(CXX) -c $< $(CXXFLAGS) $(LINK_IGRAPH)

form.o:  form.cpp form.h
	$(CXX) -c $< $(CXXFLAGS) $(LINK_IGRAPH)


Frame.o:  Frame.cpp Frame.h
	$(CXX) -c $< $(CXXFLAGS) $(LINK_IGRAPH)

igraphutil.o:  igraphutil.cpp igraphutil.h
	$(CXX) -c $< $(CXXFLAGS) $(LINK_IGRAPH)

logdebug.o:  logdebug.cpp logdebug.h
	$(CXX) -c $< $(CXXFLAGS) 

randomgen.o:  randomgen.cpp randomgen.h
	$(CXX) -c $< $(CXXFLAGS) 

dinamica.o:  dinamica.cpp dinamica.h
	$(CXX) -c $< $(CXXFLAGS) $(LINK_IGRAPH)

networkop.o:  networkop.cpp networkop.h
	$(CXX) -c $< $(CXXFLAGS) $(LINK_IGRAPH)

dataanalysis.o:  dataanalysis.cpp dataanalysis.h
	$(CXX) -c $< $(CXXFLAGS) $(LINK_IGRAPH)


Main.o:  Main.cpp 
	$(CXX) -c $< $(CXXFLAGS) $(LINK_IGRAPH)




# Now define how to link the final app - let's assume it needs image and OpenGL support
$(BIN): $(OBJ)
	$(CXX) -o $@ $(OBJ)  $(LINKFLTK_IMG) $(LINKFLTK_GL) $(LINK_IGRAPH)
	$(STRIP) $@
	$(POSTBUILD) $@  # only required on OSX, but call it anyway for portability


############### end #################