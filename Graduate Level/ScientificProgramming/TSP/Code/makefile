# Define the system and library format
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

# Define the paths to CPLEX, Concert, and CPOptimizer
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio2211/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio2211/concert
CPOPTIMIZERDIR = /opt/ibm/ILOG/CPLEX_Studio2211/cpoptimizer

# Define the compiler and options
CC  = g++ -O0 -g
CCOPT = -m64 -fPIC -fno-strict-aliasing -fexceptions -fpermissive -DNDEBUG -DIL_STD -std=c++17

# Define the library directories
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPOPTIMIZERLIBDIR = $(CPOPTIMIZERDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPOPTIMIZERBINDIR = $(CPOPTIMIZERDIR)/bin/$(SYSTEM)

# Define the linker directories and flags
CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) -L$(CPOPTIMIZERLIBDIR)
CLNDIRS   = -L$(CPLEXLIBDIR)
CCLNFLAGS = -lcp -lilocplex -lconcert -lcplex -lm -lpthread

# Define the include directories
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CPOPTIMIZERINCDIR = $(CPOPTIMIZERDIR)/include

# Define the compiler flags
CCFLAGS = $(CCOPT)  -I$(CONCERTINCDIR) -I$(CPLEXINCDIR) -I$(CPOPTIMIZERINCDIR)


# Define the directories for source and object files
SRCDIR = .
OBJDIR = obj

# Find all .cpp files in the source directory
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

# Define the output executable
OUTPUT = TSP_Solver

# Define the build rule
all: $(OUTPUT)

$(OUTPUT): $(OBJECTS)
	$(CC) $(CCFLAGS) $(CCLNDIRS) -o $(OUTPUT) $(OBJECTS) $(CCLNFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CC) -c $(CCFLAGS) $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

# Define the clean rule
clean:
	/bin/rm -rf $(OBJDIR) *.o *~ *.class $(OUTPUT)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp