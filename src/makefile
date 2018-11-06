.KEEP_STAT:

all: main

Compiler	= g++
FLAGS		= -D NDEBUG -O3 -m64
LIB		= -lz -lm -lpthread -lstdc++
SOURCE		= main.cpp GetData.cpp Mapping.cpp AlignmentCandidates.cpp AlignmentRescue.cpp KmerAnalysis.cpp nw_alignment.cpp tools.cpp bwt_index.cpp bwt_search.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o kart $(LIB)
			
static:		$(OBJECT)
			$(Compiler) $(FLAGS) -static $(OBJECT) -o kart_static $(LIB)

%.o:		%.cpp $(HEADER)
			$(Compiler) $(FLAGS) -c $<

clean:
		rm -f *.o *~
