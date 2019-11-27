.KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -Wall -D NDEBUG -O3 -m64
LIB		= -lz -lm -lpthread -lbz2 -llzma -lstdc++
SOURCE		= main.cpp GetData.cpp Mapping.cpp AlignmentCandidates.cpp AlignmentRescue.cpp KmerAnalysis.cpp nw_alignment.cpp tools.cpp bwt_index.cpp bwt_search.cpp
HEADER		= structure.h
OBJECT		= $(SOURCE:%.cpp=%.o)

%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<

main:		$(OBJECT)
			$(CXX) $(FLAGS) $(OBJECT) -o kart BWT_Index/libbwa.a htslib/libhts.a $(LIB)
			
clean:
		rm -f *.o *~
