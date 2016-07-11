.KEEP_STAT:

all: main

Compiler        = gcc
FLAGS           = -g -O2 -Wall -lz -lm
LIB             =
SOURCE          = wgsim.c
OBJECT          = $(SOURCE:%.cpp=%.o)

main:           $(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -lpthread -o wgsim $(LIB)
