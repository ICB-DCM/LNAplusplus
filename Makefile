C_SRC=$(wildcard models/BirthDeath/C/*.c)
C_OBJS=$(C_SRC:.c=.o)
CPP=g++
CC=gcc

INCLUDES=-Imodels/BirthDeath/C -Iinclude

all: BirthDeath

debug: CFLAGS+=-g
debug: CPPFLAGS+=-g
debug: BirthDeath

BirthDeath: $(C_OBJS) src/computeLinearNoise.o src/main.cpp
	$(CPP) $(CPPFLAGS) -o BirthDeath src/main.cpp $(C_OBJS) src/computeLinearNoise.o $(INCLUDES) -lsundials_cvodes -lsundials_nvecserial -lblitz -lstdc++ -lgsl -lgslcblas

src/computeLinearNoise.o: src/computeLinearNoise.cpp
	$(CPP) $(CPPFLAGS) -c -o src/computeLinearNoise.o src/computeLinearNoise.cpp $(INCLUDES)

models/BirthDeath/C/%.o:models/BirthDeath/C/%.c
	$(CC) -c $(CFLAGS) -o $@ $< -Imodels/BirthDeath/C

.PHONY: clean

clean:
	rm src/*.o models/BirthDeath/C/*.o

