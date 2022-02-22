FLAGS := -O3
target := bin/ts4

LIBS := 
LIBS += -L${CMF}/lib -lcmf
LIBS += -L${PTL}/lib -lPTL

INCL :=
INCL += -I.
INCL += -I./src
INCL += -I${CMF}/include
INCL += -I${PTL}/include


main: setup cmf
	mpicxx -g -std=c++20 ${FLAGS} ${INCL} -c src/main.cc -o obj/main.o ${LIBS}
	mpicxx -g -std=c++20 ${FLAGS} ${INCL} obj/*.o -o ${target} ${LIBS}

cmf:
	make -C ${CMF} -f makefile

setup:
	mkdir -p bin
	mkdir -p obj

clean:
	rm -f ${target}
	rm -rf bin
	rm -rf obj
