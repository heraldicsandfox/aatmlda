CC=		g++

OBJS=		strtokenizer.o dataset.o utils.o model.o
MAIN=		lda
 
all:	$(OBJS) $(MAIN).cpp
	$(CC) -O3 -std=c++11 -o $(MAIN) $(MAIN).cpp $(OBJS)
	strip $(MAIN)

strtokenizer.o:	strtokenizer.h strtokenizer.cpp
	$(CC) -O3 -std=c++11 -c -o strtokenizer.o strtokenizer.cpp

dataset.o:	dataset.h dataset.cpp
	$(CC) -O3 -std=c++11 -c -o dataset.o dataset.cpp

utils.o:	utils.h utils.cpp
	$(CC) -O3 -std=c++11 -c -o utils.o utils.cpp

model.o:	model.h model.cpp
	$(CC) -O3 -std=c++11 -c -o model.o model.cpp

test:
	

clean:
	rm $(OBJS) 
	rm $(MAIN)

