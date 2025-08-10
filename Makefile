OBJS	= main.o 
SOURCE	= main.cpp
HEADER	= 
OUT	= heat.out
CC	 = mpic++
FLAGS	 = -pg -c -std=c++14
LFLAGS	 = 

all: $(OBJS)
	$(CC) $(OBJS) -o $(OUT) $(LFLAGS)

main.o: main.cpp
	$(CC) $(FLAGS) main.cpp

cart: main3.o
	$(CC) main3.o -o cartesian.out $(LFLAGS)

main3.o: main3.cpp
	$(CC) $(FLAGS) main3.cpp


clean:
	rm -f $(OBJS) $(OUT)
	rm main3.o
