CC = g++
PROG = topThrow

SRCS = main.cpp

LIBS = -w -lglut -lGL -lGLEW -lm -lg -lGLU

all: $(PROG)

$(PROG):	$(SRCS)
	$(CC) -o $(PROG) $(SRCS) $(LIBS)

clean:
	rm -f $(PROG)