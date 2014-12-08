GRAPH = gnuplot

CC = clang
CFLAGS = -Wall -O0 -g 
LFLAGS = -O0 -g
LIBS = -lgsl  -lgslcblas -lm

star: star-controller.o star-ode.o 
	${CC} $(LFLAGS) -o $@ $^ $(LIBS)

star_data: star
	./star > star_data

star.png: star.gpl star_data
	$(GRAPH) star.gpl

clean : 
	rm -f *~
	rm -f *.o
	rm -f star

veryclean : clean
	rm -f res1 star.png
