all: src/main.cpp
	c++ -std=c++11 -I./include -O3 src/main.cpp -o raess
clean:
	rm -f raess
