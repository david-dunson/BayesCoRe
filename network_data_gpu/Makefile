CC=gcc

all: opencl

cuda: main.cpp
	g++ -std=c++11 -g main.cpp -o runner -I. -I/usr/local/include/ -L/usr/local/lib -lafcuda
	
af_to_csv: main.cpp
	g++ -std=c++11 -g af_to_csv.cpp -o runner_af_to_csv -I. -I/usr/local/include/ -L/usr/local/lib -lafcuda

opencl: main.cpp
	g++ -std=c++11 -g main.cpp -o runner -I. -I/usr/local/include/ -L/usr/local/lib -lafopencl -lOpenCL
	
cpu: main.cpp
	g++ -std=c++11 -g main.cpp -o runner -I. -I/usr/local/include/ -L/usr/local/lib -lafcpu
	
