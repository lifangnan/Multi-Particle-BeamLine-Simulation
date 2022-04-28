# SOURCE = $(addprefix obj/, $(wildcard *.o))
SOURCE = $(wildcard obj/*.o)
# cpp_source = $(wildcard *.cpp)
# target = $(basename *.cpp)
ctest:ctest.cpp
	nvcc -c $< -o ctest.o -I D:\ProgramFiles\anaconda2022\include -I D:\ProgramFiles\anaconda2022\Library\include -I./inc 
	nvcc -I D:\ProgramFiles\anaconda2022\Library\include -L D:\ProgramFiles\anaconda2022\libs -I./inc $(SOURCE) ctest.o -o $@
