# SOURCE = $(addprefix obj/, $(wildcard *.o))
SOURCE = $(wildcard obj/*.obj)
# cpp_source = $(wildcard *.cpp)
# target = $(basename *.cpp)
simulation_for_python.dll:simulation_for_python.cpp
	nvcc -c $< -o simulation_for_python.obj -I D:\ProgramFiles\anaconda2022\include -I D:\ProgramFiles\anaconda2022\Library\include -I./inc 
	nvcc -shared -I D:\ProgramFiles\anaconda2022\Library\include -L D:\ProgramFiles\anaconda2022\libs -I./inc $(SOURCE) simulation_for_python.obj -o $@
