CC=nvcc
CXX=nvcc
PROJECT_ROOT:=$(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

#EPICS_INC_FLAGS=-I/ade/epics/supTop/base/R3.14.11/include/os/Linux \
  -I/ade/epics/supTop/base/R3.14.11/include -D_POSIX_C_SOURCE=199506L \
  -D_POSIX_THREADS -D_XOPEN_SOURCE=500 -D_X86_64_ -DUNIX -D_BSD_SOURCE \
  -Dlinux -D_REENTRANT -m64 

#EPICS_LD_FLAGS=-L/ade/epics/supTop/base/R3.14.11/lib/linux-x86_64 \
  -lcas -lgdd -lasHost -ldbStaticHost -lregistryIoc -lca -lCom \
  -lpthread -lreadline -lncurses -lm -lrt -ldl -lgcc

#GL_LD_FLAGS+=-lglut -lGLU -lGLEW 
PYTHON_INC_FLAGS=-I/home/laser/anaconda2/include/python2.7 \
  -I/home/laser/anaconda2/lib/python2.7/site-packages/numpy/core/include/numpy

PYTHON_LD_FLAGS=-L/home/laser/anaconda2/lib \
  -L/home/laser/anaconda2/lib/python2.7/site-packages/numpy/core/lib -lpython2.7

CPPFLAGS+=-arch=sm_50 -Xcompiler '-fPIC' -Xcompiler '-fopenmp' -O3 \
  -g -w #-DDOUBLE_PRECISION #-D_DEBUG

#CPPFLAGS+=$(EPICS_INC_FLAGS)
CPPFLAGS+=$(PYTHON_INC_FLAGS)
#CPPFLAGS+=-Xptxas -v -Xptxas -dlcm=ca

LDFLAGS=-lcurand -lsqlite3
