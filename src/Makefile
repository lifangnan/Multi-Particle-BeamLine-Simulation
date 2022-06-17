include make_root.mk

DIR=beam/ beamline/ pywrapper/ server/ simulation/ util/
OBJS=$(foreach d, $(DIR), $(wildcard $(d)obj/*.o))
NONGRAPHICS_OBJS=$(foreach d, $(filter-out graphics/, $(DIR)), \
  $(wildcard $(d)obj/*.o))

../bin/HPSim.dll:
	nvcc -lcurand -L D:\ProgramFiles\anaconda2022\libs -I D:\ProgramFiles\anaconda2022\Lib\site-packages\numpy\core\include\numpy -I D:\ProgramFiles\anaconda2022\Library\include -shared \
		-o $@ $(NONGRAPHICS_OBJS)
#../bin/PyEPICS.so: server/obj/pyepics_module.o server/obj/ca_check.o
#	$(CC)  $(LDFLAGS) $(EPICS_LD_FLAGS) -shared -o $@ $^
# ../db/lib/libsqliteext.so: 
# 	gcc -shared -fPIC -o ../db/lib/libsqliteext.so ../db/lib/extension-functions.c -lm

# clean:
# 	@for d in $(DIR); \
# 	   do \
# 	     (cd $$d && (if [ -d obj ] ; then rm -rf obj ;fi)); \
# 	   done
# 	@if [ -d ../bin ]; then \
#            rm -rf ../bin; fi; 
# 	@(cd test && (if [ -f ./Makefile ] ; then $(MAKE) clean;fi));
# .PHONY: compile clean
