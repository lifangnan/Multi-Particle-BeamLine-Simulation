import sys
import os
# define directory to packages and append to $PATH
par_dir = "/home/laser/hpsim"
print(par_dir)
lib_dir = os.path.join(par_dir,"bin")
print(lib_dir)
sys.path.append(lib_dir)
pkg_dir = os.path.join(par_dir,"pylib")
print(pkg_dir)
sys.path.append(pkg_dir)

import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
import math
# import additional simulation packages
import hpsim as hps
import HPSim as HPSim

GPU = 0
hps.set_gpu(GPU)

db_dir = par_dir + '/db'
lib_dir = par_dir + '/db/lib/'
dbs = ['tbtd.db','dtl.db','trst.db','ccl.db']
dbconn1 = hps.DBConnection(db_dir, dbs, lib_dir, 'libsqliteext')
dbconn1.print_libs()
dbconn1.print_dbs()
dbconn1.clear_model_index()
print("*** dB connection established ***")

bl = hps.BeamLine()
beamline = hps.get_element_list()
print(beamline)
print("*** Beamline created ***\n")

################################################################################
# create table of beamline elements at lengths
pybl = pydb.Db_bl(db_dir, dbs)
py_beamline = pybl.get_bl_elem_len()
print("*** PySQLite Beamline created ***\n")