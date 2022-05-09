from ctypes import * 

mysim = CDLL("Main_for_simulation\simulation_for_python.dll")
# lib.run_simulation()
mysim.new_my_simulator()
# lib.default_init()
mysim.init_beam(c_int(10240), c_double(939.294), c_double(1.0), c_double(0.015))
mysim.set_beamTwiss(c_double(0), c_double(0.01), c_double(0.000015), c_double(0), c_double(0.01), c_double(0.000015), c_double(0), c_double(65.430429), c_double(0.05633529), c_double(0), c_double(4.611), c_double(500), c_uint(1))
file = c_char_p(b"F:/git_workspace/Multi-Particle-BeamLine-Simulation/db/clapa1.db")
mysim.init_database(file)

mysim.init_beamline_from_DB()

mysim.init_spacecharge.argtypes = [c_uint, c_uint, c_int]
mysim.init_spacecharge(32, 128, 3)

mysim.simulate_and_getEnvelope.restype = POINTER(POINTER(c_double))
res = mysim.simulate_and_getEnvelope()