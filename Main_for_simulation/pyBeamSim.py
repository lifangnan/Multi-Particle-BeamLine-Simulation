from ctypes import * 
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns

dll = CDLL("F:\git_workspace\Multi-Particle-BeamLine-Simulation\Main_for_simulation\simulation_for_python.dll")

# 定义参数和返回值的类型
dll.init_beam.argtypes = [c_int, c_double, c_double, c_double]

dll.beam_init_from_file.argtypes = [c_char_p]

dll.beam_print_to_file.argtypes =[c_char_p, c_char_p]

dll.set_beamTwiss.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_uint]

dll.getBeamAvgx.restype = c_double
dll.getBeamAvgy.restype = c_double
dll.getBeamSigx.restype = c_double
dll.getBeamSigy.restype = c_double
dll.getBeamMaxx.restype = c_double
dll.getBeamMaxy.restype = c_double

dll.init_database.argtypes = [c_char_p]

dll.load_Beamline_From_DatFile.argtypes = [c_char_p]

dll.get_Beamline_ElementNames.restype = c_char_p
dll.get_Beamline_ElementTypes.restype = c_char_p
dll.get_Beamline_ElementLengths.restype = c_char_p

dll.init_spacecharge.argtypes = [c_uint, c_uint, c_int]

dll.simulate_and_getEnvelope.argtypes = [c_bool, c_bool]
dll.simulate_and_getEnvelope.restype = POINTER(POINTER(c_double))
dll.get_envelope_size.restype = c_int

dll.set_magnet_with_index.argtypes = [c_int, c_double]

dll.set_magnet_with_name.argtypes = [c_char_p, c_double]

dll.get_good_number.restype = c_int

class BeamSimulator():
    def __init__(self):
        dll.new_my_simulator()

    def init_Default_Simulator(self):
        dll.default_init()
    
    def init_beam(self, particle_num, rest_energy, charge, current):
        dll.init_beam(particle_num, rest_energy, charge, current)
    
    def beam_init_from_file(self, filename):
        dll.beam_init_from_file(filename.encode())
    
    def beam_print_to_file(self, filePath, comment = "CLAPA"):
        dll.beam_print_to_file(filePath.encode(), comment.encode())
    
    def set_beamTwiss(self, r_ax, r_bx, r_ex, r_ay, r_by, r_ey, r_az, r_bz, r_ez, r_sync_phi, r_sync_w, r_freq, r_seed = 1):
        dll.set_beamTwiss(r_ax, r_bx, r_ex, r_ay, r_by, r_ey, r_az, r_bz, r_ez, r_sync_phi, r_sync_w, r_freq, r_seed)
    
    def save_initial_beam(self):
        dll.save_initial_beam()
    
    def restore_initial_beam(self):
        dll.restore_initial_beam()

    def getBeamAvgx(self):
        return dll.getBeamAvgx()
    
    def getBeamAvgy(self):
        return dll.getBeamAvgy()

    def getBeamSigx(self):
        return dll.getBeamSigx()

    def getBeamSigy(self):
        return dll.getBeamSigy()

    def getBeamMaxx(self):
        return dll.getBeamMaxx()

    def getBeamMaxy(self):
        return dll.getBeamMaxy()

    def free_beam(self):
        dll.free_beam()

    def load_Beamline_From_Sqlite(self, DB_Path):
        dll.init_database(DB_Path.encode())
        dll.init_beamline_from_DB()

    def load_Beamline_From_DatFile(self, filename):
        dll.load_Beamline_From_DatFile(filename.encode())

    def get_Beamline_ElementNames(self):
        names_str = (dll.get_Beamline_ElementNames()).decode()
        names_list = names_str.split(",")
        return names_list
    
    def get_Beamline_ElementTypes(self):
        types_str = (dll.get_Beamline_ElementTypes()).decode()
        types_list = types_str.split(",")
        return types_list
    
    def get_Beamline_ElementLengths(self):
        lengths_str = (dll.get_Beamline_ElementLengths()).decode()
        lengths_list = lengths_str.split(",")
        lengths_list = [float(item) for item in lengths_list]
        return np.array(lengths_list)


    def init_spacecharge(self, r_nr = 32, r_nz = 128, r_adj_bunch = 3):
        dll.init_spacecharge(r_nr, r_nz, r_adj_bunch)
    
    def simulate_and_getEnvelope(self, use_spacecharge = False, is_sig = False):
        envelope = dll.simulate_and_getEnvelope(use_spacecharge, is_sig)
        envelope_size = dll.get_envelope_size()

        np_envelope = []
        for i in range(envelope_size):
            temp = []
            for j in range(4):
                temp.append(envelope[i][j])
            np_envelope.append(temp)
        np_envelope = np.array(np_envelope)
        return np_envelope
    
    def simulate_all(self, use_spacecharge = False, is_sig = False):
        dll.simulate(use_spacecharge, is_sig)

    def set_magnet_with_index(self, magnet_index, field_or_angle):
        dll.set_magnet_with_index(magnet_index, field_or_angle)
    
    def set_magnet_with_name(self, element_name, field_or_angle):
        dll.set_magnet_with_name(element_name.encode(), field_or_angle)
        

    def plot_envelope(self, envelope):
        element_types = self.get_Beamline_ElementTypes()
        element_lengths = self.get_Beamline_ElementLengths()
        position_start = np.array([])
        position_end = np.array([])
        for i in range(element_lengths.shape[0]):
            position_start = np.append(position_start, element_lengths[:i].sum())
            position_end = np.append(position_end, position_start[i] + element_lengths[i])

        plt.figure(figsize=(13,3))
        plt.plot(envelope[:,0], envelope[:,1])
        plt.plot(envelope[:,0], envelope[:,2])

        max_enve = envelope[:, [1,2]].max().max()

        deviceType_list = ["Dipole", "Solenoid", "Quad"]
        for i in range(len(element_types)):
            if element_types[i] in deviceType_list:
                element_start = position_start[i]
                element_end = position_end[i]
                if element_types[i] == "Dipole":
                    plt.fill_between([element_start, element_end],0, max_enve,facecolor = 'pink', alpha = 0.9)
                elif element_types[i] == "Solenoid":
                    plt.fill_between([element_start, element_end],0, max_enve,facecolor = 'green', alpha = 0.3)
                elif element_types[i] == "Quad":
                    plt.fill_between([element_start, element_end],0, max_enve,facecolor = 'blue', alpha = 0.3)   

        plt.show()
    
    def get_good_num(self):
        return dll.get_good_number()

    def plot_beam(self):
        cwd = os.getcwd()
        self.beam_print_to_file(cwd + "\\temp_Beam")
        beam_data = np.loadtxt(cwd + "\\temp_Beam")

        beam_data = beam_data[beam_data[:,-2]==0, :]
        intibeam_df = pd.DataFrame(columns=['x','y'])
        intibeam_df['x'] = beam_data[:, 0]
        intibeam_df['y'] = beam_data[:, 2]
        sns.jointplot(x="x", y="y", data=intibeam_df, kind="kde", levels=50, fill=True, cmap='binary')

        plt.figure(figsize=[15,5])
        plt.subplot(1,3,1)
        plt.scatter(beam_data[:, 0], beam_data[:, 2], s=1)
        plt.axis("equal")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("x-y")

        plt.subplot(1,3,2)
        plt.scatter(beam_data[:, 0], beam_data[:, 1], s=1)
        plt.xlabel("x")
        plt.ylabel("px")
        plt.title("x-px")

        plt.subplot(1,3,3)
        plt.scatter(beam_data[:, 2], beam_data[:, 3], s=1)
        plt.xlabel("y")
        plt.ylabel("py")
        plt.title("y-py")

        plt.show()

        os.remove(cwd + "\\temp_Beam")