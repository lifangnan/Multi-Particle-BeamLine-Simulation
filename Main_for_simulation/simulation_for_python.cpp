#define EXPORT __declspec(dllexport)

#include "beam.h"
#include "beam_cu.h"
#include "beamline.h"
#include "beamline_element.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<cuda.h>
#include<string>
#include <cuda_runtime_api.h>
#include<algorithm>

#include "init.h"
#include "sql_utility.h"

#include "space_charge.h"
#include "simulation_engine.h"
#include "simulation_engine_cu.h"
#include "plot_data.h"

using namespace std;

// void modify_beamline(beamline& bl, )

Beam* new_beam(int particle_num, double rest_energy, double charge, double current){
    Beam* beam = new Beam(particle_num, rest_energy, charge, current);
    return beam;
}

class my_simulator{
    public:
        Beam* beam;
        Beam* origin_beam;
        BeamLine beamline;
        DBConnection* db;
        Scheff* spacecharge;
        SimulationEngine* simulator;
        int beamline_size;
        int envelope_size;
        int particle_size;

        my_simulator();
        void default_init_c();

        void init_beam_c(int particle_num, double rest_energy, double charge, double current);
        void beam_init_from_file_c(string file_path);
        void beam_print_to_file_c(string file_path, string comment = "CLAPA");
        void set_beamTwiss_c(double r_ax, double r_bx, double r_ex, 
                double r_ay, double r_by, double r_ey, double r_az, double r_bz, double r_ez,
                double r_sync_phi, double r_sync_w, double r_freq, unsigned int r_seed = 1);

        void init_database_c(string DB_Url, string lib_Url = "F:/git_workspace/Multi-Particle-BeamLine-Simulation/db/lib/libsqliteext");
        void init_beamline_from_DB_c();
        void init_spacecharge_c(uint r_nr = 32, uint r_nz = 128, int r_adj_bunch = 3);

        double **simulate_and_getEnvelope_c();
        // double[][] get_origin_beam();
        // double[][] get_beam();

        int get_envelope_size();
        
        void free_all_c();
};

my_simulator::my_simulator(){
    spacecharge = new Scheff(32, 128, 3);
}

void my_simulator::default_init_c(){
    init_beam_c(1024, 939.294, 1.0, 0.015);
    set_beamTwiss_c(0, 0.01, 0.000015,0, 0.01, 0.000015,0, 65.430429, 0.05633529, 0, 4.611, 500, 1);
    init_database_c( string("F:/git_workspace/Multi-Particle-BeamLine-Simulation/db/clapa1.db") );
    init_beamline_from_DB_c();
    init_spacecharge_c(32, 128, 3);
}

void my_simulator::init_beam_c(int particle_num, double rest_energy, double charge, double current){
    beam = new Beam(particle_num, rest_energy, charge, current);
    particle_size = particle_num;
}
void my_simulator::beam_init_from_file_c(string file_path){
    if(beam){
        beam->InitBeamFromFile(file_path);
    }else{
        cout<<"Please use init_beam() first."<<endl;
    }
}
void my_simulator::beam_print_to_file_c(string file_path, string comment){
    if(beam){
        beam->PrintToFile(file_path, comment);
    }else{
        cout<<"No beam exists."<<endl;
    }
}
void my_simulator::set_beamTwiss_c(double r_ax, double r_bx, double r_ex, 
    double r_ay, double r_by, double r_ey, double r_az, double r_bz, double r_ez,
    double r_sync_phi, double r_sync_w, double r_freq, unsigned int r_seed){
        if(beam) beam->InitWaterbagBeam(r_ax, r_bx, r_ex, r_ay, r_by, r_ey, r_az, r_bz, r_ez, r_sync_phi, r_sync_w, r_freq, r_seed);
}

void my_simulator::init_database_c(string DB_Url, string lib_Url){
    db = new DBConnection(DB_Url);
    db->LoadLib(lib_Url);
    db->PrintDBs();
}

void my_simulator::init_beamline_from_DB_c(){
    beamline = BeamLine();
    GenerateBeamLine(beamline, db);
    beamline_size = beamline.GetElementNames().size();
}

void my_simulator::init_spacecharge_c(uint r_nr, uint r_nz, int r_adj_bunch){
    spacecharge = new Scheff(r_nr, r_nz, r_adj_bunch);
    spacecharge->SetInterval(0.025);
    spacecharge->SetAdjBunchCutoffW(0.8);
    spacecharge->SetRemeshThreshold(0.02);
}

double **my_simulator::simulate_and_getEnvelope_c(){
    double **rt_envelope;
    int n = 10240, m = 4;
    rt_envelope = (double **)malloc(10240 * sizeof(double *));
    for(int i=0; i<n; i++){
        rt_envelope[i] = (double *)malloc(m* sizeof(double));
    }
    // double (*rt_envelope)[4] = new double[10240][4];

    SetGPU(0);

    simulator = new SimulationEngine();
    PlotData* pltdata = new PlotData(particle_size);
    simulator->InitEngine(beam, &beamline, spacecharge, false, pltdata);

    // ofstream out_file(filename.c_str());
    double position = 0.0;
    beam->UpdateSigX();
    beam->UpdateSigY();
    beam->UpdateLoss();
    // out_file<< "0.0" <<" "<<beam->GetSigX(true)<< " "<<beam->GetSigY(true)<<" "<<beam->GetLossNum()<<endl;
    rt_envelope[0][0] = 0.0; rt_envelope[0][1] = beam->GetSigX(true);
    rt_envelope[0][2] = beam->GetSigY(true); rt_envelope[0][3] = beam->GetLossNum();
    envelope_size = 1;
    
    vector<string> element_names = beamline.GetElementNames();
    for(string element: element_names){
        BeamLineElement* temp_element = beamline[element];
        position += temp_element->GetLength();
        simulator->Simulate(element, element);
        beam->UpdateSigX();
        beam->UpdateSigY();
        beam->UpdateLoss();
        // out_file<< position <<" "<<beam->GetSigX(true)<< " "<<beam->GetSigY(true)<<" "<<beam->GetLossNum()<<endl;
        rt_envelope[envelope_size][0] = position; rt_envelope[envelope_size][1] = beam->GetSigX(true);
        rt_envelope[envelope_size][2] = beam->GetSigY(true); rt_envelope[envelope_size][3] = beam->GetLossNum();
        envelope_size++;

        cout<<element<<" "<<temp_element->GetType()<<endl;
        // beam->UpdateAvgX();
        // cout<<"AvgX:" << beam->GetAvgX()<<endl;
        // beam->UpdateAvgY();
        // cout<<"AvgY:" << beam->GetAvgY()<<endl;
        // cout<<"SigX:" << beam->GetSigX(true)<<endl;
        // cout<<"SigY:" << beam->GetSigY(true)<<endl;
    }

    return rt_envelope;
}

int my_simulator::get_envelope_size(){
    return envelope_size;
}

void my_simulator::free_all_c(){
    
}


extern "C" {
    my_simulator* mysim;
    EXPORT void new_my_simulator(){
        mysim = new my_simulator();
    }

    EXPORT void default_init(){
        mysim->default_init_c();
    }

    EXPORT void init_beam(int particle_num, double rest_energy, double charge, double current){
        mysim->init_beam_c(particle_num, rest_energy, charge, current);
    }

    EXPORT void beam_init_from_file(string file_path){
        mysim->beam_init_from_file_c(file_path);
    }  

    EXPORT void beam_print_to_file(string file_path, string comment = "CLAPA"){
        mysim->beam_print_to_file_c(file_path, comment);
    }  

    EXPORT void set_beamTwiss(double r_ax, double r_bx, double r_ex, 
                double r_ay, double r_by, double r_ey, double r_az, double r_bz, double r_ez,
                double r_sync_phi, double r_sync_w, double r_freq, unsigned int r_seed = 1){
        mysim->set_beamTwiss_c(r_ax, r_bx, r_ex, r_ay, r_by, r_ey, r_az, r_bz, r_ez, r_sync_phi, r_sync_w, r_freq, r_seed);
    }  

    EXPORT void init_database(char* DB_Url){
        // cout<<"start"<<endl;
        // cout<<string(DB_Url)<<endl;
        mysim->init_database_c(string(DB_Url));
    }  

    EXPORT void init_beamline_from_DB(){
        mysim->init_beamline_from_DB_c();
    }  

    EXPORT void init_spacecharge(uint r_nr = 32, uint r_nz = 128, int r_adj_bunch = 3){
        mysim->init_spacecharge_c(r_nr, r_nz, r_adj_bunch);
    }  

    EXPORT double **simulate_and_getEnvelope(){
        return mysim->simulate_and_getEnvelope_c();
    }   

    EXPORT int get_envelope_size(){
        return mysim->get_envelope_size();
    }
}