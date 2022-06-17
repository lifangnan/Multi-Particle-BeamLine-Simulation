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

void simulate_and_save_file( Beam* beam, BeamLine beamline, int particle_num=1024, string filename = "sig_envelope"){
    vector<string> element_names = beamline.GetElementNames();
    Scheff* spacecharge = new Scheff(32, 128, 3);
    spacecharge->SetInterval(0.025);
    spacecharge->SetAdjBunchCutoffW(0.8);
    spacecharge->SetRemeshThreshold(0.02);

    SimulationEngine* simulator = new SimulationEngine();
    PlotData* pltdata = new PlotData(particle_num);
    simulator->InitEngine(beam, &beamline, spacecharge, false, pltdata);

    ofstream out_file(filename.c_str());

    double position = 0.0;
    beam->UpdateSigX();
    beam->UpdateSigY();
    beam->UpdateLoss();
    out_file<< "0.0" <<" "<<beam->GetSigX(true)<< " "<<beam->GetSigY(true)<<" "<<beam->GetLossNum()<<endl;

    for(string element: element_names){
        BeamLineElement* temp_element = beamline[element];
        position += temp_element->GetLength();
        simulator->Simulate(element, element);
        beam->UpdateSigX();
        beam->UpdateSigY();
        beam->UpdateLoss();
        out_file<< position <<" "<<beam->GetSigX(true)<< " "<<beam->GetSigY(true)<<" "<<beam->GetLossNum()<<endl;

        cout<<endl<<element<<" "<<temp_element->GetType()<<endl;
        beam->UpdateAvgX();
        cout<<"AvgX:" << beam->GetAvgX()<<endl;
        beam->UpdateAvgY();
        cout<<"AvgY:" << beam->GetAvgY()<<endl;
        cout<<"SigX:" << beam->GetSigX(true)<<endl;
        cout<<"SigY:" << beam->GetSigY(true)<<endl;
    }
}

int main(){
    int particle_num = 1024 * 32;

    Beam* beam = new Beam(particle_num, 939.294, 1.0, 0.015);
    // beam->InitDCBeam(0.095, 47.0, 0.00327,  -0.102, 60.0, 0.002514, 180.0, 0.0, 0.7518);
    beam->InitWaterbagBeam(0, 0.01, 0.000015,0, 0.01, 0.000015,0, 65.430429, 0.05633529, 0, 4.611, 500, 1);
    beam->SaveInitialBeam();
    beam->RestoreInitialBeam();

    beam->PrintToFile("initbeam","message");

    SetGPU(0);

    BeamLine beamline = BeamLine();
    // cout<< beamline.GetSize()<<endl;

    // const string url = "../db/clapa1.db";
    DBConnection* db = new DBConnection(std::string("../db/clapa1.db"));
    db->LoadLib("../db/lib/libsqliteext");
    db->PrintDBs();
    // db->PrintLibs();

    GenerateBeamLine(beamline, db);

    string filename = "beam_envelope";
    simulate_and_save_file(beam, beamline, particle_num, filename);

    // vector<uint> loss_ind = beam->GetLoss();
    // cout<<loss_ind.size()<<endl;
    // for(uint ind: loss_ind){
    //     cout<<ind<<",";
    // }


    return 0;
}