#include "beam.h"
#include "beam_cu.h"
#include "beamline.h"
#include "beamline_element.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include<cuda.h>
#include<string>
#include <cuda_runtime_api.h>
#include<algorithm>

#include "init.h"
#include "sql_utility.h"
#include "constant.h"

#include "space_charge.h"
#include "simulation_engine.h"
#include "simulation_engine_cu.h"
#include "plot_data.h"

#include <time.h>

using namespace std;


vector<double> getBeamMax(Beam* beam){
    vector<uint> loss = beam->GetLoss();
    vector<double> x = beam->GetX();
    vector<double> y = beam->GetY();
    double maxX = 0, maxY = 0;
    double tempX, tempY;

    for(int i=0; i<x.size(); i++){
        if(loss[i] == 0){
            tempX = abs(x[i]);
            tempY = abs(y[i]);
            if(tempX > maxX) maxX = tempX;
            if(tempY > maxY) maxY = tempY;
        }
    }
    vector<double> maxxy = { maxX,  maxY};
    return maxxy;
}

void simulate_and_save_file( Beam* beam, BeamLine beamline, int particle_num=1024, string filename = "sig_envelope"){
    vector<string> element_names = beamline.GetElementNames();
    Scheff* spacecharge = new Scheff(32, 128, 3);
    spacecharge->SetInterval(0.025);
    // spacecharge->SetAdjBunchCutoffW(0.8);
    spacecharge->SetRemeshThreshold(0.02);

    SimulationEngine* simulator = new SimulationEngine();
    PlotData* pltdata = new PlotData(particle_num);
    simulator->InitEngine(beam, &beamline, spacecharge, false, pltdata);
    string r_on_off = "off";
    simulator->SetSpaceCharge(r_on_off);

    ofstream out_file(filename.c_str());

    double position = 0.0;
    vector<double> maxXY;
    beam->UpdateSigX();
    beam->UpdateSigY();
    beam->UpdateLoss();
    // vector<double> vec_x, vec_y;
    // vec_x = beam->GetX();
    // vec_y = beam->GetY();
    out_file<< "position(m)" <<" "<<"Sig_x"<< " "<<"Sig_y"<<" "<<"Loss_num"<<" "<<"Type"<<endl;
    maxXY = getBeamMax(beam);
    out_file<< "0.0" <<" "<<maxXY[0]<< " "<<maxXY[1]<<" "<<beam->GetLossNum()<<" "<<"Begin"<<endl;
    // out_file<< "0.0" <<" "<<beam->GetSigX(true)<< " "<<beam->GetSigY(true)<<" "<<beam->GetLossNum()<<" "<<"Begin"<<endl;

    // 不分段的模拟
    clock_t start_time, end_time;
    start_time = clock();
    // simulator->Simulate(element_names.front(), element_names.back());
    // end_time = clock();
    // printf("Total time in function: %lf", double(end_time-start_time)/CLOCKS_PER_SEC );
    
    for(string element: element_names){
        BeamLineElement* temp_element = beamline[element];
        position += temp_element->GetLength();
        simulator->Simulate(element, element);
        beam->UpdateSigX();
        beam->UpdateSigY();
        beam->UpdateLoss();
        maxXY = getBeamMax(beam);
        // out_file<< position <<" "<<maxXY[0]<< " "<<maxXY[1]<<" "<<beam->GetLossNum()<<" "<<temp_element->GetType()<<endl;
        // out_file<< position <<" "<<beam->GetSigX(true)<< " "<<beam->GetSigY(true)<<" "<<beam->GetLossNum()<<" "<<temp_element->GetType()<<endl;

        // cout<<endl<<element<<" "<<temp_element->GetType()<<endl;
        // beam->UpdateAvgX();
        // cout<<"AvgX:" << beam->GetAvgX()<<endl;
        // beam->UpdateAvgY();
        // cout<<"AvgY:" << beam->GetAvgY()<<endl;
        // cout<<"SigX:" << beam->GetSigX(true)<<endl;
        // cout<<"SigY:" << beam->GetSigY(true)<<endl;
    }
    end_time = clock();
    printf("Total time in function: %lf", double(end_time-start_time)/CLOCKS_PER_SEC );
}

void loadBeamlineFromDatFile(string filename, BeamLine& beamline){
    ifstream infile;
	infile.open(filename, ios::in);
	if (!infile.is_open())
	{
		cout << "读取文件失败" << endl;
		return;
	}
	string buf;

    vector<string> type_vector = {"DRIFT", "EDGE", "BEND", "QUAD", "SOLENOID"};

    vector<vector<string> > input_dat;
    string current_type;

    // 记录每种类型元器件出现的次数
    map<string, int> device_count;
    device_count[string("DRIFT")] = 0;
    device_count[string("QUAD")] = 0;
    device_count[string("BEND")] = 0;
    device_count[string("SOLENOID")] = 0;

	while (getline(infile,buf))
	{
        if( buf.size() == 0 || buf.at(0) == ';' || buf.at(0) == ':' ){
            continue;
        }

        vector<string> res;
        stringstream in_str(buf);
        string temp_str;
        while(in_str >> temp_str){
            res.push_back(temp_str);
        }

        string first_str = res[0];
        if( first_str.find(":") != first_str.npos){
            int split_index = first_str.find(":");
            if(first_str.size() <= split_index+1){
                res[0] = first_str.substr(0, first_str.size()-1);
            }else{
                res[0] = first_str.substr(split_index+1, first_str.size()-split_index-1);
                res.insert(res.begin(), first_str.substr(0, split_index));
            }
        }else{
            transform(first_str.begin(), first_str.end(), first_str.begin(), ::toupper);
            if(device_count.count(first_str) == 1){
                device_count[first_str] += 1;
                res.insert(res.begin(), first_str + to_string(device_count[first_str]) );
            }else{
                res.insert(res.begin(), first_str );
            }
        }

        current_type = res[1];
        transform(current_type.begin(), current_type.end(), current_type.begin(), ::toupper);
        res[1] = current_type;

        input_dat.push_back(res);
    }

    // BeamLine beamline = BeamLine();
    for(int i = 0; i < input_dat.size(); ++i){
        vector<string> one_line = input_dat[i];
        // cout<< one_line[0]<<","<< one_line[1]<<","<< one_line[2]<<","<< one_line[3]<<endl;
        current_type = one_line[1];

        // tracewin的dat文件单位为mm、T、T/m
        // 而该多粒子模拟算法单位为m、T、T/m
        if(current_type == "DRIFT"){
            Drift* drift_p = new Drift(one_line[0]);
            drift_p->SetLength( stod(one_line[2])/1000.0 );
            drift_p->SetAperture( stod(one_line[3])/1000.0 );
            beamline.AddElement(drift_p);
        }
        else if (current_type == "QUAD")
        {
            Quad* quad_p = new Quad(one_line[0]);
            quad_p->SetLength( stod(one_line[2])/1000.0 );
            quad_p->SetGradient( stod(one_line[3]) );
            quad_p->SetAperture( stod(one_line[4])/1000.0 );
            beamline.AddElement(quad_p);
        }
        else if (current_type == "EDGE")
        {
            // 在TraceWin里，EDGE-BEND-EDGE这样的组合表示偏转铁，单位需要注意
            vector<string> Bend_line = input_dat[i+1];
            vector<string> outEdge_line = input_dat[i+2];
            if(Bend_line[1] != "BEND" || outEdge_line[1] != "EDGE"){
                cout<< "Something wrong when reading EDGE."<<endl;
                continue;
            }
            i += 2;

            Dipole* dipole_p = new Dipole(Bend_line[0]);
            dipole_p->SetRadius( stod(Bend_line[3])/1000.0 );
            dipole_p->SetAngle( -stod(Bend_line[2])*RADIAN );
            dipole_p->SetHalfGap( stod(one_line[4])/1000.0/2.0 );
            dipole_p->SetEdgeAngleIn( stod(one_line[2])*RADIAN );
            dipole_p->SetEdgeAngleOut( stod(outEdge_line[2])*RADIAN );
            // dipole_p->SetEdgeAngleIn( -22.5*RADIAN );
            // dipole_p->SetEdgeAngleOut( -22.5*RADIAN );
            dipole_p->SetK1( stod(one_line[5]) );
            dipole_p->SetK2( stod(one_line[6]) );
            // dipole_p->SetK1( 0.0 );
            // dipole_p->SetK2( 0.0 );
            dipole_p->SetFieldIndex( stod(Bend_line[4]) );
            dipole_p->SetAperture( stod(Bend_line[5])/1000.0 );
            dipole_p->SetKineticEnergy( 100.0 );
            dipole_p->SetLength( abs(dipole_p->GetRadius() * dipole_p->GetAngle()) );
            beamline.AddElement(dipole_p);
        }
        else if (current_type == "SOLENOID")
        {
            Solenoid* solenoid_p = new Solenoid( one_line[0] );
            solenoid_p->SetLength( stod(one_line[2])/1000.0 );
            solenoid_p->SetField( stod(one_line[3]) );
            solenoid_p->SetAperture( stod(one_line[4])/1000.0 );
            beamline.AddElement(solenoid_p);
        }
    }
}

int main(){
    // int particle_num = 1024 * 32;

    // Beam* beam = new Beam(particle_num, 939.294, 1.0, 0.015);
    // // beam->InitDCBeam(0.095, 47.0, 0.00327,  -0.102, 60.0, 0.002514, 180.0, 0.0, 0.7518);
    // beam->InitWaterbagBeam(0, 0.01, 0.000015,0, 0.01, 0.000015,0, 65.430429, 0.05633529, 0, 4.611, 500, 1);
    // beam->SaveInitialBeam();
    // beam->RestoreInitialBeam();

    // // beam->PrintToFile("initbeam","message");

    // SetGPU(0);

    // BeamLine beamline = BeamLine();
    // // cout<< beamline.GetSize()<<endl;

    // // const string url = "../db/clapa1.db";
    // DBConnection* db = new DBConnection(std::string("../db/clapa1.db"));
    // db->LoadLib("../db/lib/libsqliteext");
    // db->PrintDBs();
    // // db->PrintLibs();

    // // GenerateBeamLine(beamline, db);
    
    // Drift* drift1 = new Drift("drift1");
    // drift1->SetLength(1);

    // Solenoid* sole1 = new Solenoid("Solenoid1");
    // sole1->SetLength(0.45);
    // // sole1->SetAperture(10);
    // sole1->SetField(0.5);

    // Drift* drift2 = new Drift("drift2");
    // drift2->SetLength(1);

    // Solenoid* sole2 = new Solenoid("Solenoid2");
    // sole2->SetLength(0.45);
    // // sole1->SetAperture(10);
    // sole2->SetField(-0.5);

    // Drift* drift3 = new Drift("drift3");
    // drift3->SetLength(1);

    // beamline.AddElement(drift1);
    // beamline.AddElement(sole1);
    // beamline.AddElement(drift2);
    // beamline.AddElement(sole2);
    // beamline.AddElement(drift3);

    // vector<string> beamlinenames = beamline.GetElementNames();
    // for(auto name: beamlinenames){
    //     cout<< name<<endl;
    // }

    // string filename = "beam_envelope";
    // simulate_and_save_file(beam, beamline, particle_num, filename);

    // vector<uint> loss_ind = beam->GetLoss();
    // cout<<loss_ind.size()<<endl;
    // for(uint ind: loss_ind){
    //     cout<<ind<<",";
    // }



    // Part2: Area for test
    // BeamLine bl = BeamLine();
    // loadBeamlineFromDatFile("clapa2.dat", bl);
    // // loadBeamlineFromDatFile("test_dipole.dat", bl);
    // bl.Print();
    // cout<< bl.GetSize()<<endl;
    // vector<string> beamlinenames = bl.GetElementNames();
    // for(auto name: beamlinenames){
    //     cout<< name<<endl;
    // }

    // set beam
    int particle_num = 10240;
    // Beam* beam = new Beam(particle_num, 939.294, 1.0, 0.015);
    // // beam->InitDCBeam(0.095, 47.0, 0.00327,  -0.102, 60.0, 0.002514, 180.0, 0.0, 0.7518);
    // beam->InitWaterbagBeam(0, 0.0001, 0.00000002484565,0, 0.0001, 0.00000002484565, 0, 8, 3.1415926*pow(10.0,-11), 0, 100, 500, 1);
    // beam->SaveInitialBeam();
    // beam->RestoreInitialBeam();
    BeamLine bl = BeamLine();
    loadBeamlineFromDatFile("yc_model.dat", bl);
    bl.Print();
    Beam* beam1 = new Beam(particle_num, 939.294, 1.0, 0.015);
    beam1->InitWaterbagBeam(0, 0.005, 0.00015,0, 0.005, 0.00015,0, 6.5430429, 0.00005633529, 0, 3, 500, 1);
    beam1->SaveInitialBeam();

    // beam->PrintToFile("initbeam_spacecharge","message");
    SetGPU(0);
    Scheff* spacecharge = new Scheff(32, 128, 3);
    spacecharge->SetInterval(0.025);
    // spacecharge->SetAdjBunchCutoffW(0.8);
    spacecharge->SetRemeshThreshold(0.02);

    cout<<"posi1"<<endl;
    SimulationEngine* simulator = new SimulationEngine();
    PlotData* pltdata = new PlotData(particle_num);
    simulator->InitEngine(beam1, &bl, spacecharge, false, pltdata);
    string r_on_off = "off";
    simulator->SetSpaceCharge(r_on_off);
    cout<<"posi2"<<endl;
    simulator->Simulate("DRIFT1", "DRIFT3");
    cout<<"posi3"<<endl;
    beam1->PrintToFile("./beam_file/beam_1");
    simulator->Simulate("Q1", "Q1");
    beam1->PrintToFile("./beam_file/beam_Q1");
    simulator->Simulate("DRIFT4", "Q2");
    beam1->PrintToFile("./beam_file/beam_Q2");
    simulator->Simulate("DRIFT5", "Q3");
    beam1->PrintToFile("./beam_file/beam_Q3");
    simulator->Simulate("DRIFT6", "DRIFT63");
    beam1->PrintToFile("./beam_file/beam_end");

    // string filename = "beam_envelope";
    // clock_t start_time, end_time;
    // start_time = clock();
    // simulate_and_save_file(beam, bl, particle_num, filename);
    // end_time = clock();
    // printf("Total time: %lf", double(end_time-start_time)/CLOCKS_PER_SEC );
    // beam->PrintToFile("endbeam_spacecharge","message");

    return 0;
}