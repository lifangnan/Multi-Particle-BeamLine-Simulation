#define EXPORT __declspec(dllexport)

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

class my_simulator{
    public:
        Beam* beam;
        Beam* beam_init;
        BeamLine beamline;
        DBConnection* db;
        Scheff* spacecharge;
        SimulationEngine* simulator;
        int beamline_size;
        int envelope_size;
        int particle_size;

        char *rt_names, *rt_types, *rt_lengths;

        my_simulator();
        void default_init_c();

        void init_beam_c(int particle_num, double rest_energy, double charge, double current);
        void beam_init_from_file_c(string file_path);
        void beam_print_to_file_c(string file_path, string comment = "CLAPA");
        void set_beamTwiss_c(double r_ax, double r_bx, double r_ex, 
                double r_ay, double r_by, double r_ey, double r_az, double r_bz, double r_ez,
                double r_sync_phi, double r_sync_w, double r_freq, unsigned int r_seed = 1);
        void save_initial_beam_c();
        void restore_initial_beam_c();
        double getBeamAvgx_c();
        double getBeamAvgy_c();
        double getBeamSigx_c();
        double getBeamSigy_c();
        double getBeamMaxx_c();
        double getBeamMaxy_c();


        void free_beam_c();

        void init_database_c(string DB_Url, string lib_Url = "F:/git_workspace/Multi-Particle-BeamLine-Simulation/db/lib/libsqliteext");
        void init_beamline_from_DB_c();
        void load_Beamline_From_DatFile_c(string filename);
        char* get_Beamline_ElementNames_c();
        char* get_Beamline_ElementTypes_c();
        char* get_Beamline_ElementLengths_c();

        void init_spacecharge_c(uint r_nr = 32, uint r_nz = 128, int r_adj_bunch = 3);

        double **rt_envelope;
        double **simulate_and_getEnvelope_c(bool use_spacecharge = false, bool is_sig = false);
        void simulate_c(bool use_spacecharge = false, bool is_sig = false);
        // double[][] get_beam();

        int get_envelope_size();

        void set_magnet_with_index_c(int magnet_index, double field_or_angle);
        void set_magnet_with_name_c(string elementname, double field_or_angle);

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
    if(beam) delete beam;
    beam = new Beam(particle_num, rest_energy, charge, current);
    // beam_init = new Beam(particle_num, rest_energy, charge, current);
    particle_size = particle_num;
}
void my_simulator::beam_init_from_file_c(string file_path){
    if(beam){
        beam->InitBeamFromFile(file_path);
        // beam_init->InitBeamFromFile(file_path);
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
        if(beam){
            beam->InitWaterbagBeam(r_ax, r_bx, r_ex, r_ay, r_by, r_ey, r_az, r_bz, r_ez, r_sync_phi, r_sync_w, r_freq, r_seed);
        }
        // if(beam_init) beam_init->InitWaterbagBeam(r_ax, r_bx, r_ex, r_ay, r_by, r_ey, r_az, r_bz, r_ez, r_sync_phi, r_sync_w, r_freq, r_seed);
}

void my_simulator::save_initial_beam_c(){
    if(beam){
        beam->SaveInitialBeam();
    }
}
void my_simulator::restore_initial_beam_c(){
    if(beam){
        beam->RestoreInitialBeam();
    }
}

double my_simulator::getBeamAvgx_c(){
    vector<uint> loss = beam->GetLoss();
    vector<double> x = beam->GetX();
    double sumX = 0;
    int count_goodnum = 0;

    for(int i=0; i<x.size(); i++){
        if(loss[i] == 0){
            count_goodnum++;
            sumX += x[i];
        }
    }
    return sumX / count_goodnum;
}

double my_simulator::getBeamAvgy_c(){
    vector<uint> loss = beam->GetLoss();
    vector<double> y = beam->GetY();
    double sumY = 0;
    int count_goodnum = 0;

    for(int i=0; i<y.size(); i++){
        if(loss[i] == 0){
            count_goodnum++;
            sumY += y[i];
        }
    }
    return sumY / count_goodnum;
}

double my_simulator::getBeamSigx_c(){
    vector<uint> loss = beam->GetLoss();
    vector<double> x = beam->GetX();
    double sumX = 0;
    int count_goodnum = 0;

    for(int i=0; i<x.size(); i++){
        if(loss[i] == 0){
            count_goodnum++;
            sumX += x[i];
        }
    }
    double avgX = sumX / count_goodnum;
    double temp_sum = 0;
    for(int i=0; i<x.size(); i++){
        if(loss[i] == 0){
            temp_sum += (x[i] - avgX)*(x[i] - avgX);
        }
    }
    if(count_goodnum > 1){
        double stdX = sqrt(temp_sum/(count_goodnum-1));
        return stdX;
    }else{
        return 0.0;
    } 
}

double my_simulator::getBeamSigy_c(){
    vector<uint> loss = beam->GetLoss();
    vector<double> Y = beam->GetY();
    double sumY = 0;
    int count_goodnum = 0;

    for(int i=0; i<Y.size(); i++){
        if(loss[i] == 0){
            count_goodnum++;
            sumY += Y[i];
        }
    }
    double avgY = sumY / count_goodnum;
    double temp_sum = 0;
    for(int i=0; i<Y.size(); i++){
        if(loss[i] == 0){
            temp_sum += (Y[i] - avgY)*(Y[i] - avgY);
        }
    }
    if(count_goodnum > 1){
        double stdY = sqrt(temp_sum/(count_goodnum-1));
        return stdY;
    }else{
        return 0.0;
    } 
}

double my_simulator::getBeamMaxx_c(){
    vector<uint> loss = beam->GetLoss();
    vector<double> x = beam->GetX();
    double maxX = 0, tempX = 0;

    for(int i=0; i<x.size(); i++){
        if(loss[i] == 0){
            tempX = abs(x[i]);
            if(tempX > maxX) maxX = tempX;
        }
    }
    return maxX;
}

double my_simulator::getBeamMaxy_c(){
    vector<uint> loss = beam->GetLoss();
    vector<double> y = beam->GetY();
    double maxY = 0, tempY = 0;

    for(int i=0; i<y.size(); i++){
        if(loss[i] == 0){
            tempY = abs(y[i]);
            if(tempY > maxY) maxY = tempY;
        }
    }
    return maxY;
}

void my_simulator::free_beam_c(){
    if(!beam) delete beam;
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

void my_simulator::load_Beamline_From_DatFile_c(string filename){
    beamline = BeamLine();
    ifstream infile;
	infile.open(filename, ios::in);
	if (!infile.is_open())
	{
		cout << "读取文件失败" << endl;
		return;
	}
	string buf;

    vector<string> type_vector = {"DRIFT", "EDGE", "BEND", "QUAD", "SOLENOID", "ApertureCircular", "ApertureRectangular"};

    vector<vector<string> > input_dat;
    string current_type;

    // 记录每种类型元器件出现的次数
    map<string, int> device_count;
    device_count[string("DRIFT")] = 0;
    device_count[string("QUAD")] = 0;
    device_count[string("BEND")] = 0;
    device_count[string("SOLENOID")] = 0;
    device_count[string("APERTURECIRCULAR")] = 0;
    device_count[string("APERTURERECTANGULAR")] = 0;

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
            dipole_p->SetAngle( stod(Bend_line[2])*RADIAN );
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
            dipole_p->SetKineticEnergy( stod(Bend_line[7]) );
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
        else if (current_type == "APERTURERECTANGULAR")
        {
            ApertureRectangular* apertureRectangular_p = new ApertureRectangular( one_line[0] );
            apertureRectangular_p->SetIn();
            apertureRectangular_p->SetApertureXLeft( stod(one_line[2])/1000.0 );
            apertureRectangular_p->SetApertureXRight( stod(one_line[3])/1000.0 );
            apertureRectangular_p->SetApertureYBottom( stod(one_line[4])/1000.0 );
            apertureRectangular_p->SetApertureYTop( stod(one_line[5])/1000.0 );
            beamline.AddElement(apertureRectangular_p);
        }
        else if (current_type == "APERTURECIRCULAR")
        {
            ApertureCircular* apertureCircular_p = new ApertureCircular( one_line[0] );
            apertureCircular_p->SetIn();
            apertureCircular_p->SetAperture( stod(one_line[2])/1000.0 );
            beamline.AddElement(apertureCircular_p);
        }
    }
}

char* my_simulator::get_Beamline_ElementNames_c(){
    delete rt_names;
    vector<string> beamlinenames = beamline.GetElementNames();
    string names_str = "";
    for(int i=0; i<beamlinenames.size(); i++){
        names_str += beamlinenames[i];
        if(i != beamlinenames.size()-1){
            names_str += ",";
        }
    }
    int charlength = strlen(names_str.c_str()) + 1;
    rt_names = new char[charlength];
    strcpy_s(rt_names, charlength, names_str.c_str());
    return rt_names;
}

char* my_simulator::get_Beamline_ElementTypes_c(){
    delete rt_types;
    int bl_size = beamline.GetSize();
    string types_str = "";
    for(int i=0; i<bl_size; i++){
        types_str += beamline[i]->GetType();
        if(i != bl_size-1){
            types_str += ",";
        }
    }
    int charlength = strlen(types_str.c_str()) + 1;
    rt_types = new char[charlength];
    strcpy_s(rt_types, charlength, types_str.c_str());
    return rt_types;
}

char* my_simulator::get_Beamline_ElementLengths_c(){
    delete rt_lengths;
    int bl_size = beamline.GetSize();
    string lengths_str = "";
    for(int i=0; i<bl_size; i++){
        lengths_str += to_string( beamline[i]->GetLength() );
        if(i != bl_size-1){
            lengths_str += ",";
        }
    }
    int charlength = strlen(lengths_str.c_str()) + 1;
    rt_lengths = new char[charlength];
    strcpy_s(rt_lengths, charlength, lengths_str.c_str());
    return rt_lengths;
}

void my_simulator::init_spacecharge_c(uint r_nr, uint r_nz, int r_adj_bunch){
    spacecharge = new Scheff(r_nr, r_nz, r_adj_bunch);
    spacecharge->SetInterval(0.025);
    spacecharge->SetAdjBunchCutoffW(0.8);
    spacecharge->SetRemeshThreshold(0.02);
}

double **my_simulator::simulate_and_getEnvelope_c(bool use_spacecharge, bool is_sig){
    delete rt_envelope;
    int n = 10240, m = 4;
    rt_envelope = (double **)malloc(10240 * sizeof(double *));
    for(int i=0; i<n; i++){
        rt_envelope[i] = (double *)malloc(m* sizeof(double));
    }
    // double (*rt_envelope)[4] = new double[10240][4];

    SetGPU(0);

    delete simulator;
    simulator = new SimulationEngine();
    PlotData* pltdata = new PlotData(particle_size);
    simulator->InitEngine(beam, &beamline, spacecharge, false, pltdata);
    if(use_spacecharge){
        simulator->SetSpaceCharge(string("on"));
    }else{
        simulator->SetSpaceCharge(string("off"));
    }

    // ofstream out_file(filename.c_str());
    double position = 0.0;
    // vector<double> maxXY;
    // beam->UpdateSigX();
    // beam->UpdateSigY();
    beam->UpdateLoss();
    // maxXY = getBeamMax(beam);
    
    // out_file<< "0.0" <<" "<<beam->GetSigX(true)<< " "<<beam->GetSigY(true)<<" "<<beam->GetLossNum()<<endl;
    if(is_sig){
        rt_envelope[0][0] = 0.0; rt_envelope[0][1] = getBeamSigx_c();
        rt_envelope[0][2] = getBeamSigy_c(); rt_envelope[0][3] = beam->GetLossNum();
    }else{
        rt_envelope[0][0] = 0.0; rt_envelope[0][1] = getBeamMaxx_c();
        rt_envelope[0][2] = getBeamMaxy_c(); rt_envelope[0][3] = beam->GetLossNum();
    }
    
    envelope_size = 1;
    
    vector<string> element_names = beamline.GetElementNames();
    for(string element: element_names){
        BeamLineElement* temp_element = beamline[element];
        position += temp_element->GetLength();
        simulator->Simulate(element, element);
        // beam->UpdateSigX();
        // beam->UpdateSigY();
        beam->UpdateLoss();
        // maxXY = getBeamMax(beam);
        // out_file<< position <<" "<<beam->GetSigX(true)<< " "<<beam->GetSigY(true)<<" "<<beam->GetLossNum()<<endl;
        if(is_sig){
            rt_envelope[envelope_size][0] = position; rt_envelope[envelope_size][1] = getBeamSigx_c();
            rt_envelope[envelope_size][2] = getBeamSigy_c(); rt_envelope[envelope_size][3] = beam->GetLossNum();
        }else{
            rt_envelope[envelope_size][0] = position; rt_envelope[envelope_size][1] = getBeamMaxx_c();
            rt_envelope[envelope_size][2] = getBeamMaxy_c(); rt_envelope[envelope_size][3] = beam->GetLossNum();
        }
        envelope_size++;

        // cout<<element<<" "<<temp_element->GetType()<<endl;
        // beam->UpdateAvgX();
        // cout<<"AvgX:" << beam->GetAvgX()<<endl;
        // beam->UpdateAvgY();
        // cout<<"AvgY:" << beam->GetAvgY()<<endl;
        // cout<<"SigX:" << beam->GetSigX(true)<<endl;
        // cout<<"SigY:" << beam->GetSigY(true)<<endl;
    }

    delete pltdata;
    return rt_envelope;
}

void my_simulator::simulate_c(bool use_spacecharge, bool is_sig){
    SetGPU(0);

    delete simulator;
    simulator = new SimulationEngine();
    PlotData* pltdata = new PlotData(particle_size);
    simulator->InitEngine(beam, &beamline, spacecharge, false, pltdata);
    if(use_spacecharge){
        simulator->SetSpaceCharge(string("on"));
    }else{
        simulator->SetSpaceCharge(string("off"));
    }
    
    vector<string> element_names = beamline.GetElementNames();
    simulator->Simulate(element_names.front(), element_names.back());

    delete pltdata;
}

int my_simulator::get_envelope_size(){
    return envelope_size;
}

void my_simulator::set_magnet_with_index_c(int magnet_index, double field_or_angle){
    BeamLineElement* magnet = beamline[magnet_index];
    string magnet_type = magnet->GetType();
    if(magnet_type == "Dipole"){
        Dipole* dipole = dynamic_cast<Dipole*>(magnet);
        dipole->SetAngle(field_or_angle);
    }else if(magnet_type == "Quad"){
        Quad* quad = dynamic_cast<Quad*>(magnet);
        quad->SetGradient(field_or_angle);
    }else if(magnet_type == "Solenoid"){
        Solenoid* sole = dynamic_cast<Solenoid*>(magnet);
        sole->SetField(field_or_angle);
    }
}

void my_simulator::set_magnet_with_name_c(string elementname, double field_or_angle){
    BeamLineElement* magnet = beamline[elementname];
    string magnet_type = magnet->GetType();
    if(magnet_type == "Dipole"){
        Dipole* dipole = dynamic_cast<Dipole*>(magnet);
        dipole->SetAngle(field_or_angle);
    }else if(magnet_type == "Quad"){
        Quad* quad = dynamic_cast<Quad*>(magnet);
        quad->SetGradient(field_or_angle);
    }else if(magnet_type == "Solenoid"){
        Solenoid* sole = dynamic_cast<Solenoid*>(magnet);
        sole->SetField(field_or_angle);
    }
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

    EXPORT void beam_init_from_file(char* file_path){
        mysim->beam_init_from_file_c(string(file_path));
    }  

    EXPORT void beam_print_to_file(char* file_path, char* comment = "CLAPA"){
        mysim->beam_print_to_file_c(string(file_path), string(comment));
    }  

    EXPORT void set_beamTwiss(double r_ax, double r_bx, double r_ex, 
                double r_ay, double r_by, double r_ey, double r_az, double r_bz, double r_ez,
                double r_sync_phi, double r_sync_w, double r_freq, unsigned int r_seed = 1){
        mysim->set_beamTwiss_c(r_ax, r_bx, r_ex, r_ay, r_by, r_ey, r_az, r_bz, r_ez, r_sync_phi, r_sync_w, r_freq, r_seed);
    }  

    EXPORT void save_initial_beam(){
        mysim->save_initial_beam_c();
    }

    EXPORT void restore_initial_beam(){
        mysim->restore_initial_beam_c();
    }

    EXPORT double getBeamAvgx(){
        return mysim->getBeamAvgx_c();
    }

    EXPORT double getBeamAvgy(){
        return mysim->getBeamAvgy_c();
    }
    
    EXPORT double getBeamSigx(){
        return mysim->getBeamSigx_c();
    }

    EXPORT double getBeamSigy(){
        return mysim->getBeamSigy_c();
    }

    EXPORT double getBeamMaxx(){
        return mysim->getBeamMaxx_c();
    }

    EXPORT double getBeamMaxy(){
        return mysim->getBeamMaxy_c();
    }

    EXPORT void free_beam(){
        mysim->free_beam_c();
    }

    EXPORT void init_database(char* DB_Url){
        // cout<<"start"<<endl;
        // cout<<string(DB_Url)<<endl;
        mysim->init_database_c(string(DB_Url));
    }  

    EXPORT void init_beamline_from_DB(){
        mysim->init_beamline_from_DB_c();
    }  

    EXPORT void load_Beamline_From_DatFile(char* filename){
        mysim->load_Beamline_From_DatFile_c(string(filename));

    }

    EXPORT char* get_Beamline_ElementNames(){
        return mysim->get_Beamline_ElementNames_c();
    }

    EXPORT char* get_Beamline_ElementTypes(){
        return mysim->get_Beamline_ElementTypes_c();
    }

    EXPORT char* get_Beamline_ElementLengths(){
        return mysim->get_Beamline_ElementLengths_c();
    }

    EXPORT void init_spacecharge(uint r_nr = 32, uint r_nz = 128, int r_adj_bunch = 3){
        mysim->init_spacecharge_c(r_nr, r_nz, r_adj_bunch);
    }  

    EXPORT double **simulate_and_getEnvelope(bool use_spacecharge, bool is_sig){
        return mysim->simulate_and_getEnvelope_c(use_spacecharge, is_sig);
    }   

    EXPORT void simulate(bool use_spacecharge, bool is_sig){
        mysim->simulate_c(use_spacecharge, is_sig);
    }

    EXPORT int get_envelope_size(){
        return mysim->get_envelope_size();
    }

    EXPORT void set_magnet_with_index(int magnet_index, double field_or_angle){
        mysim->set_magnet_with_index_c(magnet_index, field_or_angle);
    }

    EXPORT void set_magnet_with_name(char* element_name, double field_or_angle){
        mysim->set_magnet_with_name_c(string(element_name), field_or_angle);
    }

    EXPORT int get_good_number(){
        vector<uint> transverse_loss = mysim->beam->GetLoss();
        int good_count = 0;
        for(uint state: transverse_loss){
            if(state == 0) good_count++;
        }
        return good_count;
    }

}