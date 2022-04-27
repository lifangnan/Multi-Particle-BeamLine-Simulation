#include <cstdlib>
#include <iostream>
#include "pv_observer.h"
#include "sql_utility.h"

/*!
 * \brief Constructor.
 * \param r_pv Name of the EPICS PV
 * \param r_db Name of the database
 */
PVObserver::PVObserver(std::string r_pv, std::string r_db): pv_(r_pv), db_(r_db)
{
}

/*!
 * \brief Update both the database & model for simulation.
 */
void PVObserver::Update(std::string r_val)
{
  val_ = r_val;
  UpdateDB();
  UpdateModel();
}

/*!
 * \brief Update the database.
 */
void PVObserver::UpdateDB() 
{
  char* errmsg;
  sqlite3_exec(db_conn_, "BEGIN TRANSACTION", NULL, NULL, &errmsg);
  std::string sql = "update " + db_ + ".epics_channel set value = " + val_ + 
    " where lcs_name = '" + pv_ + "'";
  sqlite3_stmt* stmt_;
  SQLCheck(sqlite3_prepare_v2(db_conn_, sql.c_str(), -1, &stmt_, NULL), 
    "sqlite3_prepare: " + sql);
  sqlite3_step(stmt_);
  SQLCheck(sqlite3_finalize(stmt_), stmt_, "sqlite3_finalize for PV:" + pv_);
  sqlite3_exec(db_conn_, "END TRANSACTION", NULL, NULL, &errmsg);
}

MasterPVObserver::MasterPVObserver(std::string r_pv, std::string r_db) 
  : PVObserver(r_pv, r_db)
{
}

void MasterPVObserver::AttachBeamLineElement(BeamLineElement* r_elem)
{
  std::cerr << "Cannot directly attach element to MasterPVObserver! " 
    << std::endl;
}

void MasterPVObserver::AttachPVObserver(PVObserver* r_pvo)
{
  pvo_.push_back(r_pvo);
}

void MasterPVObserver::UpdateModel()
{
  for(int i = 0; i < pvo_.size(); ++i)
    pvo_[i]->UpdateModel(); 
}

/*!
 * \brief Get the names of the observed beamline elements.
 */
std::vector<std::string> MasterPVObserver::GetBeamLineElementNames() const
{
  std::vector<std::string> rlt(pvo_.size(), "");
  for(int i = 0; i < pvo_.size(); ++i)
    rlt[i] = pvo_[i]->GetPV();
  return rlt;
}

QuadPVObserver::QuadPVObserver(std::string r_pv, std::string r_db) 
  : PVObserver(r_pv, r_db)
{
}

void QuadPVObserver::AttachBeamLineElement(BeamLineElement* r_elem)
{
  if (Quad* tmp_quad = dynamic_cast<Quad*>(r_elem))
    quad_.push_back(tmp_quad);
  else
  {
    std::cerr << "Cann't attach " << r_elem->GetName() 
      << " to a QuadPVObserver!" << std::endl;
    exit(-1);
  }
}

void QuadPVObserver::UpdateModel()
{
//  char* errmsg;  
//  sqlite3_exec(GetDBconn(), "BEGIN TRANSACTION", NULL, NULL, &errmsg);

  for(int i = 0; i < quad_.size(); ++i)
  {
    std::string sql = "select gradient_model from " + GetDB() + 
      ".quad where name = '" + quad_[i]->GetName() + "'";
    std::string data = GetDataFromDB(GetDBconn(), sql.c_str());
    if (data != "")
      quad_[i]->SetGradient(std::atof(data.c_str()));
    else
     std::cerr << "QuadPVObserver::UpdateModel() failed, no "
	"gradient were found for " << "quad : " << quad_[i]->GetName() 
	<< std::endl; 
  }
//sqlite3_exec(GetDBconn(), "END TRANSACTION", NULL, NULL, &errmsg);
}

std::vector<std::string> QuadPVObserver::GetBeamLineElementNames() const
{
  std::vector<std::string> rlt(quad_.size(), "");
  for(int i = 0; i < quad_.size(); ++i)
    rlt[i] = quad_[i]->GetName(); 
  return rlt; 
}

RFPhasePVObserver::RFPhasePVObserver(std::string r_pv, std::string r_db)
  : PVObserver(r_pv, r_db)
{
}

void RFPhasePVObserver::UpdateModel()
{
  for(int i = 0; i < gap_.size(); ++i)
  {
    std::string sql = "select beam_phase_shift_model from " + GetDB() + 
      ".rf_gap where name = '" + gap_[i]->GetName() + "'";
    std::string data = GetDataFromDB(GetDBconn(), sql.c_str());
    if(data != "")
      gap_[i]->SetPhaseShift(std::atof(data.c_str()));
    else
      std::cerr << "RFPhasePVObserver::UpdateModel() failed, no "
	"beam_phase_shift were found for RFGap: " << gap_[i]->GetName() 
	<< std::endl;
  }
}

void RFPhasePVObserver::AttachBeamLineElement(BeamLineElement* r_elem)
{
  if (RFGap* tmp_gap = dynamic_cast<RFGap*>(r_elem))
    gap_.push_back(tmp_gap);
  else
  {
    std::cerr << "Cann't attach " << r_elem->GetName() 
      << " to a RFPhasePVObserver!" << std::endl;
    exit(-1);
  }
}

std::vector<std::string> RFPhasePVObserver::GetBeamLineElementNames() const
{
  std::vector<std::string> rlt(gap_.size(), "");
  for(int i = 0; i < gap_.size(); ++i)
    rlt[i] = gap_[i]->GetName(); 
  return rlt; 
}

RFAmplitudePVObserver::RFAmplitudePVObserver(std::string r_pv, std::string r_db)
  : PVObserver(r_pv, r_db)
{
}

void RFAmplitudePVObserver::UpdateModel()
{
  for(int i = 0; i < gap_.size(); ++i)
  {
    std::string sql = "select amplitude_model, ref_phase_model from " 
                + GetDB() + ".rf_gap where name = '" + gap_[i]->GetName() + "'";
    std::vector<std::vector<std::string> > data = GetQueryResults(GetDBconn(), 
      sql.c_str());
    if(!data.empty())
    {
      gap_[i]->SetRFAmplitude(std::atof(data[0][0].c_str()));
      gap_[i]->SetRefPhase(std::atof(data[0][1].c_str()));
    }
    else
      std::cerr << "RFAmplitudePVObserver::UpdateModel() failed, no amplitude "
	"& ref_phase were found for RFGap: " << gap_[i]->GetName() << std::endl;
  }
}

void RFAmplitudePVObserver::AttachBeamLineElement(BeamLineElement* r_elem)
{
  if (RFGap* tmp_gap = dynamic_cast<RFGap*>(r_elem))
    gap_.push_back(tmp_gap);
  else
  {
    std::cerr << "Cann't attach " << r_elem->GetName() 
      << " to a RFAmplitudePVObserver!" << std::endl;
    exit(-1);
  }
}

std::vector<std::string> RFAmplitudePVObserver::GetBeamLineElementNames() const
{
  std::vector<std::string> rlt(gap_.size(), "");
  for(int i = 0; i < gap_.size(); ++i)
    rlt[i] = gap_[i]->GetName(); 
  return rlt; 
}

BuncherPVObserver::BuncherPVObserver(std::string r_pv, std::string r_db)
  : PVObserver(r_pv, r_db)
{
}

void BuncherPVObserver::AttachBeamLineElement(BeamLineElement* r_elem)
{
  if(Buncher* buncher = dynamic_cast<Buncher*>(r_elem))
    buncher_ = buncher;
  else
  {
    std::cerr << "Cann't attach " << r_elem->GetName() 
      << " to a BuncherPhasePVObserver!" << std::endl;
    exit(-1);
  }
}

std::vector<std::string> BuncherPVObserver::GetBeamLineElementNames() const
{
  std::vector<std::string> rlt(1, buncher_->GetName());
  return rlt; 
}

BuncherPhasePVObserver::BuncherPhasePVObserver(std::string r_pv, 
  std::string r_db) : BuncherPVObserver(r_pv, r_db)
{
}

void BuncherPhasePVObserver::UpdateModel()
{
  std::string sql = "select phase_model from " + GetDB() + 
    ".buncher where name = '" + buncher_->GetName() + "'";
  std::string data = GetDataFromDB(GetDBconn(), sql.c_str());
  if (data != "")
    buncher_->SetPhase(std::atof(data.c_str()));
  else
    std::cerr << "BuncherPhasePVObserver::UpdateModel() failed, no phase_model "
    "were found for buncher : " << buncher_->GetName() << std::endl; 
}

BuncherAmplitudePVObserver::BuncherAmplitudePVObserver(std::string r_pv, 
  std::string r_db) : BuncherPVObserver(r_pv, r_db)
{
}

void BuncherAmplitudePVObserver::UpdateModel()
{
  std::string sql = "select voltage_model from " + GetDB() + ".buncher where "
    "name = '" + buncher_->GetName() + "'";
  std::string data = GetDataFromDB(GetDBconn(), sql.c_str());
  if (data != "")
    buncher_->SetVoltage(std::atof(data.c_str()));
  else
    std::cerr << "BuncherAmplitudePVObserver::UpdateModel() failed, no "
    "voltage_model were found for buncher : " << buncher_->GetName() 
    << std::endl; 
}

BuncherOnOffPVObserver::BuncherOnOffPVObserver(std::string r_pv, 
  std::string r_db) : BuncherPVObserver(r_pv, r_db)
{
}

void BuncherOnOffPVObserver::UpdateModel()
{
  std::string sql = "select on_off from " + GetDB() + ".buncher where name = '" 
                    + buncher_->GetName() + "'";
  std::string data = GetDataFromDB(GetDBconn(), sql.c_str());
  if (data != "")
    if(std::atof(data.c_str()) > 0)
      buncher_->TurnOn();
    else
      buncher_->TurnOff();
  else
   std::cerr << "BuncherOneOffPVObserver::UpdateModel() failed, no on_off were "
    "found for buncher : " << buncher_->GetName() << std::endl; 
}

DipolePVObserver::DipolePVObserver(std::string r_pv, std::string r_db)
  : PVObserver(r_pv, r_db)
{
}

void DipolePVObserver::AttachBeamLineElement(BeamLineElement* r_elem)
{
  if(Dipole* dipole = dynamic_cast<Dipole*>(r_elem))
    dipole_.push_back(dipole);
  else if(ApertureRectangular* aper = dynamic_cast<ApertureRectangular*>(r_elem))
    aperture_r_.push_back(aper);
  else if(Drift* drift = dynamic_cast<Drift*>(r_elem))
    drift_.push_back(drift);
  else
  {
    std::cerr << "Cann't attach " << r_elem->GetName() 
      << " to a DipolePhasePVObserver!" << std::endl;
    exit(-1);
  }
}

std::vector<std::string> DipolePVObserver::GetBeamLineElementNames() const
{
  std::vector<std::string> rlt(dipole_.size(), "");
  for(int i = 0; i < dipole_.size(); ++i)
    rlt[i] = dipole_[i]->GetName(); 
  for(int i = 0; i < aperture_r_.size(); ++i)
    rlt[i] = aperture_r_[i]->GetName();
  for(int i = 0; i < drift_.size(); ++i)
    rlt[i] = drift_[i]->GetName();
  return rlt; 
}

void DipolePVObserver::UpdateModel()
{
  for(int i = 0; i < dipole_.size(); ++i)
  {
    std::string sql = "select rho_model, angle_model, edge_angle1_model, "
      "edge_angle2_model, kenergy_model from " + GetDB() + 
      ".dipole where name = '" + dipole_[i]->GetName() + "'";
    std::vector<std::vector<std::string> > data = GetQueryResults(GetDBconn(), 
      sql.c_str());
    if(!data.empty())
    {
      dipole_[i]->SetRadius(std::atof(data[0][0].c_str()));
      dipole_[i]->SetAngle(std::atof(data[0][1].c_str()));
      dipole_[i]->SetEdgeAngleIn(std::atof(data[0][2].c_str()));
      dipole_[i]->SetEdgeAngleOut(std::atof(data[0][3].c_str()));
      dipole_[i]->SetKineticEnergy(std::atof(data[0][4].c_str()));
    }
    else
      std::cerr << "DipolePVObserver::UpdateModel() failed, "
	"for dipole : " << dipole_[i]->GetName() << std::endl;
  }
  for(int i = 0; i < aperture_r_.size(); ++i)
  {
    std::string sql = "select aperture_xl_model, aperture_xr_model from " + 
       GetDB() + ".raperture where name = '" + aperture_r_[i]->GetName() + "'";
    std::vector<std::vector<std::string> > data = GetQueryResults(GetDBconn(), 
      sql.c_str());
    if(!data.empty())
    {
      aperture_r_[i]->SetApertureXLeft(std::atof(data[0][0].c_str()));
      aperture_r_[i]->SetApertureXRight(std::atof(data[0][1].c_str()));
    }
    else
      std::cerr << "DipolePVObserver::UpdateModel() failed, "
	"for rectangular aperture : " << aperture_r_[i]->GetName()<< std::endl;
  }
  for(int i = 0; i < drift_.size(); ++i)
  {
    std::string sql = "select length_model from " + GetDB() + 
      ".drift where name = '" +drift_[i]->GetName() + "'";
    std::vector<std::vector<std::string> > data = GetQueryResults(GetDBconn(), 
      sql.c_str());
    if(!data.empty())
      drift_[i]->SetLength(std::atof(data[0][0].c_str()));
    else
      std::cerr << "DipolePVObserver::UpdateModel() failed, "
	"for drift : " << drift_[i]->GetName()<< std::endl;
  }
}
