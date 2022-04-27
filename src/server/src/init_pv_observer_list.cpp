#include <iostream>
#include "init_pv_observer_list.h"

// TODO: move typedef to sql_utility
typedef std::vector<std::vector<std::string> > VV;
typedef std::map<PVObserver*, std::vector<std::string> > Mmap;

void InitPVObserverList(PVObserverList& r_pvlist, BeamLine& r_bl, 
  DBConnection& r_dbcon, bool r_verbose)
{
  std::vector<std::string> dbs = r_dbcon.dbs;
  sqlite3* db_conn = r_dbcon.db_conn;

  for(int dbs_indx = 0; dbs_indx < dbs.size(); ++dbs_indx)
  {
    std::string db = dbs[dbs_indx]; 
    std::string sql = "select lcs_name, value_type from " + db 
		      + ".epics_channel";
    VV pv_list = GetQueryResults(db_conn, sql.c_str());
    // record the masters for the current db, populate them after done with all
    // the other PVs in the current db, not sure if the effect of master PVs 
    // only exists in the current db.
    Mmap masters;
    for(int i = 0; i < pv_list.size(); ++i)
    {
      std::string pv = pv_list[i][0]; 
      std::string pv_type = pv_list[i][1];
      sql = "select model_type, name from " + db + 
	    ".channel_list where channel1='" + pv + "' or channel2='" + 
	    pv + "' or channel3 ='" + pv + "' or channel4 ='" + pv + "'";
      VV elems_info = GetQueryResults(db_conn, sql.c_str());
      for(int j = 0; j < elems_info.size(); ++j)
      {
        std::string elem_type = elems_info[j][0];
        std::string elem_name = elems_info[j][1];
        // first time seeing this pv, create the PVObserver
        // TODO: change this to using find
        if(r_pvlist[pv] == NULL) 
        { 
          if(elem_type == "quad")
            r_pvlist.AddPVObserver(pv, new QuadPVObserver(pv, db));    
          else if(elem_type == "rf_module")
          {
            if(pv_type == "rf_ph")
              r_pvlist.AddPVObserver(pv, new RFPhasePVObserver(pv, db));
            else if(pv_type == "rf_amp" || pv_type == "delay") 
              r_pvlist.AddPVObserver(pv, new RFAmplitudePVObserver(pv, db));
	    else if(pv_type == "rf_ph_master")
	    {
	      MasterPVObserver* amaster = new MasterPVObserver(pv, db);
              r_pvlist.AddPVObserver(pv, amaster);
	      masters[amaster] = std::vector<std::string>();
	    }
          }
          else if(elem_type == "buncher")
          {
            if(pv_type == "buncher_ph") 
              r_pvlist.AddPVObserver(pv, new BuncherPhasePVObserver(pv, db));
            else if(pv_type == "buncher_amp")
              r_pvlist.AddPVObserver(pv, 
				     new BuncherAmplitudePVObserver(pv, db));
            else if(pv_type == "buncher_on_off") 
              r_pvlist.AddPVObserver(pv, new BuncherOnOffPVObserver(pv, db));
          }
          else if(elem_type == "dipole")
            r_pvlist.AddPVObserver(pv, new DipolePVObserver(pv, db));
	  if(r_verbose)
	    std::cout << "---------- Add PV: " << pv << std::endl;
        } // if r_pvlist
        // Attach BeamLineElement to the PVObserver
        if(elem_type != "rf_module" && r_bl[elem_name] != NULL)
        {
            r_pvlist.AttachBeamLineElementToPVObserver(pv, r_bl[elem_name]);
	    if(r_verbose)
	      std::cout << "Attach " << elem_name << " to " << pv << std::endl;
        }
        else if (elem_type == "rf_module")
        {
	  if(pv_type != "rf_ph_master")
	  {
	    sql = "select g.name from " + db + ".rf_gap g join " + db + 
		  ".rf_module m on m.id = g.module_id where m.name = '" + 
		  elem_name + "'";  
	    VV gap_names = GetQueryResults(db_conn, sql.c_str());
	    for(int gp = 0; gp < gap_names.size(); ++gp)
	    {
	      r_pvlist.AttachBeamLineElementToPVObserver(pv, 
						 	r_bl[gap_names[gp][0]]);
	      if(r_verbose)
		std::cout << "Attach " << gap_names[gp][0] << " to " 
		  << pv << std::endl;
	    }
	  }
	  else
	    masters[r_pvlist[pv]].push_back(elem_name);
        }
        else
          std::cerr << "InitPVObserverList error: Cannot find " << elem_name 
	    << " in beamline!" << std::endl;
      }// for j elems_info
    }// for i pv_list
    // Now deal with the master PVs
    Mmap::iterator it = masters.begin();
    for(; it != masters.end(); ++it)
    {
      std::string mpv = it->first->GetPV();
      MasterPVObserver* mob = dynamic_cast<MasterPVObserver*>(it->first);
      std::vector<std::string>& elems = it->second;
      for(int i = 0; i < elems.size(); ++i)
      {
	sql = "select lcs_name from " + db + ".epics_channel e join " + db + 
	      ".rf_module m where m.name = '" + elems[i] + "' and "
	      "m.phase_channel = e.id;";	
	VV pch = GetQueryResults(db_conn, sql.c_str());
	std::string pch_pv= pch[0][0];
	mob->AttachPVObserver(r_pvlist[pch_pv]);
      }
    }
  }// for dbs_indx
  r_pvlist.SetDBconn(db_conn);
}
