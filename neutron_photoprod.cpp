#include <iostream>
#include <cassert>
#include <string>
#include <fstream>
#include <chrono>

#include "TString.h"

#include "gsl/gsl_math.h"

#include "BH_Deu.h"

using namespace std;

//root .x neutron_photoprod.cpp

int neutron_photoprod()
{
	auto starttime = chrono::system_clock::now();
	auto endtime = starttime;
	auto durtime = chrono::duration_cast<chrono::minutes>(endtime - starttime);

	TString filename="/w/work/clas12/tyson/data_repo/BH_deut/neutron_photoprod_10p6/nphoto10p6_";
	int NFiles = 1;
    Long64_t evpfile = 10000,NTotEvents=0;

	// electron: lepType -> 0      muon: lepType -> 1
	int lepType = 0;
    int t_type=2112;
    double m=0.939565;
    int leppid=11;
    if(lepType==1){leppid=13;}

    double Ebeam = 10.6;

	double var[10];
	double weight = 0.0;
	TLorentzVector pgamFv, pelecFV;
	TLorentzVector pout[3];
	TLorentzVector kf[2];

	// initialize the model for Bethe-Heitler process
	BH_deuteron::setModel(lepType, Ebeam,t_type);

	// for Bremsstrahlung Photon
	double kmin = 5.0;
	double kmax = Ebeam;
	incidentPhoton::SetBremsstrahlung();

	

	for (int j = 0; j < NFiles; j++) 
	{
		cout << "FileNb: "<< j << endl;

        FILE* f;
	    TString name;
		name = filename + Form("%.4d.dat", j);
		f = fopen(name.Data(), "w");

        Long64_t physEveNum=0;
		while(physEveNum<evpfile)
		{
			// 15cm LD2 target,   0.01 factor is included in the distribution function
			weight = incidentPhoton::BremsstrahlungPhoton(&pgamFv, kmin, kmax, Ebeam,5.0,769.1);
			weight *= BH_deuteron::GetBHdeu(&pgamFv, pout, var);

			if (weight > 0.0) 
			{
				physEveNum ++;
                NTotEvents++;

                //header
				fprintf(f, "%d %d %d %d %d %d %.6E %d %.6E %.6E\n",4,1,1,0,0,11,Ebeam,t_type,m, weight);
                //e-
				fprintf(f, "%d %d %d %d %d %d %.6E %.6E %.6E %.6E %.6E %.6E %.6E %.6E\n",1, 0, 1,leppid,0,0, pout[1].X(), pout[1].Y(), pout[1].Z(), pout[1].E(),pout[1].M(),0.0,0.0,0.0);
				//e+
				fprintf(f, "%d %d %d %d %d %d %.6E %.6E %.6E %.6E %.6E %.6E %.6E %.6E\n",2, 0, 1,-leppid,0,0, pout[2].X(), pout[2].Y(), pout[2].Z(), pout[2].E(),pout[2].M(),0.0,0.0,0.0);
                //n
				fprintf(f, "%d %d %d %d %d %d %.6E %.6E %.6E %.6E %.6E %.6E %.6E %.6E\n",3, 0, 1,t_type,0,0, pout[0].X(), pout[0].Y(), pout[0].Z(), pout[0].E(),pout[0].M(),0.0,0.0,0.0);
                //gamma -- don't actually want to simulate but could be good info
                //number 0 before pid makes geant4 not propagate it
				fprintf(f, "%d %d %d %d %d %d %.6E %.6E %.6E %.6E %.6E %.6E %.6E %.6E\n",4, 0, 0,22,0,0, pgamFv.X(), pgamFv.Y(), pgamFv.Z(), pgamFv.E(),pgamFv.M(),0.0,0.0,0.0);
                
			}
		}
		fclose(f);
	}
	cout << "Generated "<<NTotEvents<<" events in"<<NFiles<<" files."<< endl;
	endtime = chrono::system_clock::now();
	durtime = chrono::duration_cast<chrono::minutes>(endtime - starttime);
	cout << "Execution time: " << durtime.count() << " minutes " << endl;

	return 0;
}