#include <iostream>
#include <cassert>
#include <string>
#include <fstream>
#include <chrono>

#include "TString.h"

#include "gsl/gsl_math.h"

#include "BH_Deu.h"

using namespace std;

//root .x proton_electroprod.cpp

int proton_electroprod()
{
	auto starttime = chrono::system_clock::now();
	auto endtime = starttime;
	auto durtime = chrono::duration_cast<chrono::minutes>(endtime - starttime);

	//TString filename="/w/work/clas12/tyson/data_repo/BH_deut/proton_electroprod_10p6/pelec10p6_";
    TString filename="/w/work/clas12/tyson/data_repo/BH_deut/proton_electroprod_10p2/pelec10p2_";
	int NFiles = 100;
    Long64_t evpfile = 10000,NTotEvents=0;
   //49 mins for 10M events
   //5 min for 1M events

	// electron: lepType -> 0      muon: lepType -> 1
	int lepType = 0;
    int t_type=2212;
    double m=0.938272;

    int leppid=11;
    if(lepType==1){leppid=13;}

    double Ebeam = 10.2;

	double var[10];
	double weight = 0.0;
	TLorentzVector pgamFv, pelecFV;
	TLorentzVector pout[3];
	TLorentzVector kf[2];

	// initialize the model for Bethe-Heitler process
	BH_deuteron::setModel(lepType, Ebeam,t_type);

	// for electroproduction
	incidentPhoton::perange[0] = 0.0;
	incidentPhoton::perange[1] = Ebeam - 5.0;
    //applying a range of cos theta leads to more efficient sampling
    //cross section drops off with polar angle
    double degtorad = M_PI / 180.0;
    incidentPhoton::cthrange[0] = cos(10.0 * degtorad);
	incidentPhoton::cthrange[1] = cos(0.0 * degtorad);

	pelecFV.SetXYZM(0., 0., Ebeam, leptonTensor::melec);

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

            weight = incidentPhoton::VirtualPhoton(&pelecFV, kf);
			pgamFv.SetPxPyPzE(kf[1].Px(), kf[1].Py(), kf[1].Pz(), kf[1].E());
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
                //p
				fprintf(f, "%d %d %d %d %d %d %.6E %.6E %.6E %.6E %.6E %.6E %.6E %.6E\n",3, 0, 1,t_type,0,0, pout[0].X(), pout[0].Y(), pout[0].Z(), pout[0].E(),pout[0].M(),0.0,0.0,0.0);
                //e'
				fprintf(f, "%d %d %d %d %d %d %.6E %.6E %.6E %.6E %.6E %.6E %.6E %.6E\n",4, 0, 1,11,0,0, kf[0].X(), kf[0].Y(), kf[0].Z(), kf[0].E(),kf[0].M(),0.0,0.0,0.0);
                
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