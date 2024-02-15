#pragma once

#include <iostream>
#include <cassert>
#include <string>
#include <fstream>

#include "gsl/gsl_math.h"
#include "gsl/gsl_spline.h"

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"

using namespace std;

namespace formFun
{
    //a struct to define constant parameter
    struct paraconst
    {
        const double GMp = 2.79284356;
        const double GMn = -1.91304272;
        const double proton_mass = 0.938272;                 //GeV
        const double neutron_mass = 0.93956542;
        const double promass = 0.9389187;                  //   (proton mass + neutron mass )/2
        const double deumass = 1.875613;
        const double mJPmass = 3.096916;
        const double alp = 1.0 / 137.036;
    };

    // nucleon form factor from  DOI: 10.1016/j.physletb.2017.11.023  Zhihong Ye et al.
    double getGform(const int kID, const double kQ2)
    {
        double GNGD_Fit;
        // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
        if (kID < 1 || kID>4) {
            cerr << "*** ERROR***, kID is not any of [1->GEp, 2->GMp, 3->GEn, 4->GMn]" << endl;
            GNGD_Fit = -1000.0;
            return GNGD_Fit;
        }

        // z-Expansion Parameters for Form Factor Values

        const double GN_Coef_Fit[4][13] = {
            {0.239163298067, -1.10985857441, 1.44438081306, 0.479569465603, -2.28689474187,  1.12663298498,
            1.25061984354,-3.63102047159, 4.08221702379,  0.504097346499,  -5.08512046051,  3.96774254395,-0.981529071103}, /*GEp*/
            {0.264142994136, -1.09530612212, 1.21855378178, 0.661136493537, -1.40567892503, -1.35641843888,
            1.44702915534, 4.2356697359, -5.33404565341, -2.91630052096,    8.70740306757, -5.70699994375, 1.28081437589}, /*GMp*/
            {0.048919981379,-0.064525053912,-0.240825897382,0.392108744873, 0.300445258602,-0.661888687179,
            -0.175639769687, 0.624691724461,-0.077684299367,-0.236003975259, 0.090401973470, 0.0, 0.0}, /*GEn*/
            {0.257758326959,-1.079540642058, 1.182183812195,0.711015085833,-1.348080936796,-1.662444025208,
            2.624354426029, 1.751234494568,-4.922300878888, 3.197892727312,-0.712072389946, 0.0, 0.0} /*GMn*/
        };


        // Apply the z-expansion formula for form factor
        const double tcut = 0.0779191396;
        const double t0 = -0.7;
        double z = (sqrt(tcut + kQ2) - sqrt(tcut - t0)) / (sqrt(tcut + kQ2) + sqrt(tcut - t0));
        double GNQ2 = 0.0;
        for (int i = 0; i < 13; i++)
        {
            GNQ2 += GN_Coef_Fit[kID - 1][i] * pow(z, i);
        }

        GNGD_Fit = GNQ2;

        //double GDip= pow(1./(1. + kQ2/0.71), 2);
        //GNGD_Fit[0] = GNQ2 / GDip;       //Note that Coef_Fit have been divided by mu_p or mu_n

        return GNGD_Fit;
    }

    double proF1(const double Q2)
    {
        paraconst Gpara;

        double Gmp = Gpara.GMp;
        double mproton = Gpara.proton_mass;

        const int idE = 1;
        const int idM = 2;
        double f1;

        //Note that Coef_Fit have been divided by mu_p or mu_n
        f1 = (getGform(idE, Q2) + Q2 / (4.0 * pow(mproton, 2)) * getGform(idM, Q2) * Gmp) / (1 + Q2 / (4.0 * pow(mproton, 2)));
        return f1;
    }

    double proF2(const double Q2)
    {
        paraconst Gpara;

        double Gmp = Gpara.GMp;
        double mproton = Gpara.proton_mass;

        const int idE = 1;
        const int idM = 2;
        double f2;

        //Note that Coef_Fit have been divided by mu_p or mu_n
        f2 = (getGform(idM, Q2) * Gmp - getGform(idE, Q2)) / (1 + Q2 / (4.0 * pow(mproton, 2)));
        return f2;
    }

    double neuF1(const double Q2)
    {
        paraconst Gpara;

        double Gmn = Gpara.GMn;
        double mneutron = Gpara.neutron_mass;

        const int idE = 3;
        const int idM = 4;
        double f1;

        f1 = (getGform(idE, Q2) + Q2 / (4.0 * pow(mneutron, 2)) * getGform(idM, Q2) * Gmn) / (1 + Q2 / (4.0 * pow(mneutron, 2)));
        return f1;
    }

    double neuF2(const double Q2)
    {
        paraconst Gpara;

        double Gmn = Gpara.GMn;
        double mneutron = Gpara.neutron_mass;

        const int idE = 3;
        const int idM = 4;
        double f2;

        f2 = (getGform(idM, Q2) * Gmn - getGform(idE, Q2)) / (1 + Q2 / (4.0 * pow(mneutron, 2)));
        return f2;
    }


    class deuVertexFun
    {
    public:
        //formClass();
        void forminterpInit();
        void forminterpFree();

        double deufd(const double p);
        double deugd(const double p);
        double deuhd(const double p);
        double deuid(const double p);

    private:
        // directory of the data table for deuteron vertex function
        string datafile = "formtable/";
        enum { wfdatanum = 61, Tlabnum = 57, deuformdatanum = 1001 };

        const gsl_interp_type* usetype = gsl_interp_cspline;
        gsl_interp_accel* deuformaccel[4];
        gsl_spline* deuforminterp[4];

        double deuformtab[5][deuformdatanum];
    };


	void deuVertexFun::forminterpInit()
	{
		//input data table for deuteron vertex function
		ifstream deuforminput(datafile + "deuformdata_WJC_2.dat");
		assert(deuforminput.is_open());

		for (int i = 0; i < deuformdatanum; i++)
		{
			for (int j = 0; j < 5; j++)
			{
				deuforminput >> deuformtab[j][i];
			}
		}
		deuforminput.close();
		cout << "read deuteron vertex function data" << endl;

		for (int i = 0; i < 4; i++)
		{
			deuformaccel[i] = gsl_interp_accel_alloc();
			deuforminterp[i] = gsl_spline_alloc(usetype, deuformdatanum);
			gsl_spline_init(deuforminterp[i], deuformtab[0], deuformtab[i + 1], deuformdatanum);
		}
		cout << "initialize the interpolation function of deuteron vertex form factor " << endl;
	}

    void deuVertexFun::forminterpFree()
    {
        for (int i = 0; i < 4; i++)
        {
            gsl_interp_accel_free(deuformaccel[i]);
            gsl_spline_free(deuforminterp[i]);
        }
    }

    double deuVertexFun::deufd(const double p)
    {
        return gsl_spline_eval(deuforminterp[0], p, deuformaccel[0]);
    }

    double deuVertexFun::deugd(const double p)
    {
        return gsl_spline_eval(deuforminterp[1], p, deuformaccel[1]);
    }

    double deuVertexFun::deuhd(const double p)
    {
        return gsl_spline_eval(deuforminterp[2], p, deuformaccel[2]);
    }

    double deuVertexFun::deuid(const double p)
    {
        return gsl_spline_eval(deuforminterp[3], p, deuformaccel[3]);
    }

}

namespace leptonTensor
{
    const double melec = 0.511e-3;                // electron mass
    const double mmuon = 0.105658;                // muon mass
    const double mdeuteron = 1.875613;
    const double alp = 1.0 / 137.036;

    double mlepton = melec;
    //double thetae, phie;
    //double sll, spn, tqsq;                 // t=q^2   q:  virtual photon momentum
    //double Egam, pgam;

    // kFv: incident photon,  p3Fv: electron,  p4Fv: positron
    double kFv[4];
    double p3Fv[4];
    double p4Fv[4];
    double lepmomFv[4][4];
    double lepcoe[4];

    //four moemtum dot product
    inline double Fvprod(double mom1[4], double mom2[4])
    {
        return mom1[0] * mom2[0] - mom1[1] * mom2[1] - mom1[2] * mom2[2] - mom1[3] * mom2[3];
    }

    void setlepmomVar(double Egama, double pgama, double thetaea, double phiea, double slla, double spna, double tqsqa)
    {
        double pcnorm;
        double p34norm, theta34, alphaq;
        // Lorentz gamma     gamma*beta
        double Lorgam, Lorgamb;

        // nu and qnorm defined in lepton tensor part
        double nul, qnorml;

        double ull, tll;
        double lcoe33p, lcoe44p, lcoek3, lcoek4;
        double lcoegmn, lcoekk, lcoe33, lcoe44;

        pcnorm = sqrt(slla / 4. - pow(mlepton, 2));
        nul = (spna - tqsqa - pow(mdeuteron, 2)) / (2. * mdeuteron);
        qnorml = sqrt(pow(nul, 2) - tqsqa);

        p34norm = sqrt(pow(Egama - nul, 2) - slla);
        // polar angle between p34 and p_gam   p34 = p3 + p4
        theta34 = acos((pow(pgama, 2) + pow(p34norm, 2) - pow(qnorml, 2)) / (2. * pgama * p34norm));
        // polar angle between q and p_gam
        alphaq = acos((pow(pgama, 2) + pow(qnorml, 2) - pow(p34norm, 2)) / (2. * pgama * qnorml));

        //Lorentz transformation
        Lorgam = sqrt(1. + pow(p34norm, 2) / slla);
        Lorgamb = p34norm / sqrt(slla);

        kFv[0] = Egama;
        kFv[1] = pgama * sin(alphaq);
        kFv[2] = 0.;
        kFv[3] = pgama * cos(alphaq);

        // Lorentz transformation and rotation with respect to y axis
        // transform to the nuclear system laboratory frame with q along z direction
        p3Fv[0] = Lorgam * sqrt(pow(mlepton, 2) + pow(pcnorm, 2)) + Lorgamb * pcnorm * cos(thetaea);
        p3Fv[1] = Lorgamb * sqrt(pow(mlepton, 2) + pow(pcnorm, 2)) * sin(theta34 + alphaq) + Lorgam * pcnorm * cos(thetaea) * sin(theta34 + alphaq)
            + pcnorm * cos(phiea) * cos(theta34 + alphaq) * sin(thetaea);
        p3Fv[2] = pcnorm * sin(phiea) * sin(thetaea);
        p3Fv[3] = Lorgamb * sqrt(pow(mlepton, 2) + pow(pcnorm, 2)) * cos(theta34 + alphaq) + Lorgam * pcnorm * cos(theta34 + alphaq) * cos(thetaea)
            - pcnorm * cos(phiea) * sin(theta34 + alphaq) * sin(thetaea);

        p4Fv[0] = Lorgam * sqrt(pow(mlepton, 2) + pow(pcnorm, 2)) - Lorgamb * pcnorm * cos(thetaea);
        p4Fv[1] = Lorgamb * sqrt(pow(mlepton, 2) + pow(pcnorm, 2)) * sin(theta34 + alphaq) - Lorgam * pcnorm * cos(thetaea) * sin(theta34 + alphaq)
            - pcnorm * cos(phiea) * cos(theta34 + alphaq) * sin(thetaea);
        p4Fv[2] = -pcnorm * sin(phiea) * sin(thetaea);
        p4Fv[3] = Lorgamb * sqrt(pow(mlepton, 2) + pow(pcnorm, 2)) * cos(theta34 + alphaq) - Lorgam * pcnorm * cos(theta34 + alphaq) * cos(thetaea)
            + pcnorm * cos(phiea) * sin(theta34 + alphaq) * sin(thetaea);

        for (int i = 0; i < 4; i++)
        {
            lepmomFv[0][i] = 0.;
            lepmomFv[1][i] = kFv[i];
            lepmomFv[2][i] = p3Fv[i];
            lepmomFv[3][i] = p4Fv[i];
        }

        // approximation:
        // for quasi-real scattering, we replace k by q*, but keep ampliutde for real photon scattering unchanged
        // tll=(k-p3)^2     ull=(k-p4)^2
        tll = pow(mlepton, 2) - 2. * Fvprod(kFv, p3Fv) + pow(Egama, 2) - pow(pgama, 2);
        ull = pow(mlepton, 2) - 2. * Fvprod(kFv, p4Fv) + pow(Egama, 2) - pow(pgama, 2);

        lcoe33p = 8. / (pow(mlepton, 2) - tll);
        lcoe44p = 8. / (pow(mlepton, 2) - ull);
        lcoek3 = 4. * (5. * pow(mlepton, 4) + (slla + tll) * ull - pow(mlepton, 2) * (slla + 3. * (tll + ull)))
            / ((pow(mlepton, 2) - tll) * pow(pow(mlepton, 2) - ull, 2));
        lcoek4 = 4. * (5. * pow(mlepton, 4) + (slla + ull) * tll - pow(mlepton, 2) * (slla + 3. * (tll + ull)))
            / ((pow(mlepton, 2) - ull) * pow(pow(mlepton, 2) - tll, 2));

        lcoe33 = lcoe33p + 2. * lcoek3;
        lcoe44 = lcoe44p + 2. * lcoek4;
        lcoekk = -16. * pow(mlepton, 2) / ((pow(mlepton, 2) - tll) * (pow(mlepton, 2) - ull));
        lcoegmn = -2. * (18. * pow(mlepton, 8) - 4. * pow(mlepton, 6) * (3. * slla + 7 * (tll + ull)) +
            tll * ull * (2. * pow(slla, 2) + pow(tll, 2) + pow(ull, 2) + 2. * slla * (tll + ull)) +
            pow(mlepton, 4) * (2. * pow(slla, 2) + 14. * slla * (tll + ull) + 15. * pow(tll + ull, 2)) -
            pow(mlepton, 2) * (2. * pow(slla, 2) * (tll + ull) + 3. * pow(tll + ull, 3) + 4. * slla * (pow(tll, 2) + 3. * tll * ull + pow(ull, 2))))
            / (pow(pow(mlepton, 2) - tll, 2) * pow(pow(mlepton, 2) - ull, 2));

        lepcoe[0] = lcoegmn;
        lepcoe[1] = lcoekk;
        lepcoe[2] = lcoe33;
        lepcoe[3] = lcoe44;

    }

}

namespace nuclearTensor
{
    double mpro = 0.938272; //default proton                  
    const double mdeu = 1.875613;
    const double alp = 1.0 / 137.036;
    const double plim = 2.0; //GeV    the up limit of three internal integral momentum
    int p_type=2212;  //2212 for protons, 2112 for neutrons                          

    const string strpnCoe[3] = { "pwpw_coetab/","pwfsi_coetab/","fsifsi_coetab/" };
    const string strpnCoedim[3] = { "pwpw_coedim/","pwfsi_coedim/","fsifsi_coedim/" };

    const string strtenName[3] = { "H1_","Hc","H2" };
    const string strdimfileExt = "dim.txt";
    const string strfileExt = ".txt";
    const string strfilegmnExt = "gmn.txt";

    const string strdeuForm[4] = { "fd","gd","hd","id" };
    const string strdeuFormc[4] = { "fdc","gdc","hdc","idc" };

    const int npp = 4;
    const int momtenpwNum = 6;
    const string strmompw[npp - 1] = { "P0","p2","q" };

    // dimension matrix for coefficients exponent matrix for a specific deuteron form factor product
    int pwCoeDimmat[4][4][npp][npp];

    // coefficients exponent matrix for a specific deuteron form factor product and a specific momentum tensor
    // pointer array for two dimentional matrix of coefficient exponent table
    // the dimensions for exponent vector matrix are from   Dimmat  , so declare dynamical allocated matrix
    // momtennum plus 1  for g_mn tensor
    int** pwCoeVecmat[4][4][momtenpwNum + 1];

    //coefficient matrix saving actual double value coefficient   polynomials of momentum four product
    double pwCoemat[4][4][momtenpwNum + 1];

    formFun::deuVertexFun* myform = new formFun::deuVertexFun;

    //use qnorm, nu, and p1norm to determine theta with energy conservation 
    //p1norm : the norm of p1 three momentum
    double p1norm, qnorm, nu;
    double theta1, phi1;
    double p2norm;

    //three four momentum for plane wave contribution
    double p1Fv[4];
    double P0Fv[4];
    double p2Fv[4];
    double qFv[4];

    double pwmomFv[npp][4];
    // -2  p2dotp2   P0dotP0    +4  m md fo1 fo2
    double pwMomdotvar[momtenpwNum - 2 + 4];

    void setPType(int in_p_type){
        p_type=in_p_type;
        if(p_type==2112){
            mpro=0.939565;
        }

    }

    int indexFun(int row, int col, int momnum)
    {
        return (row - 1) * (momnum - 1 + (momnum - 1 - (row - 1) + 1)) / 2 + col - (row - 1);
    }

    int** allocateMatrix(int row, int col)
    {
        int** matrix;
        matrix = new int* [row];

        for (int i = 0; i < row; i++)
        {
            matrix[i] = new int[col];
        }

        return matrix;
    }
    void freematrix(int row, int** matrix)
    {
        for (int i = 0; i < row; i++)
        {
            delete[] matrix[i];
        }

        delete[] matrix;
    }

    double polyval(int row, int col, double var[], int** vecmat, int coemomdim, int varprenum)
    {
        // varprenum = 4 for proton-proton term     varprenum = 6 for proton-neutron cross term
        int pdotpnum;

        double coevalmono = 1.0;
        double coeval = 0.0;

        for (int ii = 0; ii < row; ii++)
        {
            coevalmono = 1.0;
            for (int jj = 0; jj < varprenum; jj++)
            {
                coevalmono = coevalmono * pow(var[jj], vecmat[ii][jj]);
            }
            //the num of momentum four product is related to the dimension of mass exponent in the denominator
            pdotpnum = coemomdim - (vecmat[ii][0] + vecmat[ii][1]) / 2;
            for (int jj = varprenum; jj < varprenum + pdotpnum; jj++)
            {
                //minus 1 for the index starts from 0 in c++, 
                //the momdotvar position vector is imported from mathematica, index starts from 1
                coevalmono = coevalmono * var[vecmat[ii][jj] - 1];
            }
            coevalmono = coevalmono * (double)(vecmat[ii][col - 1]);

            coeval = coeval + coevalmono;
        }
        return coeval;
    }

    //contraction function for lepton pair production
    void contractFun(double lcoe[4], double lmom[4][4], double nuclmom[][4], int nuclmomNum, double contrmat[])
    {
        //calculate the contraction of lepton tensor and nuclear tensor
        double contrelem;
        int index;

        //contraction between nuclear part gmn and lepton tensor
        //nu and qnorm must be setted before the function call
        contrelem = 0.;
        contrelem = contrelem + lcoe[0] * ((1. - pow(nu / qnorm, 2)) * 1. + 1. + 1.);

        for (int i = 1; i < 4; i++)
        {
            contrelem = contrelem + lcoe[i]
                * (pow(lmom[i][0] - lmom[i][3] * nu / qnorm, 2) * 1.
                    + pow(lmom[i][1], 2) * (-1.)
                    + pow(lmom[i][2], 2) * (-1.));
        }

        contrmat[0] = contrelem;

        //contraction between nuclear part ( p_mu p_nu + p_nu p_mu )/2  and lepton tensor
        for (int k = 1; k < nuclmomNum; k++)
        {
            for (int n = k; n < nuclmomNum; n++)
            {
                //indexFun gives index from 1
                index = indexFun(k, n, nuclmomNum);
                contrelem = 0.;

                contrelem = contrelem + lcoe[0] * (
                    (1. - pow(nu / qnorm, 2)) * nuclmom[k][0] * nuclmom[n][0] - nuclmom[k][1] * nuclmom[n][1] - nuclmom[k][2] * nuclmom[n][2]);

                for (int i = 1; i < 4; i++)
                {
                    contrelem = contrelem + lcoe[i] * (
                        pow(lmom[i][0] - lmom[i][3] * nu / qnorm, 2) * nuclmom[k][0] * nuclmom[n][0]
                        + pow(lmom[i][1], 2) * nuclmom[k][1] * nuclmom[n][1]
                        + pow(lmom[i][2], 2) * nuclmom[k][2] * nuclmom[n][2]
                        - 2. * lmom[i][1] * (lmom[i][0] - lmom[i][3] * nu / qnorm)
                        * (nuclmom[k][0] * nuclmom[n][1] + nuclmom[n][0] * nuclmom[k][1]) / 2.
                        - 2. * lmom[i][2] * (lmom[i][0] - lmom[i][3] * nu / qnorm)
                        * (nuclmom[k][0] * nuclmom[n][2] + nuclmom[n][0] * nuclmom[k][2]) / 2.
                        + 2. * lmom[i][1] * lmom[i][2]
                        * (nuclmom[k][1] * nuclmom[n][2] + nuclmom[n][1] * nuclmom[k][2]) / 2.);
                }

                contrmat[index] = contrelem;

            }
        }
    }

    void setnuclmomVar(double p1a, double phi1a, double spna, double tqsqa)
    {
        //use p1 to determine proton scattering angle theta1, 
        //for a certain theta1, there may be two solutions for proton momentum p1

        p1norm = p1a;
        phi1 = phi1a;
        nu = (spna - tqsqa - pow(mdeu, 2)) / (2. * mdeu);
        qnorm = sqrt(pow(nu, 2) - tqsqa);

        /*nu = nua;
        qnorm = qa;*/

        theta1 = acos((pow(qnorm, 2) + 2. * (nu + mdeu) * sqrt(pow(mpro, 2) + pow(p1norm, 2)) - pow(nu + mdeu, 2)) / (2. * qnorm * p1norm));
        p2norm = sqrt(pow(nu + mdeu - sqrt(pow(mpro, 2) + pow(p1norm, 2)), 2) - pow(mpro, 2));

        p1Fv[0] = sqrt(pow(mpro, 2) + pow(p1norm, 2));
        p1Fv[1] = p1norm * sin(theta1) * cos(phi1);
        p1Fv[2] = p1norm * sin(theta1) * sin(phi1);
        p1Fv[3] = p1norm * cos(theta1);

        P0Fv[0] = mdeu;
        P0Fv[1] = 0.;
        P0Fv[2] = 0.;
        P0Fv[3] = 0.;

        p2Fv[0] = sqrt(pow(mpro, 2) + pow(p2norm, 2));
        p2Fv[1] = -p1norm * sin(theta1) * cos(phi1);
        p2Fv[2] = -p1norm * sin(theta1) * sin(phi1);
        p2Fv[3] = qnorm - p1norm * cos(theta1);

        qFv[0] = nu;
        qFv[1] = 0.;
        qFv[2] = 0.;
        qFv[3] = qnorm;

    }

    void setpwMomdotvar()
    {
        double Q2 = pow(qnorm, 2) - pow(nu, 2);
        pwMomdotvar[0] = mpro;
        pwMomdotvar[1] = mdeu;
        if(p_type==2112){
            pwMomdotvar[2] = formFun::neuF1(Q2);
            pwMomdotvar[3] = formFun::neuF2(Q2);
            //std::cout<<"using neutron funct"<<std::endl;
        } else{
            pwMomdotvar[2] = formFun::proF1(Q2);
            pwMomdotvar[3] = formFun::proF2(Q2);
            //std::cout<<"using proton funct"<<std::endl;
        }

        pwMomdotvar[4] = leptonTensor::Fvprod(P0Fv, p2Fv);
        pwMomdotvar[5] = leptonTensor::Fvprod(P0Fv, qFv);
        pwMomdotvar[6] = leptonTensor::Fvprod(p2Fv, qFv);
        pwMomdotvar[7] = leptonTensor::Fvprod(qFv, qFv);

        for (int i = 0; i < 4; i++)
        {
            pwmomFv[0][i] = 0.0;
            pwmomFv[1][i] = P0Fv[i];
            pwmomFv[2][i] = p2Fv[i];
            pwmomFv[3][i] = qFv[i];

        }

    }

    void importpwcoeMat()
    {
        int rowdim;
        const int coemomdimgmn = 2;
        const int coemomdimp = 1;
        const int varprenum = 4;
        //denominator's largest mass dimension  m^8   8/2=4       
        //nucleon electromagnetic form factor 2   deuteron form factor  4  deuteron polarization summation  2      
        const int demaxdim = 4;

        const int coldimgmn = varprenum + demaxdim + coemomdimgmn + 1;
        const int coldimp = varprenum + demaxdim + coemomdimp + 1;
        int index;
        for (int i = 0; i < 4; i++)
        {
            for (int j = i; j < 4; j++)
            {
                ifstream readDimmat(strpnCoe[0] + strpnCoedim[0] + strtenName[0] + strdeuForm[i] + strdeuForm[j] + strdimfileExt);
                assert(readDimmat.is_open());
                for (int ii = 0; ii < npp; ii++)
                {
                    for (int jj = 0; jj < npp; jj++)
                    {
                        readDimmat >> pwCoeDimmat[i][j][ii][jj];
                    }
                }
                readDimmat.close();

                // import exponent vector matrix for g_mu_nu tensor
                ifstream readexpVecmat(strpnCoe[0] + strtenName[0] + strdeuForm[i] + strdeuForm[j] + strfilegmnExt);
                assert(readexpVecmat.is_open());

                rowdim = pwCoeDimmat[i][j][0][0];

                index = 0;
                pwCoeVecmat[i][j][index] = allocateMatrix(rowdim, coldimgmn);

                for (int kk = 0; kk < rowdim; kk++)
                {
                    for (int nn = 0; nn < coldimgmn; nn++)
                    {
                        readexpVecmat >> pwCoeVecmat[i][j][index][kk][nn];
                    }
                }
                readexpVecmat.close();

                // import exponent vector matrix for rank two momentum tensor
                for (int ii = 1; ii < npp; ii++)
                {
                    for (int jj = ii; jj < npp; jj++)
                    {
                        ifstream readexpVecmat(strpnCoe[0] + strtenName[0] + strdeuForm[i] + strdeuForm[j]
                            + strmompw[ii - 1] + strmompw[jj - 1] + strfileExt);
                        assert(readexpVecmat.is_open());
                        rowdim = pwCoeDimmat[i][j][ii][jj];

                        //the index for a certain momentum rank 2 tensor
                        //index for g_mu_nu  is  0
                        index = indexFun(ii, jj, npp);
                        pwCoeVecmat[i][j][index] = allocateMatrix(rowdim, coldimp);
                        for (int kk = 0; kk < rowdim; kk++)
                        {
                            for (int nn = 0; nn < coldimp; nn++)
                            {
                                readexpVecmat >> pwCoeVecmat[i][j][index][kk][nn];
                            }
                        }
                        readexpVecmat.close();
                    }
                }


            }
        }

        cout << "import exponent vector matrix for plane wave approximation" << endl;
    }

    void setpwcoeVal()
    {
        int rowdim;
        //  +4  m md fo1 fo2      +4   denominator's largest dimension  m^8   8/2=4
        //  +1  p2_a q_b   coe     dimension of coefficient  for plane wave is 1             +1    coefficient value
        //  the coldim  for gmn term   should add another 1, as the dimension of gmn is 0
        const int coemomdimgmn = 2;
        const int coemomdimp = 1;
        const int varprenum = 4;
        const int demaxdim = 4;        //denominator's largest mass dimension  m^8   8/2=4

        const int coldimgmn = varprenum + demaxdim + coemomdimgmn + 1;
        const int coldimp = varprenum + demaxdim + coemomdimp + 1;
        int index;
        //int pdotpnum;

        //double coevalmono = 1.0;
        //double coeval = 0.0;

        for (int i = 0; i < 4; i++)
        {
            for (int j = i; j < 4; j++)
            {
                rowdim = pwCoeDimmat[i][j][0][0];
                index = 0;

                pwCoemat[i][j][index] = polyval(rowdim, coldimgmn, pwMomdotvar, pwCoeVecmat[i][j][index], coemomdimgmn, varprenum);

                //coeval = 0.0;
                //coevalmono = 1.0;

                for (int k = 1; k < npp; k++)
                {
                    for (int n = k; n < npp; n++)
                    {
                        rowdim = pwCoeDimmat[i][j][k][n];
                        //index = (k - 1) * (npp - 1 + (npp - 1 - k + 2)) / 2 + n + 1 - k;
                        index = indexFun(k, n, npp);
                        if (rowdim == 0)
                        {
                            pwCoemat[i][j][index] = 0.0;
                        }
                        else
                        {
                            pwCoemat[i][j][index] = polyval(rowdim, coldimp, pwMomdotvar, pwCoeVecmat[i][j][index], coemomdimp, varprenum);

                        }

                    }
                }
            }
        }
    }

    double pwNume()
    {
        double deuform[4] = { myform->deufd(p2norm), myform->deugd(p2norm), myform->deuhd(p2norm), myform->deuid(p2norm) };
        int index;
        double pwnume = 0.;
        double deucoe;
        double contrmat[momtenpwNum + 1];

        contractFun(leptonTensor::lepcoe, leptonTensor::lepmomFv, pwmomFv, npp, contrmat);

        for (int i = 0; i < 4; i++)
        {
            for (int j = i; j < 4; j++)
            {
                deucoe = deuform[i] * deuform[j];

                index = 0;
                pwnume = pwnume + deucoe * pwCoemat[i][j][index] * contrmat[index];

                for (int k = 1; k < npp; k++)
                {
                    for (int n = k; n < npp; n++)
                    {
                        index = indexFun(k, n, npp);
                        pwnume = pwnume + deucoe * pwCoemat[i][j][index] * contrmat[index];
                    }
                }
            }
        }
        return pwnume;
    }

    //plane wave contribution to total differential cross section
    double pwdcsfun(double Egama, double pgama, double thetaea, double phiea, double slla, double tqsqa, double spna, double p1a, double phi1a)
    {
        double nume, deno;
        double pwampsq, pwdcs;

        double pcnorm;

        leptonTensor::setlepmomVar(Egama, pgama, thetaea, phiea, slla, spna, tqsqa);
        setnuclmomVar(p1a, phi1a, spna, tqsqa);

        //set the momentum dot variables  and the momentum tensor coefficient
        setpwMomdotvar();
        setpwcoeVal();

        nume = pwNume();
        deno = pow(mdeu, 2) - 2. * mdeu * sqrt(pow(mpro, 2) + pow(p2norm, 2));

        //std::cout<<"mpro "<<mpro<<std::endl;

        pwampsq = nume / (pow(2. * mpro, 2) * pow(deno, 2));

        pcnorm = sqrt(slla / 4. - pow(leptonTensor::mlepton, 2));

        // differential variable d p_1  d_phi_1  d s_ll  d s_pn  d t  d theta_lep  d phi_lep  d phi_p34
        // azimuthal angle with respect to p_gam remains
        // sin(thetaea) for lepton pair phase space
        pwdcs = 1. / 3. * pow(alp, 3) / (2. * pow(4. * M_PI, 5)) * p1norm / (pow(mdeu, 2) * sqrt(pow(mpro, 2) + pow(p1norm, 2)) * qnorm) * pwampsq
            * pcnorm / (pow(Egama, 2) * sqrt(slla) * pow(tqsqa, 2)) * pow(2. * mpro, 2) * sin(thetaea);

        return pwdcs;
    }

}


namespace incidentPhoton
{
    // code for bremsstrahlung photon and equivalent virtual photon 
    // from Tianbo Liu   "Lcore.h"
    TRandom3 random(0);

    TF1* TF_fBremsstrahlung;

    double cthrange[2] = { -1.0, 1.0 };
    double perange[2] = { 0.0, 10.0 };

    double Bremsstrahlung(const double* y, const double* par) 
    {//ds/dy approximate expression
      //E0: electron beam energy; k: photon energy
        if (y[0] < 0.01) 
        {// Infrared cut
            std::cerr << "Out of range in Bremsstrahlung!" << std::endl;
            return -1.0;
        }
        double result = (4.0 / 3.0 - 4.0 / 3.0 * y[0] + y[0] * y[0]) / y[0];
        return result;
    }

    double BremsstrahlungPhoton(TLorentzVector* q, const double kmin, const double kmax, const double E,const double l, const double X0) 
    { //Generate a Bremsstrahlung photon ! 
      //q: photon; E: electron beam energy; [kmin, kmax]: photon energy range
        double ymin = kmin / E;
        double ymax = kmax / E;
        double y = TF_fBremsstrahlung->GetRandom(ymin, ymax);
        q->SetXYZT(0.0, 0.0, y * E, y * E);

        // random number is normalized, so multiply a normalization factor
        return 0.5*(l/X0) * (4.0 / 3.0 * log(ymax / ymin) - 4.0 / 3.0 * (ymax - ymin) + 1.0 / 2.0 * (ymax * ymax - ymin * ymin));
    }

    int SetBremsstrahlung() 
    {
        TF_fBremsstrahlung = new TF1("fBremsstrahlung", Bremsstrahlung, 0.01, 1.0, 0);
        TF_fBremsstrahlung->SetNpx(1000);
        return 0;
    }

    double VirtualPhoton(const TLorentzVector* ki, TLorentzVector kf[2]) 
    {   //ki: e   kf: e', gamma
        const double couple = 4.0 * M_PI * leptonTensor::alp;
        const double m = leptonTensor::melec;

        double Egamin, Eelecin;
        
        double Pe = random.Uniform(perange[0], perange[1]);
        double cth = random.Uniform(cthrange[0], cthrange[1]);
        double sth = sqrt(1.0 - cth * cth);
        double phi = random.Uniform(-M_PI, M_PI);
        kf[0].SetXYZM(Pe * sth * cos(phi), Pe * sth * sin(phi), Pe * cth, m);//e'
        kf[1] = *ki - kf[0];//virtual photon

        // Egamin energy of virtual photon
        Egamin = kf[1].E();
        Eelecin = ki->E();

        double Q2 = -kf[1] * kf[1];//Q^2 = -q^2
        double fluxRatio = Egamin / Eelecin;
        double amp = (2.0 * Q2 - 4.0 * m * m) / (Q2 * Q2);
        double phase = kf[0].P() * kf[0].P() / (2.0 * kf[0].E() * pow(2.0 * M_PI, 3));
        double volume = 2.0 * M_PI * (perange[1] - perange[0]) * (cthrange[1] - cthrange[0]);
        double y = fluxRatio;
        double gamy = sqrt(Q2) / Egamin;
        double epsilon = (1.0 - y - 0.25 * gamy * gamy) / (1.0 - y + 0.5 * y * y + 0.25 * gamy * gamy);

        return couple * fluxRatio * amp * phase * volume / (1.0 - epsilon);
    }
}

namespace BH_deuteron
{
    double Md = nuclearTensor::mdeu;
    double Mp = 0.938272;
    double maxEnergy = 8.5;

    const double anglow = 5.0 * M_PI / 180.0;
    const double angup = 40.0 * M_PI / 180.0;

    // 2212 is proton, 2112 is neutron
    int target_type=2212;

    TRandom3 myrandom(0);

    //int eventNum;

    void setModel(int lepType, double maxene, int in_target_type)
    {
        //eventNum = eventNum_in;
        target_type=in_target_type;

        maxEnergy = maxene;               // max incident energy

        if (lepType == 0) {
            leptonTensor::mlepton = leptonTensor::melec;
        }
        else if (lepType == 1) {
            leptonTensor::mlepton = leptonTensor::mmuon;
        }
        else {
            cout << "input wrong number for the final state lepton type !" << endl;
        }

        nuclearTensor::setPType(target_type);
        nuclearTensor::importpwcoeMat();
        nuclearTensor::myform->forminterpInit();

        if(target_type==2112){
            Mp = 0.939565;
        } 
        
    }

    // set the approximate maximum range of differential variables for a certain incident photon energy pgam
    // when combining with virtual photon or Bremsstrahlung photon scattering, input largest momentum of gamma photon
    void setVarRan(double ranlow[8], double ranup[8])
    {
        double pgam = maxEnergy;
        double sllmin = 4. * pow(leptonTensor::mlepton, 2);
        double sllmax = pow(3.35, 2);
        double spnmin = 4. * pow(Mp, 2);
        // deuteron target rest frame, invariant mass squred for real photon
        double sGamD = pow(Md, 2) + 2. * Md * pgam;

        //std::cout<<"Mp "<<Mp<<std::endl;

        // sll 
        ranlow[0] = sllmin;
        ranup[0] = sllmax;
        //ranup[0] = pow(sqrt(sGamD) - sqrt(spnmin), 2);
        
        // spn
        ranlow[1] = spnmin;
        ranup[1] = pow(sqrt(sGamD) - sqrt(sllmin), 2);
        // t
        ranlow[2] = (-(pow(Md, 3) * pgam) + pow(Md, 2) * (-2. * pow(pgam, 2) + sllmin) + Md * pgam * (sllmin + spnmin) 
            - sqrt(pow(Md, 2) * pow(pgam, 2) * (pow(pow(Md, 2) + 2. * Md * pgam - sllmin, 2) - 2. * (pow(Md, 2) + 2. * Md * pgam + sllmin) * spnmin 
                + pow(spnmin, 2)))) / (Md * (Md + 2. * pgam));
        ranup[2] = 0.0;
        // p1norm
        ranlow[3] = 0.0;
        ranup[3] = 5.0;
        //ranup[3] = sqrt(pow(pgam + Md - 2. * leptonTensor::mlepton - Mp, 2) - pow(Mp, 2));
        
        // phi_p
        ranlow[4] = -M_PI;
        ranup[4] = M_PI;
        // theta_lep
        ranlow[5] = 0.0;
        ranup[5] = M_PI;
        // phi_lep
        ranlow[6] = -M_PI;
        ranup[6] = M_PI;
        // phi_gam   azimuthal angle with respect to p_gam, remains for the generator 
        ranlow[7] = -M_PI;
        ranup[7] = M_PI; 
    }

    bool physRegion(double Egam_in, double pgam_in, double sll_in, double t_in, double spn_in, double p1_in)
    {
        // t is negative
        //double spnmax, tmin, tmax, p1min, p1max;
        double nu_in, qnorm_in, cos_thetap1, p34norm, cos_thetap34;

        /*tmin = (4 * pow(Mp, 2) * Egam_in - Md * Egam_in * (Md + 2 * Egam_in) + (Md + Egam_in) * sll_in
            - sqrt(pow(Egam_in, 2) * (16 * pow(Mp, 4) + pow(pow(Md, 2) + 2 * Md * Egam_in - sll_in, 2)
                - 8 * pow(Mp, 2) * (pow(Md, 2) + 2 * Md * Egam_in + sll_in)))) / (Md + 2 * Egam_in);
        tmax = (4 * pow(Mp, 2) * Egam_in - Md * Egam_in * (Md + 2 * Egam_in) + (Md + Egam_in) * sll_in
            + sqrt(pow(Egam_in, 2) * (16 * pow(Mp, 4) + pow(pow(Md, 2) + 2 * Md * Egam_in - sll_in, 2)
                - 8 * pow(Mp, 2) * (pow(Md, 2) + 2 * Md * Egam_in + sll_in)))) / (Md + 2 * Egam_in);
        spnmax = (2. * Md * Egam_in - sll_in + t_in) * (Md * sll_in - t_in * (Md + 2. * Egam_in)) / (2. * Egam_in * (sll_in - t_in));*/

        if ((pow(Egam_in, 2) - pow(pgam_in, 2) + pow(Md, 2) + 2. * Md * Egam_in) > (pow(2. * leptonTensor::mlepton + 2. * Mp, 2)))
        {
            nu_in = (spn_in - t_in - pow(Md, 2)) / (2. * Md);
            qnorm_in = sqrt(pow(nu_in, 2) - t_in);

            if (nu_in < Egam_in && (pow(Egam_in - nu_in, 2) >= sll_in))
            {
                p34norm = sqrt(pow(Egam_in - nu_in, 2) - sll_in);
                cos_thetap34 = (pow(pgam_in, 2) + pow(p34norm, 2) - pow(qnorm_in, 2)) / (2. * pgam_in * p34norm);

                if (cos_thetap34 <= 1. && cos_thetap34 >= -1.)
                {
                    // nu + Md = sqrt( Mp^2 + p1norm^2 ) + sqrt( Mp^2 + ( vec{q} - vec{p1} )^2 )
                    cos_thetap1 = (pow(qnorm_in, 2) + 2. * (nu_in + Md) * sqrt(pow(Mp, 2) + pow(p1_in, 2)) - pow(nu_in + Md, 2)) / (2. * qnorm_in * p1_in);

                    if (cos_thetap1 <= 1. && cos_thetap1 >= -1.)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }

    double GetBHdeu(TLorentzVector* q_in, TLorentzVector p_out[3], double var[10])
    {
        // q_in: incident photon momentum
        // p_out[0]: proton momentum, p_out[1]: lepton momentum, p_out[2]: anti-lepton momentum
        double Egam_in = q_in->E();
        double pgam_in = q_in->P();
        double sll_in, spn_in, t_in, p1_in, phi1_in, thetalep_in, philep_in, phigam_in;
        double phaseVol = 1.;
        double ranlow[8], ranup[8];
        double dsigma = 0.;
        
        // p34 = p3 + p4       p3: lepton momentum, p4: anti-lepton momentum
        double nu_in, qnorm_in, theta_p34, theta_q, p34norm, theta_p1, pclepnorm;
        double p34Vz;
        
        setVarRan(ranlow, ranup);

        for (int i = 0; i < 8; i++)
        {
            phaseVol = phaseVol * (ranup[i] - ranlow[i]);
        }

        sll_in =      myrandom.Uniform(ranlow[0], ranup[0]);
        spn_in =      myrandom.Uniform(ranlow[1], ranup[1]);
        t_in =        myrandom.Uniform(ranlow[2], ranup[2]);
        p1_in =       myrandom.Uniform(ranlow[3], ranup[3]);                    // proton momentum
        phi1_in =     myrandom.Uniform(ranlow[4], ranup[4]);                    // azimuthal angle of proton along q direction 
        thetalep_in = myrandom.Uniform(ranlow[5], ranup[5]);                    // polar angle of lepton in lepton pair cm frame
        philep_in =   myrandom.Uniform(ranlow[6], ranup[6]);                    // azimuthal angle of lepton in lepton pair cm frame
        phigam_in =   myrandom.Uniform(ranlow[7], ranup[7]);                    // azimuthal angle of incident photon, azimuthal symmetry

        if (physRegion(Egam_in, pgam_in, sll_in, t_in, spn_in, p1_in))
        {
            // calculate momentum 
            nu_in = (spn_in - t_in - pow(Md, 2)) / (2. * Md);
            qnorm_in = sqrt(pow(nu_in, 2) - t_in);

            pclepnorm = sqrt(sll_in / 4. - pow(leptonTensor::mlepton, 2));
            theta_p1 = acos((pow(qnorm_in, 2) + 2. * (nu_in + Md) * sqrt(pow(Mp, 2) + pow(p1_in, 2)) - pow(nu_in + Md, 2)) / (2. * qnorm_in * p1_in));
            p34norm = sqrt(pow(Egam_in - nu_in, 2) - sll_in);
            // polar angle between p34 and p_gam
            theta_p34 = acos((pow(pgam_in, 2) + pow(p34norm, 2) - pow(qnorm_in, 2)) / (2. * pgam_in * p34norm));
            // polar angle between q and p_gam
            theta_q = acos((pow(pgam_in, 2) + pow(qnorm_in, 2) - pow(p34norm, 2)) / (2. * pgam_in * qnorm_in));

            p34Vz = p34norm / (Egam_in - nu_in);

            //final momentum in the frame where  pgam_x = p_gam_y = 0, phi_gam = 0
            p_out[0].SetXYZM(p1_in * sin(theta_p1) * cos(phi1_in), p1_in * sin(theta_p1) * sin(phi1_in), p1_in * cos(theta_p1), Mp);
            p_out[0].RotateY(-theta_q);

            p_out[1].SetXYZM(pclepnorm * sin(thetalep_in) * cos(philep_in), pclepnorm * sin(thetalep_in) * sin(philep_in),
                pclepnorm * cos(thetalep_in), leptonTensor::mlepton);
            p_out[2].SetXYZM(-pclepnorm * sin(thetalep_in) * cos(philep_in), -pclepnorm * sin(thetalep_in) * sin(philep_in),
                -pclepnorm * cos(thetalep_in), leptonTensor::mlepton);

            p_out[1].Boost(0., 0., p34Vz);
            p_out[2].Boost(0., 0., p34Vz);

            p_out[1].RotateY(theta_p34);
            p_out[2].RotateY(theta_p34);

            for (int i = 0; i < 3; i++)
            {
                // rotate to a random phi angle in a frame where the quasi-real photon momentum is along z axis
                p_out[i].RotateZ(phigam_in);

                // rotate to the original electron laboratory frame
                p_out[i].RotateY(q_in->Theta());
                p_out[i].RotateZ(q_in->Phi());
            }

            // limit the polar angle of e^+ and e^-
            if ((p_out[1].Theta() < anglow) || (p_out[1].Theta() > angup) || (p_out[2].Theta() < anglow) || (p_out[2].Theta() > angup) || (p_out[0].Theta() < anglow) || (p_out[0].Theta() > angup))
            {
                dsigma = -1.0;
            }
            else
            {
                
                dsigma = nuclearTensor::pwdcsfun(Egam_in, pgam_in, thetalep_in, philep_in, sll_in, t_in, spn_in, p1_in, phi1_in);
                
                dsigma = dsigma * phaseVol;
                //dsigma = dsigma * phaseVol / eventNum;

                var[0] = Egam_in;
                var[1] = thetalep_in;
                var[2] = philep_in;
                var[3] = sll_in;
                var[4] = t_in;
                var[5] = spn_in;
                var[6] = p1_in;
                var[7] = phi1_in;
                var[8] = phigam_in;
                var[9] = pgam_in;
            }
        }
        else
        {
            dsigma = -1.0;

            for (int i = 0; i < 3; i++)
            {
                p_out[i].SetPxPyPzE(0., 0., 0., 0.);
            }

            for (int i = 0; i < 10; i++)
            {
                var[i] = 0.;
            }
        }

        return dsigma;
    }
}

