//
// Created by afrotscher on 1/7/19.
//

#include <TF1.h>
#include <TGraph.h>
#include <TMinuit.h>
#include <Math/Vector3D.h>
#include "minos.h"
#include "TCanvas.h"
#include "TMinuitMinimizer.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

using std::vector, std::cerr, std::cout, std::endl, std::min, std::max;

TMinosPass minosana::analyze() {
    /// MINOS 2. Modify with hough-transformation
    vector<bool> clusterringbool;
    vector<int> clusternbr, clusterpads;
    
    TMinosPass MR;

    if(!filled) return MR;  // return if no TPC has been hit

    for (int iteration = 0; (Xpad.size() > 9 && iteration < 20); iteration++) {
        const int filter_result = Obertelli_filter(
                Xpad, Ypad, Qpad, Xpadnew,Ypadnew, Qpadnew,
                clusterringbool);

        const vector<int> temp1(filter_result, iteration);
        const vector<int> temp2(filter_result, filter_result);
        clusternbr.insert(clusternbr.end(), temp1.begin(), temp1.end());
        clusterpads.insert(clusterpads.end(), temp2.begin(), temp2.end());

        if (filter_result > 10 && clusterringbool.back()) MR.trackNbr++;
    }
    
    /// Bonus for Tracknumber between 1 and 4
    if(MR.trackNbr < 1 || MR.trackNbr > 4) return MR;
    
    for (unsigned long i = 0; i < Xpadnew.size(); i++) {
        fitdata.add(Xpadnew[i], Ypadnew[i], -1E4, -1E4,
                    Qpadnew[i], clusternbr[i], clusterpads[i], 0);
        Zpadnew.push_back(-1E4);
    }

    /// MINOS 3: Fitting the taken pads for Qmax and Ttrig information /
    //padsleft -= Xpadnew.size();
    for (unsigned long i = 0; i < minoscalibvalues.size(); i++) {

        const double x_mm = minostrackxy.at(i).at(0);
        const double y_mm = minostrackxy.at(i).at(1);
        bool fitbool = false;
        int indexfill = 0;

        for (unsigned long j = 0; j < Xpadnew.size(); j++) {
            if (abs(Xpadnew[j] - x_mm) < 0.01 && abs(Ypadnew[j] - y_mm) < 0.01) {
                fitbool = true;
                indexfill = (int)j;
                break;
            }
        }
        // Check if new channel is of interest
        if (!fitbool) continue;

        // if so, we read Q(t), and fill vectors w/ t&Q info after fitting E(t)
        TH1F hfit(Form("hfit%i",threadno), Form("hfit%i",threadno),
                  512, 0, 512);
        for (unsigned long k = 0; k < minoscalibvalues.at(i).size(); k++) {
            if (minoscalibvalues.at(i).at(k) >= 0)
                hfit.SetBinContent(hfit.FindBin(minostime.at(i).at(k)),
                                            minoscalibvalues.at(i).at(k) + 250);
        }
        // Fitting the hfit histogram of last ch. if not empty
        if (hfit.GetSumOfWeights() == 0) continue;

        hfit.GetXaxis()->SetRange(0, 510);
        const double hfit_max = hfit.GetMaximum();
        const double hfit_max_T = hfit.GetMaximumBin();

        // Find T_min and T_max limits for non-0 signals
        double T_min = hfit.FindLastBinAbove(250);
        double T_max = hfit.FindFirstBinAbove(0) - 1;

        // Take only 1.5*shaping time before max if other signals before
        if (hfit_max_T - 3.5 * (Tshaping / TimeBinElec) > T_min)
            T_min = hfit_max_T - 2 * Tshaping / TimeBinElec;
        if ((hfit_max_T + 10) < T_max || T_max == -1) T_max = hfit_max_T + 10;

        T_min = std::max(T_min, 0.);
        T_max = std::min(T_max, 510.);

        // Set fit parameters
        TF1 fit_function(Form("fit_function%i",threadno),
                [](double *x, double *p){ /// Check for boundaries of x
            if(x[0] < p[1] || x[0] > 512) return 250.;
            else return p[0] * exp(-3.*(x[0]-p[1])/p[2])*sin((x[0]-p[1])/p[2]) *
                        pow((x[0]-p[1])/p[2], 3) + 250;}, 0, 511, 3);

        fit_function.SetParameters( hfit_max - 250,
                   hfit_max_T - Tshaping / TimeBinElec, Tshaping / TimeBinElec);
        fit_function.SetParNames("Amplitude", "trigger time", "shaping time");
        fit_function.SetParLimits(0, 0, 1E5);
        fit_function.SetParLimits(1, -20, 512);
        fit_function.SetParLimits(2, 0, 512);

        // --> parameter n (no store no draw) crucial for multithread <--
        const int fit2DStatus = hfit.Fit(&fit_function, "QN", "", T_min, T_max);

        double fit_function_max = 0, fit_function_Tpad = 0, Chi2 = 0;

        if (!fit2DStatus) {
            Chi2 = fit_function.GetChisquare();
            fit_function_max = fit_function.GetMaximum();
            fit_function_Tpad = fit_function.GetParameter(1);
        }

        // attribute q_pad and z_mm value
        double q_pad = 0, z_mm = 0, t_pad = 0;
        if (fit2DStatus || fit_function_max <= 20 ||
            fit_function_max > 1E5 || fit_function_Tpad < .15 ||
            fit_function_Tpad > 512 || fit_function.GetParameter(2) < .15 ||
            fit_function.GetParameter(2) > 512) {
            // No correct fit here :(
            q_pad = hfit_max - 250;
            z_mm = -1E4;
        }
        else { // Add to the variables the fit parameter
            t_pad = fit_function_Tpad;
            z_mm = (t_pad * TimeBinElec - DelayTrigger) * VDrift;
            q_pad = fit_function_max - 250;
        }

        fitdata.replace(Xpadnew[indexfill], Ypadnew[indexfill],
                        t_pad * TimeBinElec, z_mm, q_pad, clusternbr[indexfill],
                        clusterpads[indexfill], Chi2, indexfill);
        Zpadnew[indexfill] = z_mm;
        Qpadnew[indexfill] = q_pad;
    } // End of all entries

    /// 4: MINOS Filtering the tracks off possible noise with Hough3D //

    unsigned long padsleft2 = Xpadnew.size();
    vector<double> xin, yin, zin, qin, xout, yout, zout, qout;
    int cluster_temp = 0, array_final = 0;

    vector<TGraph> grxz(MR.trackNbr + 1), gryz(MR.trackNbr + 1);

    for (unsigned long i = 0; i < padsleft2; i++) {
        // if-Loop executed only at the End
        if(cluster_temp == clusternbr[i] && i==(Xpadnew.size()-1) &&
           clusterpads[i]>=10 && clusterringbool[i] == true &&
           Zpadnew[i] >-10000 && Zpadnew[i] <=320){
            xin.push_back(Xpadnew[i]);
            yin.push_back(Ypadnew[i]);
            zin.push_back(Zpadnew[i]);
            qin.push_back(Qpadnew[i]);
        }

        // x vector not empty AND (new track OR last entry of the loop)
        if (!xin.empty() && ((cluster_temp +1 == clusternbr[i] && i) ||
                           i == (Xpadnew.size() - 1))) {
            int ringsum = 0;
            double zmax = 0;
            vector<int> ringtouch(18, 0);
            Hough_filter(xin, yin, zin, qin, xout, yout, zout, qout);

            for (unsigned long k = 0; k < xout.size(); k++) {
                zmax = max(zmax, zout[k]);
                ringtouch.at(max(0, int((sqrt(pow(xout[k], 2.) +
                                           pow(yout[k], 2.)) - 45.2) / 2.1)))++;
            }
            for (auto &j: ringtouch) if (j > 0) ringsum++;
            if (zmax > 290) ringsum = 16; // particles escaping through the back

            // Decide whether we have a track, and add to 2D-plane TGraphs
            if (xout.size() > 10 && ringsum >= 15) {
                MR.trackNbr_final++;
                //cluster1 = cluster_temp; delete if no use
                for (unsigned long l = 0; l < xout.size(); l++) {
                    dataresult.add(xout[l], yout[l], zout[l],
                                   qout[l], MR.trackNbr_final, xout.size(), zmax);
                    array_final++;
                    if(MR.trackNbr_final > (int)grxz.size()) continue;
                    grxz.at(MR.trackNbr_final - 1).
                     SetPoint(grxz.at(MR.trackNbr_final - 1).GetN(), zout[l], xout[l]);
                    gryz.at(MR.trackNbr_final - 1).
                     SetPoint(gryz.at(MR.trackNbr_final - 1).GetN(), zout[l], yout[l]);
                }
            }

            for (auto a:{&xin, &yin, &zin, &qin, &xout, &yout, &zout, &qout})
                a->clear();
        }

        cluster_temp = clusternbr[i];
        if (clusterpads[i] >= 10 && clusterringbool[i] == true &&
              Zpadnew[i] > -1E4 && Zpadnew[i] <= 320){  // 320 => 520 debug
            xin.push_back(Xpadnew[i]);
            yin.push_back(Ypadnew[i]);
            zin.push_back(Zpadnew[i]);
            qin.push_back(Qpadnew[i]);
        }
        //else continue;
    } // end of loop on pads

    if(MR.trackNbr_final == 0 || MR.trackNbr_final > 3) return MR;

    /// 5. Fitting filtered tracks in 3D (weight by charge, TMinuit) //
    vector<double> pStart_1{0,1,0,1}, pStart_2{0,1,0,1},pStart_3{0,1,0,1},
                   parFit_1(4), err_1(4), parFit_2(4), parFit_3(4), err_2(4),
                   err_3(4), chi1(2), chi2(2), chi3(2), arglist(10),covxy(3);
    vector<int> fitStatus(2);

    int iflag, nvpar, nparx;
    double amin, edm, errdef; //, chi2res1, chi2res2;
    arglist.at(0) = 3;
    // Get fit for xz-yz plane, for the reconstruction
    FindStart(pStart_1, chi1, fitStatus, grxz.at(0), gryz.at(0));

    auto funcmin = [this](const int mod){
        /// Wrapper for minimization function, sets track number to mod
        mode = mod;
        tmr = minosana::getTMinosResult();
        return [](int &, double *, double &sum, double *par, int){
            /// Minimization: calculates sum of distances between data points
            /// and regression

            sum = 0;
            double qtot =0;
            for(unsigned long i=0; i<tmr.x_mm.size(); i++){
                if(tmr.n_Cluster.at(i) == mode){
                    double d = distancelinepoint(tmr.x_mm.at(i), tmr.y_mm.at(i),
                                                 tmr.z_mm.at(i), par);
                    sum  += d*tmr.Chargemax.at(i);
                    qtot +=   tmr.Chargemax.at(i);
                }
            }
            sum /=qtot;
        };
    };
    
    TMinuit min(4);
    min.SetPrintLevel(-1);
    min.SetFCN(funcmin(1));
    // Set starting values and step sizes for parameters
    auto minimize = [&min, &iflag, &arglist, &amin, &edm, &errdef, &nvpar, &nparx]
            (vector<double> a){
        /// Set all parameters for 3D-minimization
        assert(a.size() == 4);

        min.mnparm(0, "x0", a.at(0), 0.1, -500, 500, iflag);
        min.mnparm(1, "Ax", a.at(1), 0.1, -10, 10, iflag);
        min.mnparm(2, "y0", a.at(2), 0.1, -500, 500, iflag);
        min.mnparm(3, "Ay", a.at(3), 0.1, -10, 10, iflag);

        arglist.at(0) = 100; // Number of function calls
        arglist.at(1) = 1E-6;// tolerance

        min.mnexcm("MIGRAD", arglist.data(), 2, iflag);

        // get current status of minimization
        min.mnstat(amin, edm, errdef, nvpar, nparx, iflag);
    };
    minimize(pStart_1);

    for (unsigned long i = 0; i < parFit_1.size(); i++)
        min.GetParameter(i, parFit_1.at(i), err_1.at(i));
    double temp[4][4];
    min.mnemat(&temp[0][0], 4);
    covxy.at(0) = temp[1][3]; // entry i=0, j=2

    if (MR.trackNbr_final == 1) parFit_2 = {0, 0, 0, 0};
    else {
        FindStart(pStart_2, chi2, fitStatus, grxz.at(1), gryz.at(1));
        min.SetFCN(funcmin(2));
        minimize(pStart_2);
        min.mnemat(&temp[0][0], 4);
        covxy.at(1) = temp[1][3]; // entry i=0, j=2

        for (unsigned long i = 0; i < parFit_2.size(); i++)
            min.GetParameter(i, parFit_2.at(i), err_2.at(i));
    }
    if(MR.trackNbr_final > 2){
        FindStart(pStart_3, chi3, fitStatus, grxz.at(2), gryz.at(2));
        min.SetFCN(funcmin(3));
        minimize(pStart_3);
        min.mnemat(&temp[0][0], 4);
        covxy.at(2) = temp[1][3]; // entry i=0, j=2

        for (unsigned long i = 0; i < parFit_3.size(); i++)
            min.GetParameter(i, parFit_3.at(i), err_3.at(i));
    }
    else parFit_3 = {0,0,0,0};

    //Rotate from MINOS to beamline for all tracks
    double rot = 30*TMath::Pi()/180;
    vector<double> parFit_1r = rotatesp(rot, parFit_1);
    vector<double> parFit_2r = rotatesp(rot, parFit_2);
    vector<double> parFit_3r = rotatesp(rot, parFit_3);

    /// 6. Get x,y,z reaction vertex from fitted parameters and further variables
    vertex(parFit_1r, parFit_2r, MR.x_vertex, MR.y_vertex, MR.z_vertex);
    
    MR.r_vertex = sqrt(pow(MR.x_vertex, 2) + pow(MR.y_vertex, 2));

    // cos theta = v1*v2/(|v1|*|v2|), v1 = (0,0,1)^T, v2 = (m1,m3,1)^T
    auto cthet = [](vector<double> a){
        assert(a.size() == 4);
        return acos(1/pow(1 + pow(a[1],2) + pow(a[3],2),.5)) *
               180./TMath::Pi();
    };

    auto thetaerror = [](vector<double> pf, vector<double> err, double covxy){
        /// Calculates the Gaussian error of the theta angle for each track
        assert(pf.size() == 4);
        double norm = pf[1]*pf[1]+pf[3]*pf[3]+1;
        return 180./TMath::Pi()*sqrt(pow(pf[1]*err[1]/(pow(norm,3/2.) *
                                                   sqrt(1-1/(norm))),2.) +
                                 pow(pf[3]*err[3]/(pow(norm,3/2.) *
                                                   sqrt(1-1/(norm))),2.) +
                                 2*pf[1]*pf[3]*covxy/(pow(norm,3)*(1-1/norm)));
    };
    
    MR.theta = {cthet(parFit_1r), cthet(parFit_2r), cthet(parFit_3r)};
    
    MR.thetaerr = {thetaerror(parFit_1r, err_1, covxy.at(0)),
                   thetaerror(parFit_2r, err_2, covxy.at(1)),
                   thetaerror(parFit_3r, err_3, covxy.at(2))};

    // Only calculate angles for real tracks
    while(MR.trackNbr_final < (int)MR.theta.size()){
        MR.theta.pop_back();
        MR.thetaerr.pop_back();
    }
    sort(MR.theta.begin(), MR.theta.end()); // Sort angles

  // angles in xy-plane for all tracks
    vector<double> phi2d_unsorted{
        atan2(parFit_1r.at(3), parFit_1r.at(1))*180/TMath::Pi()};
    
    if(MR.trackNbr_final > 1){
        phi2d_unsorted.push_back(atan2(parFit_2r.at(3),
                                       parFit_2r.at(1)) * 180 / TMath::Pi());
        if(MR.trackNbr_final > 2)  // Add one angle for each additional track
            phi2d_unsorted.push_back(atan2(parFit_3r.at(3),
                                           parFit_3r.at(1)) * 180 / TMath::Pi());
    }

    // Transform from negative angles to positive angles
    for(auto &i:phi2d_unsorted) if(i < 0) i += 360;
    sort(phi2d_unsorted.begin(), phi2d_unsorted.end());
    
    if(MR.trackNbr_final == 2){
        MR.phi2d = {phi2d_unsorted.at(1) - phi2d_unsorted.at(0),
                    phi2d_unsorted.at(0) + 360 - phi2d_unsorted.at(1)};
    }
    else if(MR.trackNbr_final == 3){
        MR.phi2d = {phi2d_unsorted.at(2) - phi2d_unsorted.at(1),
                    phi2d_unsorted.at(1) - phi2d_unsorted.at(0),
                    phi2d_unsorted.at(0) + 360 - phi2d_unsorted.at(2)};
    }

    // Delete largest angle
    if(!MR.phi2d.empty())
        MR.phi2d.erase(std::max_element(MR.phi2d.begin(), MR.phi2d.end()));
    sort(MR.phi2d.begin(), MR.phi2d.end());
    
    // Look at some events
    if(MR.trackNbr_final == 2 && MR.phi2d.at(0) < 150 && MR.phi2d.at(0) > 140 )
        debug();

    auto lambdainter = [](auto p1, auto p2){
        /// This method calculates the relative angles between two protons
        assert(p1.size() ==4 && p2.size() == 4);
        return acos( (p1[1] * p2[1] + p1[3] * p2[3] + 1)/
                    (sqrt(pow(p1[1], 2) + pow(p1[3], 2) + 1) *
                     sqrt(pow(p2[1], 2) + pow(p2[3], 2) + 1))) *
               180. / TMath::Pi();
    };

    // Calculate all interangles between the protons
    if(MR.trackNbr_final > 1){
        MR.lambda.push_back(lambdainter(parFit_1r, parFit_2r));
        if(MR.trackNbr_final > 2){
            MR.lambda.push_back(lambdainter(parFit_1r, parFit_3r));
            MR.lambda.push_back(lambdainter(parFit_3r, parFit_2r));
        }
    }

    double radin = 40, radout = 95;
    auto param = [](vector<double> pf, double radin){
        /// This function calculates the t for which the track crosses the TPC border
        assert(pf.size() == 4);
        return -(pf[1]*pf[0]+pf[3]*pf[2])/(pow(pf[1],2)+pow(pf[3],2)) +
               sqrt(pow((pf[1]*pf[0]+pf[3]*pf[2])/(pow(pf[1],2)+pow(pf[3],2)),2)-
          (pow(pf[0],2)+pow(pf[2],2)-pow(radin,2))/(pow(pf[1],2)+pow(pf[3],2)));
    };

    for(int i=1; i <= MR.trackNbr_final; i++){
        MR.chargeweight.push_back(0);
        for(unsigned long j=0; j<tmr.n_Cluster.size(); j++){
            if(i == tmr.n_Cluster.at(j)) MR.chargeweight.back() +=
                                         tmr.Chargemax.at(j);
        }
        /// Get length of track
        double dt = 0;

        if     (i==1){dt = abs(param(parFit_1r,radin)-param(parFit_1r,radout));
                      dt *= sqrt(1 + pow(parFit_1r[1],2)+ pow(parFit_1r[3],2));}
        else if(i==2){dt = abs(param(parFit_2r,radin)-param(parFit_2r,radout));
                      dt *= sqrt(1 + pow(parFit_2r[1],2)+ pow(parFit_2r[3],2));}
        else if(i==3){dt = abs(param(parFit_3r,radin)-param(parFit_3r,radout));
                      dt *= sqrt(1 + pow(parFit_3r[1],2)+ pow(parFit_3r[3],2));}

        if(dt > 0) MR.chargeweight.back() /= dt;
    }
    
    if(MR.trackNbr_final == 2)
        MR.vertexdist = {distancelineline(parFit_1r, parFit_2r)};
    else if(MR.trackNbr_final == 3){
        MR.vertexdist = {distancelineline(parFit_1r, parFit_2r),
                         distancelineline(parFit_1r, parFit_3r),
                         distancelineline(parFit_2r, parFit_3r)};
    }

    // Estimate the angular uncertainty in lambda
    auto phierror = [](vector<double> pf, vector<double> err, double covxy){
        /// Calculates the Gaussian error of the lambda angle for each track
        /// projection
        assert(pf.size() == 4);
        return 180./TMath::Pi()/(pow(pf[1],2)+ pow(pf[3],2))*sqrt(
                pow(pf[1]*err[3],2) + pow(-pf[3]*err[1],2) - 2*pf[1]*pf[3]*covxy);
    };
    
    MR.phi2dE = {
            phierror(parFit_1r, err_1, covxy.at(0)),
            phierror(parFit_2r, err_2, covxy.at(1)),
            phierror(parFit_3r, err_3, covxy.at(2))
    };
    // Delete all errors not related to a real track
    while(MR.trackNbr_final < (int)MR.phi2dE.size()) MR.phi2dE.pop_back();
    
    return MR;
}

int minosana::Obertelli_filter(vector<double> &x, vector<double> &y,
                               vector<double> &q, vector<double> &x_out,
                               vector<double> &y_out, vector<double> &q_out,
                               vector<bool> &ringbool){
    assert(x.size() == y.size() && x.size() == q.size());
    
    constexpr double binwidth_theta_1 = 2., binwidth_theta_2 = 2.;
    constexpr int theta_max = 360, theta_min = 0;
    constexpr int n_bin_theta1 = (theta_max - theta_min) / binwidth_theta_1,
                  n_bin_theta2 = (theta_max - theta_min) / binwidth_theta_2;

    const double PI = TMath::Pi();
    constexpr double radius_internal = 45.2,
                     radius_external = 45.2 + 18 * 2.1;

    TH2F hp_xy("hp_xy", "hp_xy", n_bin_theta1, theta_min, theta_max,
               n_bin_theta2, theta_min, theta_max);
    TH2F hpDiag_xy("hpDiag_xy", "hpDiag_xy", n_bin_theta1, theta_min,
                   theta_max, n_bin_theta2, theta_min, theta_max);
/*	TH2F *hc_xy = new TH2F("hc_xy","Track in xy plane",100,-85,85,100,-85,85);
//	TH2F *hcnew_xy = new TH2F("hcnew_xy","Cluster in xy plane AFTER Obertelli"
//                                        "transform",100,-85,85,100,-85,85);
//	TF1* line_xy = new TF1("line_xy","[0] + [1]*x",-85,85); */
    for(unsigned long i=0; i<x.size(); i++){
        //Fill coordinate space histograms for plots
        //hc_xy->Fill(x->at(i),y->at(i),q->at(i));

        //Loop of indices
        for(int j=0; j < n_bin_theta1; j++){
            const double
                theta1 = (j+0.5)*binwidth_theta_1 + theta_min,
                xt1 = radius_internal * TMath::Cos(theta1 * PI / 180.),
                yt1 = radius_internal * TMath::Sin(theta1 * PI / 180.),
                line1 = (yt1 - y.at(i))/(xt1 - x.at(i)),
                line0 = yt1 - xt1 * line1,
                AA = 1 + line1 * line1,
                BB = 2 * line0 * line1,
                CC = line0 * line0 - radius_external * radius_external,
                delta = BB * BB - 4 * AA * CC;
            
            if(delta < 0) continue;  // Skip negative determinant(?)
            
            double xt2 = -(BB + sqrt(delta))/(2 * AA),
                   yt2 = line0 + line1*xt2,
                   theta2 = 0;
            if(xt2 <= 0)	      theta2 = 180 - asin(yt2 / radius_external) * 180 / PI;
            else{
                if     (yt2 >  0) theta2 = asin(yt2 / radius_external) * 180 / PI;
                else if(yt2 <= 0) theta2 = 360 + asin(yt2 / radius_external) * 180 / PI;
            }

            if((xt1*x.at(i) + yt1*y.at(i))>=0 &&
               (xt2*x.at(i) + yt2*y.at(i))>=0 && (xt1*xt2+yt1*yt2)>=0){
                hp_xy.Fill(theta1,theta2);
                if(abs(theta1-theta2) <= 10) hpDiag_xy.Fill(theta1,theta2);
            }
            else{
                if(delta == 0) continue;

                xt2 = (-BB + sqrt(delta))/(2 * AA);
                yt2 = line0 + line1 * xt2;
                if(xt2 <= 0)	      theta2 = 180 - asin(yt2 / radius_external) * 180 / PI;
                else if(xt2 > 0){
                    if(yt2 > 0)	      theta2 = asin(yt2 / radius_external) * 180 / PI;
                    else if(yt2 <=0 ) theta2 = 360 + asin(yt2 / radius_external) * 180 / PI;
                }

                if((xt1*x.at(i) + yt1*y.at(i)) >= 0 &&
                   (xt2*x.at(i) + yt2*y.at(i)) >= 0 && (xt1*xt2+yt1*yt2) >= 0){
                    hp_xy.Fill(theta1, theta2);
                    if(abs(theta1-theta2) <= 10) hpDiag_xy.Fill(theta1,theta2);
                }
            }
        }
    }
    
    double max_xy;
    if(hpDiag_xy.GetMaximum() >= 10) max_xy = hpDiag_xy.GetMaximum();
    else max_xy = hp_xy.GetMaximum();
    
    double maxtheta1 = 0., maxtheta2=0.;
    for(int ii=0; ii < n_bin_theta1; ii++){
        for(int jj=0; jj < n_bin_theta2; jj++){
            if(hp_xy.GetBinContent(ii+1, jj+1) == max_xy){
                maxtheta1 = (ii+0.5)*binwidth_theta_1 + theta_min;
                maxtheta2 = (jj+0.5)*binwidth_theta_2 + theta_min;
                break;
                //cout << "xy: theta max are " << maxtheta1 << " , " << maxtheta2 << endl;
            }
        }
    }

    const double xmax1 = radius_internal * TMath::Cos(maxtheta1 * PI / 180.),
                 ymax1 = radius_internal * TMath::Sin(maxtheta1 * PI / 180.),
                 xmax2 = radius_external * TMath::Cos(maxtheta2 * PI / 180.),
                 ymax2 = radius_external * TMath::Sin(maxtheta2 * PI / 180.);

    // xy PEAK
    const double par1 = (ymax2-ymax1)/(xmax2-xmax1),
                 par0 = (ymax1 - xmax1*par1);

    /*cout<<"xmax1 "<<xmax1<<" ymax1 "<<ymax1<<" xmax2 "<<xmax2<<" ymax2 "<<ymax2<<endl;
    line_xy->SetParameter(0,par0);
    line_xy->SetParameter(1,par1);
    hc_xy->GetListOfFunctions()->Add(line_xy);
	line_xy->SetLineWidth(1);*/
    
    const vector<double> xTemp(x), yTemp(y), qTemp(q);
    for(auto *i:{&x, &y, &q}) i->clear();
    
    int ringsum = 0,        // hits within the first 5 rings
        filter_result = 0;  // Number of valid pads hit
    
    //Selection of x,y points IN the maxmean+/-1 found in Obertelli transform of xy plane
    for(unsigned long i=0; i<xTemp.size(); i++){
        if( (abs(par1*xTemp[i]-yTemp[i]+par0)/sqrt(1+par1*par1))<= 6 &&
            ((xmax1*xTemp[i] + ymax1*yTemp[i]) >= 0) &&
            ((xmax2*xTemp[i] + ymax2*yTemp[i]) >= 0) &&
            ((xmax1*xmax2 + ymax1*ymax2) >= 0)){
            /*cout << "Taken points= " << xTemp[i] << " , " << yTemp[i] << " , " << zTemp[i] << endl;
            hcnew_xy->Fill(xTemp[i],yTemp[i],qTemp[i]);*/
            x_out.push_back(xTemp[i]);
            y_out.push_back(yTemp[i]);
            q_out.push_back(qTemp[i]);
            filter_result++;
            const double r_mm = sqrt(xTemp[i]*xTemp[i]+yTemp[i]*yTemp[i]);
            if(r_mm < (45.2+5*2.1)) ringsum++;
        }
        else{
            x.push_back(xTemp[i]);
            y.push_back(yTemp[i]);
            q.push_back(qTemp[i]);
        }
    }

    vector<bool> ringbooladd(filter_result, ringsum > 2);
    ringbool.insert(ringbool.end(), ringbooladd.begin(), ringbooladd.end());

    /*for(int ip=0; ip<filter_result; ip++){
        if(ringsum>2) ringbool.push_back(true);
        else ringbool.push_back(false);
    }
	c1->Divide(3,1);
	// Coordinate space
	c1->cd(1);
	hc_xy->Draw("colz");
	// Hough space
	c1->cd(2);
	hp_xy->Draw("colz");
	// Coordinate space : New plot
	c1->cd(3);
	hcnew_xy->Draw("colz");

	c1->Update();
*/
    return filter_result;
}

void minosana::Hough_filter(vector<double> &x, vector<double> &y,
                            vector<double> &z, vector<double> &q,
                            vector<double> &x_out, vector<double> &y_out,
                            vector<double> &z_out, vector<double> &q_out) {
    constexpr int nt_xy = 180, nt_xz = 180, nt_yz = 180,
            nr_xy = 45,  nr_xz = 300, nr_yz = 300;
    
    TH2F hp_xy("hp_xy", "hp_xy", nt_xy, 0, 180, nr_xy, -1*nr_xy, nr_xy);
    TH2F hp_xz("hp_xz", "hp_xz", nt_xz, 0, 180, nr_xz, -1*nr_xz, nr_xz);
    TH2F hp_yz("hp_yz", "hp_yz", nt_yz, 0, 180, nr_yz, -1*nr_yz, nr_yz);
    
    /*TH2F *hc_xy = new TH2F("hc_xy","Track in xy plane",100,-85,85,100,-85,85);
    TH2F *hc_xz = new TH2F("hc_xz","Track in xz plane",250,-50,450,100,-85,85);
    TH2F *hc_yz = new TH2F("hc_yz","Track in yz plane",250,-50,450,100,-85,85);

   	TH2F *hcnew_xy = new TH2F("hcnew_xy","Track in xy plane AFTER Hough transform",100,-85,85,100,-85,85);
   	TH2F *hcnew_xz = new TH2F("hcnew_xz","Track in xz plane AFTER Hough transform",250,-50,450,100,-85,85);
   	TH2F *hcnew_yz = new TH2F("hcnew_yz","Track in yz plane AFTER Hough transform",250,-50,450,100,-85,85);

    int npeaks_xy, npeaks_xz, npeaks_yz;

    TF1* line_xy = new TF1("line_xy","[0] + [1]*x",-85,85);
    TF1* line_xz = new TF1("line_xz","[0] + [1]*x",-50,450);
    TF1* line_yz = new TF1("line_yz","[0] + [1]*x",-50,450);*/
    
    constexpr int nt = std::max(nt_xy, std::max(nt_xz, nt_yz));
    //constexpr int nr = std::max(nr_xy, std::max(nr_xz, nr_yz));
    
    auto radius = [](auto &x, auto &y, auto &thet){
        return x * TMath::Cos(thet*TMath::Pi()/180.) +
               y * TMath::Sin(thet*TMath::Pi()/180.);
    };
    
    for(unsigned long i=0; i < x.size(); i++){
        //Fill coordinate space histograms for plots
        /*hc_xy.Fill(x.at(i),y.at(i),q.at(i));
        hc_xz.Fill(z.at(i),x.at(i),q.at(i));
        hc_yz.Fill(z.at(i),y.at(i),q.at(i));*/
        
        //Loop of indices and fill Histograms
        for(int j=0; j < nt; j++){
            //xy
            const double theta_xy = j * 180./nt_xy,
                    rho_xy = radius(x.at(i), y.at(i), theta_xy);
            
            if(abs(theta_xy) < 180. && abs(rho_xy) < nr_xy){
                hp_xy.Fill(theta_xy, rho_xy);
            }
            
            //xz
            const double theta_xz = j * 180./nt_xz,
                    rho_xz = radius(z.at(i), x.at(i), theta_xz);
            
            if(abs(theta_xz) < 180. && abs(rho_xz) < nr_xz){
                hp_xz.Fill(theta_xz, rho_xz);
            }
            
            //yz
            const double theta_yz = j * 180./nt_yz,
                    rho_yz = radius(z.at(i), y.at(i), theta_yz);
            
            if(abs(theta_yz) < 180. && abs(rho_yz) < nr_yz) {
                hp_yz.Fill(theta_yz, rho_yz);
            }
        }
    }
    
    int locmaxx, locmaxy, locmaxz;
    hp_xy.GetMaximumBin(locmaxx, locmaxy, locmaxz);
    
    if(locmaxy > nr_xy) cerr << "ERRRRRR... locmaxy: " << locmaxy << " nr_xy: " << nr_xy << endl;
    const double rmean_xy =     (locmaxy + 0.5) * 2 - nr_xy;
    const double thetamean_xy = (locmaxx + 0.5) * nt_xy/nt;
    
    hp_xz.GetMaximumBin(locmaxx, locmaxy, locmaxz);
    const double rmean_xz =     (locmaxy + 0.5) * 2 - nr_xz;
    const double thetamean_xz = (locmaxx + 0.5) * nt_xz/nt;
    
    hp_yz.GetMaximumBin(locmaxx, locmaxy, locmaxz);
    const double rmean_yz =     (locmaxy + 0.5) * 2 - nr_yz;
    const double thetamean_yz = (locmaxx + 0.5) * nt_yz/nt;
    
    // Selection of x,y,z points COMMON to the 3 maxmean+/-1 found in Hough
    // spaces for xy, xz and yz spaces
    constexpr double bint_xy = 2., bint_xz = 2., bint_yz = 2.,
            binr_xy = 3., binr_xz = 3., binr_yz = 3.;
    
    for(unsigned long i=0; i<x.size(); i++){
        const double r0_xy = radius(x.at(i), y.at(i), thetamean_xy);
        
        double tmin = thetamean_xy - bint_xy;
        double tmax = thetamean_xy + bint_xy;
        if(tmin < 0)   tmin += 180.;
        if(tmax > 180) tmax -= 180.;
        const double rmin_xy = radius(x.at(i), y.at(i), tmin);
        const double rmax_xy = radius(x.at(i), y.at(i), tmax);
        
        double rinf = min( rmean_xy - binr_xy, rmean_xy + binr_xy);
        double rsup = max( rmean_xy - binr_xy, rmean_xy + binr_xy);
        if((r0_xy >= rinf || rmin_xy >= rinf || rmax_xy >= rinf) &&
           (r0_xy <= rsup || rmin_xy <= rsup || rmax_xy <= rsup)){
            const double r0_xz = radius(z.at(i), x.at(i), thetamean_xz);
            
            tmin = thetamean_xz - bint_xz;
            tmax = thetamean_xz + bint_xz;
            if((tmin) < 0)   tmin += 180.;
            if((tmax) > 180) tmax -= 180.;
            const double rmin_xz = radius(z.at(i), x.at(i), tmin);
            const double rmax_xz = radius(z.at(i), x.at(i), tmax);
            
            rinf = min( rmean_xz - binr_xz, rmean_xz + binr_xz);
            rsup = max( rmean_xz - binr_xz, rmean_xz + binr_xz);
            
            if((r0_xz >= rinf || rmin_xz >= rinf || rmax_xz >= rinf) &&
               (r0_xz <= rsup || rmin_xz <= rsup || rmax_xz <= rsup)){
                const double r0_yz = radius(z.at(i), y.at(i), thetamean_yz);
                
                tmin = thetamean_yz - bint_yz;
                tmax = thetamean_yz + bint_yz;
                if((tmin) < 0)   tmin += 180.;
                if((tmax) > 180) tmax -= 180.;
                const double rmin_yz = radius(z.at(i), y.at(i), tmin);
                const double rmax_yz = radius(z.at(i), y.at(i), tmax);
                
                rinf = min( rmean_yz - binr_yz, rmean_yz + binr_yz);
                rsup = max( rmean_yz - binr_yz, rmean_yz + binr_yz);
                
                if((r0_yz >= rinf || rmin_yz >= rinf || rmax_yz >= rinf) &&
                   (r0_yz <= rsup || rmin_yz <= rsup || rmax_yz <= rsup)){
                    // cout << "Taken points= " << x->at(i) << " , " << y->at(i)
                    //      << " , " << z->at(i) << endl;
                    // hcnew_xy->Fill(x->at(i),y->at(i),q->at(i));
                    // hcnew_xz->Fill(z->at(i),x->at(i),q->at(i));
                    // hcnew_yz->Fill(z->at(i),y->at(i),q->at(i));
                    x_out.push_back(x.at(i));
                    y_out.push_back(y.at(i));
                    z_out.push_back(z.at(i));
                    q_out.push_back(q.at(i));
                }
            }
        }
    }
    /*
     c1->Divide(3,3);
     // Coordinate space
     c1->cd(1);
     hc_xy->Draw("colz");
     c1->cd(2);
     hc_xz->Draw("colz");
     c1->cd(3);
     hc_yz->Draw("colz");

     // Hough space
     c1->cd(4);
     hp_xy->Draw("colz");
     c1->cd(5);
     hp_xz->Draw("colz");
     c1->cd(6);
     hp_yz->Draw("colz");

     // Coordinate space : New plots
     c1->cd(7);
     hcnew_xy->Draw("colz");
     c1->cd(8);
     hcnew_xz->Draw("colz");
     c1->cd(9);
     hcnew_yz->Draw("colz");

     c1->Update();*/
}

void minosana::FindStart(vector<double> &pStart, vector<double> &chi,
                         vector<int> &fitstatus, TGraph &grxz, TGraph &gryz){
    /// This function performs an x-z plane (grxz) and y-z plane (gryz) 2d fit,
    /// to get the initial parameters (stored in &pStart)

    TF1 myfit1(Form("fit%i",threadno),
               [](double *x, double *p){ return p[0]+p[1]*x[0];},
               -100,500,2);
    myfit1.SetParameters(0,10);
    fitstatus = {0,0};
    grxz.Fit(&myfit1, "RQMN");
    chi.at(0) = myfit1.GetChisquare();
    pStart.at(0) = myfit1.GetParameter(0);
    pStart.at(1) = myfit1.GetParameter(1);
    gryz.Fit(&myfit1, "RQMN");
    chi.at(1) = myfit1.GetChisquare();
    pStart.at(2) = myfit1.GetParameter(0);
    pStart.at(3) = myfit1.GetParameter(1);
}


void minosana::debug(){
    /// Look at some raw spectrums, which contain minos tracks in xyz
    if(minossingleevent.size() < 5 && filled){
        minossingleevent.emplace_back(
                Form("t%iEvt%lu",threadno,minossingleevent.size()),
                Form("Thread %i, Event %lu",threadno,minossingleevent.size()),
                100,-100,100,100,-100,100);
        for(unsigned long i=0; i<Xpadnew.size();i++) minossingleevent.back().Fill(
                                    Xpadnew.at(i),Ypadnew.at(i));
        minossingleevent.back().GetXaxis()->SetTitle("X [mm]");
        minossingleevent.back().GetYaxis()->SetTitle("Y [mm]");
        minossingleevent.back().SetMarkerStyle(3);
    }
}

vector<double> minosana::rotatesp(double &rot, vector<double> &initialvector){
    /// This function can rotate the results slopes derived from the 2D planes
    assert(initialvector.size() == 4);

    return vector<double> {
            cos(rot)*initialvector.at(0)-sin(rot)*initialvector.at(2),
            cos(rot)*initialvector.at(1)-sin(rot)*initialvector.at(3),
            sin(rot)*initialvector.at(0)+cos(rot)*initialvector.at(2),
            sin(rot)*initialvector.at(1)+cos(rot)*initialvector.at(3)
    };
}

double distancelinepoint(double x, double y, double z, double *p){
    /// Calculation of the distance between line point
    ROOT::Math::XYZVector xp(x,y,z);
    ROOT::Math::XYZVector x0(p[0],p[2],0.);
    ROOT::Math::XYZVector x1(p[0]+1.*p[1],p[2]+1.*p[3],1.);
    ROOT::Math::XYZVector u =(x1-x0).Unit();
    return ((xp-x0).Cross(u)).Mag2();
}

void minosana::vertex(const vector<double> &p, const vector<double> &pp,
                      double &xv, double &yv, double &zv){
    /// Calculates the vertex of two lines (p,pp) [closest distance]
    /// and stores x, y, and z in (xv, yv, zv)
    assert(p.size() == 4 && pp.size() == 4);

    double alpha, beta, A, B, C;
    alpha = (pp.at(1)*(p.at(0)-pp.at(0))+pp.at(3)*(p.at(2)-pp.at(2)))/
            (pow(pp.at(1),2) + pow(pp.at(3),2) + 1);
    beta = (pp.at(1)*p.at(1)+pp.at(3)*p.at(3) + 1)/
           (pow(pp.at(1),2) + pow(pp.at(3),2) + 1);

    A = beta*(pow(pp.at(1),2) + pow(pp.at(3),2) +1) -
        (pp.at(1)*p.at(1) + pp.at(3)*p.at(3) + 1);
    B = (pow(p.at(1),2) + pow(p.at(3),2) + 1) -
        beta*(pp.at(1)*p.at(1) + pp.at(3)*p.at(3) + 1);
    C = beta*(pp.at(1)*(pp.at(0)-p.at(0))+ pp.at(3)*(pp.at(2)-p.at(2))) -
        (p.at(1)*(pp.at(0)-p.at(0))+ p.at(3)*(pp.at(2)-p.at(2)));

    double sol1, solf1, x, y, z, xp, yp, zp;

    sol1 = -(A*alpha+C)/(A*beta + B);
    solf1 = alpha + beta*sol1;

    x = p.at(0) + p.at(1)*sol1;
    y = p.at(2) + p.at(3)*sol1;
    z = sol1;
    xp = pp.at(0) + pp.at(1)*solf1;
    yp = pp.at(2) + pp.at(3)*solf1;
    zp = solf1;

    xv = (x+xp)/2;
    yv = (y+yp)/2;
    zv = (z+zp)/2;
}

double minosana::distancelineline(vector<double> &l1, vector<double> &l2){
    /// Calculates the distance between two lines
    // vector:  (x_0, m_x, y_0, m_y)^T

    assert(l1.size() == 4 && l2.size() == 4);

    // 1st: vector perpendicular to both
    vector<double> n{
      l1[3]-l2[3],
      l2[1]-l1[1],
      l1[1]*l2[3]-l1[3]*l2[1]};
    double n_mag = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

    // 2nd: construct offset of both vectors
    vector<double> r{l1[0]-l2[0], l1[2]-l2[2], 0};

    return std::abs((n[0]*r[0]+n[1]*r[1]+n[2]*r[2])/n_mag);
}