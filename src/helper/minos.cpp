//
// Created by afrotscher on 1/7/19.
//

#include <TF1.h>
#include <TGraph.h>
#include <TMinuit.h>
#include <Math/Vector3D.h>
#include "minos.h"
#include <numeric>
#include "TCanvas.h"
#include "TMinuitMinimizer.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

using std::vector, std::cerr, std::cout, std::endl, std::min, std::max,
      std::accumulate, std::placeholders::_1,std::placeholders::_2,
      std::placeholders::_3, std::placeholders::_4, std::placeholders::_5;

TMinosPass minosana::analyze() {
    /// MINOS 2. Modify with hough-transformation
    //int padsleft = Xpad.size();
    vector<bool> clusterringbool;
    vector<int> clusternbr, clusterpads;
    int trackNbr = 0, trackNbr_FINAL = 0;

    if (filled) {
        for (int iteration = 0; (Xpad.size() > 9 && iteration < 20); iteration++) {
            int filter_result = Obertelli_filter(Xpad, Ypad, Qpad, Xpadnew,
                                                 Ypadnew, Qpadnew,
                                                 clusterringbool);

            vector<int> temp1(filter_result, iteration);
            vector<int> temp2(filter_result, filter_result);
            clusternbr.insert(clusternbr.end(), temp1.begin(), temp1.end());
            clusterpads.insert(clusterpads.end(), temp2.begin(), temp2.end());

            if (filter_result > 10 && clusterringbool.back()) trackNbr++;
        }
    }

    for (int i = 0; i < Xpadnew.size(); i++) {
        fitdata.add(Xpadnew[i], Ypadnew[i], -1E4, -1E4, Qpadnew[i],
                    clusternbr[i], clusterpads[i], 0);
        Zpadnew.push_back(-1E4);
    }

    if(trackNbr <1 || trackNbr > 4)
        return TMinosPass(r_vertex, theta, phi_vertex, trackNbr, trackNbr_FINAL,
                          z_vertex, {}, {}, {}, {});

    /// Bonus for Tracknumber between 1 and 4
    if (!filled) cerr << "Trackno.:" << trackNbr << " but no evts." << endl;

    /// MINOS 3: Fitting the taken pads for Qmax and Ttrig information /
    //padsleft -= Xpadnew.size();
    for (int i = 0; i < minoscalibvalues.size(); i++) {

        double x_mm = minostrackxy.at(i).at(0);
        double y_mm = minostrackxy.at(i).at(1);
        bool fitbool = false;
        int indexfill = 0;

        for (int j = 0; j < Xpadnew.size(); j++) {
            if (abs(Xpadnew[j] - x_mm) < 0.01 && abs(Ypadnew[j] - y_mm) < 0.01) {
                fitbool = true;
                indexfill = j;
                break;
            }
        }
        // Check if new channel is of interest
        if (!fitbool) continue;

        // if so, we read Q(t), and fill vectors w/ t&Q info after fitting E(t)
        TH1F hfit(Form("hfit%i",threadno), Form("hfit%i",threadno), 512, 0, 512);
        for (int k = 0; k < minoscalibvalues.at(i).size(); k++) {
            if (minoscalibvalues.at(i).at(k) >= 0)
                hfit.SetBinContent(hfit.FindBin(minostime.at(i).at(k)),
                                            minoscalibvalues.at(i).at(k) + 250);
        }
        // Fitting the hfit histogram of last ch. if not empty
        if (hfit.GetSumOfWeights() == 0) continue;

        hfit.GetXaxis()->SetRange(0, 510);
        double hfit_max = hfit.GetMaximum();
        double hfit_max_T = hfit.GetMaximumBin();

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
        int fit2DStatus = hfit.Fit(&fit_function, "QN", "", T_min, T_max);

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

    int padsleft2 = Xpadnew.size();
    vector<double> xin, yin, zin, qin, xout, yout, zout, qout;
    int cluster_temp = 0, array_final = 0;

    vector<TGraph> grxz(trackNbr+1), gryz(trackNbr+1);

    for (int i = 0; i < padsleft2; i++) {
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

            for (int k = 0; k < xout.size(); k++) {
                zmax = max(zmax, zout[k]);
                ringtouch.at(max(0, int((sqrt(pow(xout[k], 2.) +
                                           pow(yout[k], 2.)) - 45.2) / 2.1)))++;
            }
            for (auto &j: ringtouch) if (j > 0) ringsum++;
            if (zmax > 290) ringsum = 16; // particles escaping through the back

            // Decide whether we have a track, and add to 2D-plane TGraphs
            if (xout.size() > 10 && ringsum >= 15) {
                trackNbr_FINAL++;
                //cluster1 = cluster_temp; delete if no use
                for (int l = 0; l < xout.size(); l++) {
                    dataresult.add(xout[l], yout[l], zout[l],
                                   qout[l], trackNbr_FINAL, xout.size(), zmax);
                    array_final++;
                    if(trackNbr_FINAL > grxz.size()) continue;
                    grxz.at(trackNbr_FINAL-1).
                     SetPoint(grxz.at(trackNbr_FINAL-1).GetN(),zout[l],xout[l]);
                    gryz.at(trackNbr_FINAL-1).
                     SetPoint(gryz.at(trackNbr_FINAL-1).GetN(),zout[l],yout[l]);
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

    if(trackNbr_FINAL == 0 ||  trackNbr_FINAL > 3){
        return TMinosPass(r_vertex, theta, phi_vertex, trackNbr, trackNbr_FINAL,
                          z_vertex, {}, {}, {}, {});
    }

  // Look at some events
  if(trackNbr_FINAL ==2) debug();

    /// 5. Fitting filtered tracks in 3D (weight by charge, TMinuit) //
    vector<double> pStart_1{0,1,0,1}, pStart_2{0,1,0,1},pStart_3{0,1,0,1},
                   parFit_1(4), err_1(4), parFit_2(4), parFit_3(4), err_2(4),
                   err_3(4), chi1(2), chi2(2), chi3(2), arglist(10),covxy(3);
    vector<int> fitStatus(2);

    int iflag, nvpar, nparx;
    double amin, edm, errdef, chi2res1, chi2res2;
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
            for(int i=0; i<tmr.x_mm.size(); i++){
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

    //minos5.lock();
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

    for (int i = 0; i < parFit_1.size(); i++)
        min.GetParameter(i, parFit_1.at(i), err_1.at(i));
    double temp[4][4];
    min.mnemat(&temp[0][0], 4);
    covxy.at(0) = temp[1][3]; // entry i=0, j=2

    if (trackNbr_FINAL == 1) parFit_2 = {0, 0, 0, 0};
    else {
        FindStart(pStart_2, chi2, fitStatus, grxz.at(1), gryz.at(1));
        min.SetFCN(funcmin(2));
        minimize(pStart_2);
        min.mnemat(&temp[0][0], 4);
        covxy.at(1) = temp[1][3]; // entry i=0, j=2

        for (int i = 0; i < parFit_2.size(); i++)
            min.GetParameter(i, parFit_2.at(i), err_2.at(i));
    }
    if(trackNbr_FINAL > 2){
        FindStart(pStart_3, chi3, fitStatus, grxz.at(2), gryz.at(2));
        min.SetFCN(funcmin(3));
        minimize(pStart_3);
        min.mnemat(&temp[0][0], 4);
        covxy.at(2) = temp[1][3]; // entry i=0, j=2

        for (int i = 0; i < parFit_3.size(); i++)
            min.GetParameter(i, parFit_3.at(i), err_3.at(i));
    }
    else parFit_3 = {0,0,0,0};

    //minos5.unlock();
    //tmr = {}; // reset out-of-class data structure

    //chi2res1 = (chi1.at(0) + chi1.at(1)) / grxz.at(0).GetN();
    //if(grxz.size() > 1) chi2res2 = (chi2.at(0) + chi2.at(1)) / grxz.at(1).GetN();

    //Rotate from MINOS to beamline for all tracks
    double rot = 30*TMath::Pi()/180;
    vector<double> parFit_1r = rotatesp(rot, parFit_1);
    vector<double> parFit_2r = rotatesp(rot, parFit_2);
    vector<double> parFit_3r = rotatesp(rot, parFit_3);

    /// 6. Get x,y,z reaction vertex from fitted parameters and further variables
    vertex(parFit_1r, parFit_2r, x_vertex, y_vertex, z_vertex);

    r_vertex = sqrt(pow(x_vertex, 2) + pow(y_vertex, 2));

    // cos theta = v1*v2/(|v1|*|v2|), v1 = (0,0,1)^T, v2 = (m1,m3,1)^T
    auto cthet = [](vector<double> a){
        assert(a.size() == 4);
        return acos(1/pow(1+pow(a[1],2) + pow(a[3],2),.5)) * 180./TMath::Pi();};

    theta.at(0) = cthet(parFit_1r);
    theta.at(1) = cthet(parFit_2r);
    theta.at(2) = cthet(parFit_3r);

    // Only calculate angles for real tracks
    while(trackNbr_FINAL < theta.size()){
        theta.pop_back();
    }
    sort(theta.begin(), theta.end()); // Sort angles

  // angles in xy-plane for all tracks
    vector<double> lambda2d{
        atan2(parFit_1r.at(3), parFit_1r.at(1))*180/TMath::Pi(),
        atan2(parFit_2r.at(3), parFit_2r.at(1))*180/TMath::Pi(),
        atan2(parFit_3r.at(3), parFit_3r.at(1))*180/TMath::Pi()
    };

    // Transform fomr neg. angles to positive
    for(auto &i:lambda2d) if(i<0) i = 360 +i;
    sort(lambda2d.begin(), lambda2d.end());

    vector<double> lambda2dc{
        lambda2d.at(2)-lambda2d.at(1),
        lambda2d.at(1)-lambda2d.at(0),
        lambda2d.at(0)+360-lambda2d.at(2)
    };

    // Delete largest angle
    lambda2dc.erase(std::max_element(lambda2dc.begin(), lambda2dc.end()));
    sort(lambda2dc.begin(), lambda2dc.end());

    auto phiinter = [](auto p1, auto p2){
        /// This method calculates the relative angles between two protons
        assert(p1.size() ==4 && p2.size() == 4);
        return acos( (p1[1] * p2[1] + p1[3] * p2[3] + 1)/
                    (sqrt(pow(p1[1], 2) + pow(p1[3], 2) + 1) *
                     sqrt(pow(p2[1], 2) + pow(p2[3], 2) + 1))) *
               180. / TMath::Pi();
    };

    // Calculate all interangles between the protons
    // require events to be simulation-like
    const bool simulation_pass = lambda2dc.at(1) > 160-7./9.*lambda2dc.at(0);

    if(trackNbr_FINAL > 1 && simulation_pass){
        phi_vertex.push_back(phiinter(parFit_1r,parFit_2r));
        if(trackNbr_FINAL > 2){
            phi_vertex.push_back(phiinter(parFit_1r,parFit_3r));
            phi_vertex.push_back(phiinter(parFit_3r,parFit_2r));
        }
    }

    if(trackNbr_FINAL != 3) lambda2dc = {}; // Filter out other events

    double radin = 40, radout = 95;
    auto param = [](vector<double> pf, double radin){
        /// This function calculates the t for which the track crosses the TPC border
        assert(pf.size() == 4);
        return -(pf[1]*pf[0]+pf[3]*pf[2])/(pow(pf[1],2)+pow(pf[3],2)) +
               sqrt(pow((pf[1]*pf[0]+pf[3]*pf[2])/(pow(pf[1],2)+pow(pf[3],2)),2)-
          (pow(pf[0],2)+pow(pf[2],2)-pow(radin,2))/(pow(pf[1],2)+pow(pf[3],2)));
    };

    for(int i=1; i<=trackNbr_FINAL; i++){
        chargeweight.push_back(0);
        for(int j=0; j<tmr.n_Cluster.size(); j++){
            if(i == tmr.n_Cluster.at(j)) chargeweight.back() +=
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

        if(dt > 0) chargeweight.back() /= dt;
    }

    vector<double> verticedist;
    if(trackNbr_FINAL == 2)
        verticedist.push_back(distancelineline(parFit_1r, parFit_2r));
    else if(trackNbr_FINAL == 3){
        verticedist.push_back(distancelineline(parFit_1r, parFit_2r));
        verticedist.push_back(distancelineline(parFit_1r, parFit_3r));
        verticedist.push_back(distancelineline(parFit_2r, parFit_3r));
    }

    // Estimate the angular uncertainty in lambda
    auto lambdaerror = [](vector<double> pf, vector<double> err, double covxy){
        /// Calculates the Gaussian error of the lambda angle for each track
        assert(pf.size() == 4);
        return 180./TMath::Pi()/(pow(pf[1],2)+ pow(pf[3],2))*sqrt(
                pow(pf[1]*err[3],2) + pow(-pf[3]*err[1],2) - 2*pf[1]*pf[3]*covxy);
    };

    vector<double> lambdaE{
        lambdaerror(parFit_1r, err_1, covxy.at(0)),
        lambdaerror(parFit_2r, err_2, covxy.at(1)),
        lambdaerror(parFit_3r, err_3, covxy.at(2))
    };
    // Delete all errors not related to a real track
    while(trackNbr_FINAL < lambdaE.size()) lambdaE.pop_back();

    return TMinosPass(r_vertex, theta, phi_vertex, trackNbr, trackNbr_FINAL,
                      z_vertex, lambda2dc, chargeweight, verticedist, lambdaE);
}

int minosana::Obertelli_filter(vector<double> &x, vector<double> &y,
                               vector<double> &q, vector<double> &x_out,
                               vector<double> &y_out, vector<double> &q_out,
                               vector<bool> &ringbool){
    double bint1 = 2., bint2 = 2.;
    int maxt = 360, mint = 0;
    int nt1 = (maxt-mint)/bint1, nt2 = (maxt-mint)/bint1;

    double PI = TMath::Pi();
    double Rint = 45.2, Rext = 45.2 + 18*2.1;
    int filter_result = 0;

    TH2F hp_xy("hp_xy", "hp_xy", nt1, mint, maxt, nt2, mint, maxt);
    TH2F hpDiag_xy("hpDiag_xy", "hpDiag_xy", nt1, mint, maxt, nt2, mint, maxt);
//	TH2F *hc_xy = new TH2F("hc_xy","Track in xy plane",100,-85,85,100,-85,85);
//	TH2F *hcnew_xy = new TH2F("hcnew_xy","Cluster in xy plane AFTER Obertelli"
//                                        "transform",100,-85,85,100,-85,85);

    double max_xy;
//	TF1* line_xy = new TF1("line_xy","[0] + [1]*x",-85,85);

    vector<double> xTemp, yTemp, qTemp;

    double theta1, theta2 = 0, xt1, yt1, xt2, yt2, line0, line1, delta, AA, BB,
           CC, maxtheta1 = 0., maxtheta2=0., xmax1, ymax1, xmax2, ymax2,
           par0, par1, r_mm;
    int ringsum = 0;
    bool maxfound = false;

    for(unsigned int i=0;i<x.size();i++){
        xTemp.push_back(x.at(i));
        yTemp.push_back(y.at(i));
        qTemp.push_back(q.at(i));

        //Fill coordinate space histograms for plots
//		hc_xy->Fill(x->at(i),y->at(i),q->at(i));

        //Loop of indices
        for(int j=0; j<nt1; j++){
            theta1 = (j+0.5)*bint1 + mint;
            xt1 = Rint * TMath::Cos(theta1*PI/180.);
            yt1 = Rint * TMath::Sin(theta1*PI/180.);
            line1 = (yt1 - y.at(i))/(xt1 - x.at(i));
            line0 = yt1 - xt1 * line1;
            AA = 1 + line1*line1;
            BB = 2*line0*line1;
            CC = line0*line0 - Rext*Rext;

            delta = BB*BB - 4*AA*CC;
            if(delta >= 0){
                xt2 = -(BB + sqrt(delta))/(2*AA);
                yt2 = line0 + line1*xt2;
                if(xt2 <= 0)	        theta2 = 180 - asin(yt2/Rext)*180/PI;
                else{
                    if     (yt2 >  0)   theta2 = asin(yt2/Rext)*180/PI;
                    else if(yt2 <= 0)	theta2 = 360 + asin(yt2/Rext)*180/PI;
                }

                //if(yt2>0){theta2 = 180./PI*acos(xt2/Rext);}
                //else{theta2=360. - 180./PI*acos(xt2/Rext);}

                if((xt1*x.at(i) + yt1*y.at(i))>=0 &&
                   (xt2*x.at(i) + yt2*y.at(i))>=0 && (xt1*xt2+yt1*yt2)>=0){
                    hp_xy.Fill(theta1,theta2);
                    if(abs(theta1-theta2) <= 10) hpDiag_xy.Fill(theta1,theta2);
                }
                else{
                    if(delta != 0){
                        xt2 = (-BB + sqrt(delta))/(2*AA);
                        yt2 = line0 + line1*xt2;
                        if(xt2 <= 0)	      theta2 = 180 - asin(yt2/Rext)*180/PI;
                        else if(xt2 > 0){
                            if(yt2 > 0)	      theta2 = asin(yt2/Rext)*180/PI;
                            else if(yt2 <=0 ) theta2 = 360 + asin(yt2/Rext)*180/PI;
                        }
                        //if(yt2>0){theta2 = 180./PI*acos(xt2/Rext);}
                        //else{theta2=360. - 180./PI*acos(xt2/Rext);}
                        if( (xt1*x.at(i) + yt1*y.at(i)) >= 0 &&
                            (xt2*x.at(i) + yt2*y.at(i)) >= 0 &&
                            (xt1*xt2+yt1*yt2) >= 0){
                            hp_xy.Fill(theta1,theta2);
                            if(abs(theta1-theta2) <= 10) hpDiag_xy.Fill(theta1,theta2);
                        }
                    }
                }
            }
        }
    }
    x.clear();
    y.clear();
    q.clear();

    if(hpDiag_xy.GetMaximum() >= 10) max_xy = hpDiag_xy.GetMaximum();
//		cout << "Max taken in diag... withh value=" << max_xy << endl;
    else max_xy = hp_xy.GetMaximum();

    for(int ii=0; ii<nt1; ii++){
        if(maxfound) break;
        for(int jj=0; jj<nt2; jj++){
            if(hp_xy.GetBinContent(ii+1, jj+1) == max_xy){
                maxtheta1 = (ii+0.5)*bint1 + mint;
                maxtheta2 = (jj+0.5)*bint2 + mint;
                maxfound = true;
                //cout << "xy: theta max are " << maxtheta1 << " , " << maxtheta2 << endl;

            }
            if(maxfound) break;
        }
    }

    xmax1 = Rint * TMath::Cos(maxtheta1*PI/180.);
    ymax1 = Rint * TMath::Sin(maxtheta1*PI/180.);
    xmax2 = Rext * TMath::Cos(maxtheta2*PI/180.);
    ymax2 = Rext * TMath::Sin(maxtheta2*PI/180.);

    // xy PEAK
    par1 = (ymax2-ymax1)/(xmax2-xmax1);
    par0 = (ymax1 - xmax1*par1);

    /*cout<<"xmax1 "<<xmax1<<" ymax1 "<<ymax1<<" xmax2 "<<xmax2<<" ymax2 "<<ymax2<<endl;
    line_xy->SetParameter(0,par0);
    line_xy->SetParameter(1,par1);
    hc_xy->GetListOfFunctions()->Add(line_xy);
	line_xy->SetLineWidth(1);*/

    //Selection of x,y points IN the maxmean+/-1 found in Obertelli transform of xy plane
    for(unsigned int i=0; i<xTemp.size(); i++){
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
            r_mm = sqrt(xTemp[i]*xTemp[i]+yTemp[i]*yTemp[i]);
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
    int nt_xy = 180, nt_xz = 180, nt_yz = 180,
        nr_xy = 45,  nr_xz = 300, nr_yz = 300;
    double bint_xy = 2., bint_xz = 2., bint_yz = 2.,
           binr_xy = 3., binr_xz = 3., binr_yz = 3.;
    int nt = nt_xy,nr = nr_xy;
    double PI = TMath::Pi();

    double rho_xy, rho_xz, rho_yz;
    double theta_xy, theta_xz, theta_yz;

    TH2F hp_xy("hp_xy", "hp_xy", nt_xy, 0, 180, nr_xy, -1*nr_xy, nr_xy);
    TH2F hp_xz("hp_xz", "hp_xz", nt_xz, 0, 180, nr_xz, -1*nr_xz, nr_xz);
    TH2F hp_yz("hp_yz", "hp_yz", nt_yz, 0, 180, nr_yz, -1*nr_yz, nr_yz);

    /*TH2F *hc_xy = new TH2F("hc_xy","Track in xy plane",100,-85,85,100,-85,85);
    TH2F *hc_xz = new TH2F("hc_xz","Track in xz plane",250,-50,450,100,-85,85);
    TH2F *hc_yz = new TH2F("hc_yz","Track in yz plane",250,-50,450,100,-85,85);

   	TH2F *hcnew_xy = new TH2F("hcnew_xy","Track in xy plane AFTER Hough transform",100,-85,85,100,-85,85);
   	TH2F *hcnew_xz = new TH2F("hcnew_xz","Track in xz plane AFTER Hough transform",250,-50,450,100,-85,85);
   	TH2F *hcnew_yz = new TH2F("hcnew_yz","Track in yz plane AFTER Hough transform",250,-50,450,100,-85,85);

    	int npeaks_xy, npeaks_xz, npeaks_yz; */

    vector<double> thetapeaks_xy, rpeaks_xy, thetapeaks_xz, rpeaks_xz,
                   thetapeaks_yz, rpeaks_yz;
    double max_xy, max_xz, max_yz, rmean_xy=0, thetamean_xy=0, rmean_xz=0,
           thetamean_xz=0, rmean_yz=0, thetamean_yz=0;
    /*TF1* line_xy = new TF1("line_xy","[0] + [1]*x",-85,85);
    TF1* line_xz = new TF1("line_xz","[0] + [1]*x",-50,450);
    TF1* line_yz = new TF1("line_yz","[0] + [1]*x",-50,450);*/

    double r0_xy, r0_xz, r0_yz, rmin_xy, rmin_xz, rmin_yz,
            rmax_xy, rmax_xz, rmax_yz,
            tmin, tmax,
            rinf, rsup;

    if(nt<nt_xz)nt=nt_xz;
    if(nr<nr_xz)nr=nr_xz;
    if(nt<nt_yz)nt=nt_yz;
    if(nr<nr_yz)nr=nr_yz;

    //cout<<"In filter !!!"<< endl << "size: "<<x.size()<<endl;

    auto radius = [](auto &x, auto &y, auto &thet){
        return x * TMath::Cos(thet*TMath::Pi()/180.) +
               y * TMath::Sin(thet*TMath::Pi()/180.);
    };

    for(unsigned int i=0; i < x.size(); i++){
        //Fill coordinate space histograms for plots
        /*hc_xy.Fill(x.at(i),y.at(i),q.at(i));
        hc_xz.Fill(z.at(i),x.at(i),q.at(i));
        hc_yz.Fill(z.at(i),y.at(i),q.at(i));*/

        //Loop of indices and fill Histograms
        for(int j=0; j < nt; j++){
            //xy
            theta_xy = j * 180./nt_xy;
            rho_xy = radius(x.at(i), y.at(i), theta_xy);

            if(abs(theta_xy) < 180. && abs(rho_xy) < nr_xy){
                //if(i%40==0) cout<<"i="<<i<<" xy "<<rho_xy<<" "<<theta_xy<<endl;
                hp_xy.Fill(theta_xy, rho_xy);
            }

            //xz
            theta_xz = j * 180./nt_xz;
            rho_xz = radius(z.at(i), x.at(i), theta_xz);

            if(abs(theta_xz) < 180. && abs(rho_xz) < nr_xz){
                //if(i%40==0) cout<<"i="<<i<<" xz "<<rho_xz<<" "<<theta_xz<<endl;
                hp_xz.Fill(theta_xz, rho_xz);
            }

            //yz
            theta_yz = j * 180./nt_yz;
            rho_yz = radius(z.at(i), y.at(i), theta_yz);

            if(abs(theta_yz) < 180. && abs(rho_yz) < nr_yz) {
                //if(i%40==0) cout<<"i="<<i<<" yz "<<rho_yz<<" "<<theta_yz<<endl;
                hp_yz.Fill(theta_yz, rho_yz);
            }
        }
    }

    max_xy = hp_xy.GetMaximum();
    max_xz = hp_xz.GetMaximum();
    max_yz = hp_yz.GetMaximum();

    for(int ii=0; ii<nt; ii++){
        for(int jj=0; jj<nr; jj++){
            if(hp_xy.GetBinContent(ii+1, jj+1) == max_xy && jj<nr_xy){
                thetapeaks_xy.push_back((ii+0.5) * nt_xy/nt);
                rpeaks_xy.push_back((jj+0.5) * 2 - nr_xy);
                rmean_xy     += rpeaks_xy.back();
                thetamean_xy += thetapeaks_xy.back();
                // cout << "xy: " << thetapeaks_xy.back() << " , "
                //      << rpeaks_xy.back() << endl;
            }
            if(hp_xz.GetBinContent(ii+1, jj+1) == max_xz){
                thetapeaks_xz.push_back((ii+0.5) * nt_xz/nt);
                rpeaks_xz.push_back((jj+0.5)*2 - nr_xz);
                rmean_xz     += rpeaks_xz.back();
                thetamean_xz += thetapeaks_xz.back();
                // cout << "xz: " << thetapeaks_xz.back() << " , "
                //      << rpeaks_xz.back() << endl;
            }
            if(hp_yz.GetBinContent(ii+1, jj+1) == max_yz){
                thetapeaks_yz.push_back((ii+0.5)*nt_yz/nt);
                rpeaks_yz.push_back((jj+0.5)*2 - nr_yz);
                rmean_yz     += rpeaks_yz.back();
                thetamean_yz += thetapeaks_yz.back();
                // cout << "yz: " << thetapeaks_yz.back() << " , "
                //      << rpeaks_yz.back() << endl;
            }
        }
    }

    /*cout << "Number of max found :::     IN xy = " << rpeaks_xy.size()
         << " ,     IN xz = " << rpeaks_xz.size() << " ,     IN yz = "
         << rpeaks_yz.size() << endl;

    // xy PEAK
    rmean_xy = rmean_xy / rpeaks_xy.size();
    thetamean_xy = thetamean_xy / thetapeaks_xy.size();
    line_xy.SetParameter(0,rmean_xy / (TMath::Sin(thetamean_xy*PI/180)));
    line_xy.SetParameter(1,( -(TMath::Cos(thetamean_xy*PI/180)) /
                               (TMath::Sin(thetamean_xy*PI/180)) ));
    hc_xy.GetListOfFunctions()->Add(line_xy);

    // xz PEAK
    rmean_xz = rmean_xz / rpeaks_xz.size();
    thetamean_xz = thetamean_xz / thetapeaks_xz.size();
    line_xz.SetParameter(0,rmean_xz / (TMath::Sin(thetamean_xz*PI/180)));
    line_xz.SetParameter(1,( -(TMath::Cos(thetamean_xz*PI/180)) /
                                 (TMath::Sin(thetamean_xz*PI/180)) ));
    hc_xz.GetListOfFunctions()->Add(line_xz);

    // yz PEAK
    rmean_yz = rmean_yz / rpeaks_yz.size();
    thetamean_yz = thetamean_yz / thetapeaks_yz.size();
    line_yz.SetParameter(0,rmean_yz / (TMath::Sin(thetamean_yz*PI/180)));
    line_yz.SetParameter(1,( -(TMath::Cos(thetamean_yz*PI/180)) /
                               (TMath::Sin(thetamean_yz*PI/180)) ));
    hc_yz.GetListOfFunctions()->Add(line_yz); */

    rmean_xy = rpeaks_xy[0];
    thetamean_xy = thetapeaks_xy[0];
    rmean_xz = rpeaks_xz[0];
    thetamean_xz = thetapeaks_xz[0];
    rmean_yz = rpeaks_yz[0];
    thetamean_yz = thetapeaks_yz[0];

    /*line_xy.SetLineWidth(1);
    line_xz.SetLineWidth(1);
    line_yz.SetLineWidth(1);*/

    // Selection of x,y,z points COMMON to the 3 maxmean+/-1 found in Hough
    // spaces for xy, xz and yz spaces
    for(unsigned int i=0; i<x.size(); i++){
        r0_xy = radius(x.at(i), y.at(i), thetamean_xy);

        tmin = thetamean_xy - bint_xy;
        tmax = thetamean_xy + bint_xy;
        if((tmin) < 0)   tmin = tmin + 180.;
        if((tmax) > 180) tmax = tmax - 180.;
        rmin_xy = radius(x.at(i), y.at(i), tmin);
        rmax_xy = radius(x.at(i), y.at(i), tmax);

        rinf = min( rmean_xy - binr_xy, rmean_xy + binr_xy);
        rsup = max( rmean_xy - binr_xy, rmean_xy + binr_xy);
        if((r0_xy >= rinf || rmin_xy >= rinf || rmax_xy >= rinf) &&
           (r0_xy <= rsup || rmin_xy <= rsup || rmax_xy <= rsup)){
            r0_xz = radius(z.at(i), x.at(i), thetamean_xz);

            tmin = thetamean_xz - bint_xz;
            tmax = thetamean_xz + bint_xz;
            if((tmin) < 0)   tmin = tmin + 180.;
            if((tmax) > 180) tmax = tmax - 180.;
            rmin_xz = radius(z.at(i), x.at(i), tmin);
            rmax_xz = radius(z.at(i), x.at(i), tmax);

            rinf = min( rmean_xz - binr_xz, rmean_xz + binr_xz);
            rsup = max( rmean_xz - binr_xz, rmean_xz + binr_xz);

            if((r0_xz >= rinf || rmin_xz >= rinf || rmax_xz >= rinf) &&
               (r0_xz <= rsup || rmin_xz <= rsup || rmax_xz <= rsup)){
                r0_yz = radius(z.at(i), y.at(i), thetamean_yz);

                tmin = thetamean_yz - bint_yz;
                tmax = thetamean_yz + bint_yz;
                if((tmin) < 0)   tmin = tmin + 180.;
                if((tmax) > 180) tmax = tmax - 180.;
                rmin_yz = radius(z.at(i), y.at(i), tmin);
                rmax_yz = radius(z.at(i), y.at(i), tmax);

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

     c1->Update();
     */
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
        for(int i=0; i<Xpadnew.size();i++) minossingleevent.back().Fill(
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

    // vector:
    // (x_0, m_x, y_0, m_y)^T

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