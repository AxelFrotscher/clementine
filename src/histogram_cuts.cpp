//
// Created by axel on 05.06.18.
//

#include "histogram_cuts.hh"
#include "TMath.h"

using namespace std;

const bool closeness(const vector<double> &d, double sigma){
    // Standard deviation for vector d
    double sum = accumulate(d.begin(),d.end(),0.0);
    double mean = sum /d.size();
    vector<double> diff(d.size());

    transform(d.begin(),d.end(), diff.begin(),
              [mean](double x) { return x - mean; });

    double sq_sum= inner_product(diff.begin(),diff.end(),diff.begin(),0);
    return sigma * mean > sqrt(sq_sum / d.size());
}

double linfit(double *x, double *par){
    // Linear fit function
    return par[0] + par[1]*x[0];
}

const double slope(const vector<double> &x, const vector<double> &y){
    // Simple linear regression (explicit)
    if(!(x.size())==y.size())
        __throw_invalid_argument("Argument number mismatch!\n");
    if(x.size()<2) __throw_invalid_argument("Slope points too few!\n");

    const auto n= x.size();
    const auto s_x = accumulate(x.begin(),x.end(), 0.0);
    const auto s_y = accumulate(y.begin(),y.end(), 0.0);
    const auto s_xx = inner_product(x.begin(),x.end(),x.begin(), 0.0);
    const auto s_xy = inner_product(x.begin(),x.end(),y.begin(), 0.0);
    const auto a = (n*s_xy-s_x*s_y)/(n*s_xx-s_x*s_x);
    return a;
}