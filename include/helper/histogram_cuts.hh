//
// Created by axel on 05.06.18.
//

#pragma once

#include <mutex>
#include <vector>

struct calibpar{
    calibpar() = default;
    calibpar(double d, double e, double f, double g, double h, double i,
             double j, double k, double l, double m, double n) :
            F7absF5X(d), F7linF5X(e), F7linF5A(f), F7linF3X(g), F7absF5X0(h),
            F11absF9X(i), F11linF9X(j),F11linF9A(k), F11linF11X(l), F11linF11A(m),
            F11absF9X0(n){
    }
    // Nice encapsulation of variables for the correction
    double F7absF5X;
    double F7linF5X;
    double F7linF5A;
    double F7linF3X;
    double F7absF5X0;

    double F11absF9X;
    double F11linF9X;
    double F11linF9A;
    double F11linF11X;
    double F11linF11A;
    double F11absF9X0;
};

// Definition of global variables here. Keep as short as possible!
inline calibpar p1;

double linfit(double *x, double *par);
bool closeness(const std::vector<double> &d, double sigma=0.1);
double slope(const std::vector<double> &x, const std::vector<double> &y);
double constfit(double *x, double *par);