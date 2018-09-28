//
// Created by afrotscher on 6/26/18.
//

#pragma once
#include "histogram_cuts.hh"

namespace runinfo{
    const std::vector<int> transsize = {395'267,
                                        0,
                                        0,
                                        0};    // Runs containing this number are transmission
    const std::vector<int> emptysize = {513'225,
                                        0,
                                        0,
                                        0};    // Empty-target measurement
    const std::vector<int> fulldata  = {36'004'149,
                                        0,
                                        0,
                                        0}; // // full physics run
    const std::string prefix = "/d/d02-1/ag_ob/SEASTAR2_DATA/root/";

    const std::vector<std::vector<std::vector<int>>> plasticrange{
            {{480,620},{700,920},{220,330},{270,1510}},
            {{337,520},{450,750},{220,330},{330,1200}},
            {{0,0},{0,0},{0,0},{0,0}},
            {{0,0},{0,0},{0,0},{0,0}}
    };

    const std::vector<std::vector<uint>> pidZrange{
        {200,36,46},{200,28,38},{0,0,0},{0,0,0}}; // y boundaries

    const std::vector<std::vector<std::string>> reactionmodes{
        {"111NbP2P","110NbP2P","110MoP3P","111MoP3P","112MoP3P",
         "113MoP3P","112TcP3P","113TcP3P","114TcP3P"},
        {"90SeP2P"},
        {},
        {}};
}

namespace nancy{
    //Variables and constants used for 111Nb->110Zr setting
    // Cut Isotopes for higherorder correction, F3-7, and F8-11

    // center x, centery, radius x, radius y
    const std::vector<std::vector<std::vector<double>>> cutval{ // for corrections we use 85Ge
            {{2.6449,42.0,0.008,0.6}, {2.6437,41.86,0.008,0.5}},
            {{2.6509,34.0,0.006,0.6}, {2.6608,33.26,0.011,0.5}},
            {{0,0,0,0},{0,0,0,0}},
            {{0,0,0,0},{0,0,0,0}}
    };

    const std::vector<calibpar> hoparame{
        { 2.644,             // F7absF5X,     110Nb
         -1.453E-5,          // F7linF5X
         -0.000'252'8,       // F7linF5A
          0.000'103,         // F7linF3X
          cutval[0][0][0],   // F7absF5X0
          2.647,             // F11absF9X
         -5.871E-5,          // F11linF9X
          0.000'290'5,       // F11linF9A
          0,                 // F11linF11X
          8.052E-5,          // F11linF11A
          cutval[0][1][0]},  // F11absF9X0
        { 2.651,             // F7absF5X +     88Ge
         -1.951E-5,          // F7linF5X +
         -7.127E-5,          // F7linF5A (*-1) +
          0.000'271'5,         // F7linF3X +
          cutval[1][0][0],   // F7absF5X0
          2.672,             // F11absF9X +
         -1.334E-4,          // F11linF9X +
          0.000'292'2,       // F11linF9A +
          0,                 // F11linF11X+
          1.079E-4,          // F11linF11A+
          cutval[1][1][0]},  // F11absF9X0
        { 2.644,             // F7absF5X        94Se
         -1.453E-5,          // F7linF5X
         -0.000'252'8,       // F7linF5A
          0.000'103,         // F7linF3X
          cutval[2][0][0],   // F7absF5X0
          2.647,             // F11absF9X
         -5.871E-5,          // F11linF9X
          0.000'290'5,       // F11linF9A
          0,                 // F11linF11X
          8.052E-5,          // F11linF11A
          cutval[2][1][0]},  // F11absF9X0
        { 2.644,             // F7absF5X         100Kr
         -1.453E-5,          // F7linF5X
         -0.000'252'8,       // F7linF5A
          0.000'103,         // F7linF3X
          cutval[3][0][0],   // F7absF5X0
          2.647,             // F11absF9X
         -5.871E-5,          // F11linF9X
          0.000'290'5,       // F11linF9A
          0,                 // F11linF11X
          8.052E-5,          // F11linF11A
          cutval[3][1][0]},  // F11absF9X0
          };

    // For the PID-plot ratios are needed. Boundaries for inc and outg. defined
    const std::vector<double> incval111Nb{
            2.7125, // center x
            41.0,  // center y
            0.007, // radius x
            0.6    // radius y
    };

    const std::vector<double> targetval110Nb{2.683, 40.95, incval111Nb.at(2),
                                             incval111Nb.at(3)};
    const std::vector<double> targetval109Nb{2.66, 40.9, incval111Nb.at(2),
                                             incval111Nb.at(3)};
    const std::vector<double> targetval110Zr{2.75, 40.0, 0.01, 0.5};
    const std::vector<double> incval110Nb{2.6882, 40.9, 0.006, 0.6};

    const std::vector<double> targetval109Zr{2.725, 39.9, 0.01, 0.5};

    const std::vector<double> incval111Mo{2.6482, 41.9, 0.006, 0.6};
    const std::vector<double> incval112Mo{2.672, 41.9, 0.006, 0.6};
    const std::vector<double> incval113Tc{2.6334, 42.86, 0.004, 0.6};
    const std::vector<double> incval112Tc{2.6103, 42.9, 0.004, 0.6};
    const std::vector<double> incval114Tc{2.65676, 42.882, 0.004, 0.5};

    const std::vector<double> targetval111Nb{2.707, 40.9, 0.01, 0.5};
    const std::vector<double> targetval112Nb{2.729, 41.0, 0.01, 0.5};

    const std::vector<double> incval110Mo{2.6247, 41.9, 0.006, 0.6};

    const std::vector<double> targetval108Zr{2.6992, 39.91, 0.01, 0.4};

    const std::vector<double> incval113Mo{2.696, 41.92, 0.004, 0.6};

    const std::vector<double> targetval111Zr{2.775, 40, 0.01, 0.4};


    const std::vector<double> incval90Se{2.651, 34.0, 0.006, 0.6};
    const std::vector<double> targetval89As{2.711, 32.3, 0.01, 0.4};
}

namespace nancytrans{
    // Variables and constants used for the empty-target run 111Nb->110Zr setting
    // Cut Isotopes for higherorder correction, F3-7, and F8-11
    const std::vector<std::vector<double>> cutval{ // for corrections we use 85Ge
        {2.6429,       // center x
         41.91,         // center y
         0.008,        // radius x
         0.6     },    // radius y
        {2.6460,       // center x
         40.343,        // center y
         0.008,        // radius x
         0.5     }     // radius y
    };

    const calibpar hoparame(2.644,            // F7absF5X
                            -1.453E-5,        // F7linF5X
                            -0.000'252'8,       // F7linF5A
                            0.000'103,         // F7linF3X
                            cutval[0][0],     // F7absF5X0
                            2.647,            // F11absF9X
                            -5.871E-5,        // F11linF9X
                            0.000'290'5,        // F11linF9A
                            0,                // F11linF11X
                            8.052E-5,         // F11linF11A
                            cutval[1][0]);    // F11absF9X0

    // For the PID-plot ratios are needed. Boundaries for inc and outg. defined
    // center x centery radius x radius y
    const std::vector<std::vector<double>> incval{
        {2.706, 41.0, 0.009, 0.6},
        {0.,0.,0.,0.},
        {0.,0.,0.,0.},
        {0.,0.,0.,0.}
    };

    const std::vector<std::vector<double>> targetval{
        {2.706, 41.0, incval.at(0).at(2), incval.at(0).at(3)},
        {0, 0, incval.at(1).at(2), incval.at(1).at(3)},
        {0, 0, incval.at(2).at(2), incval.at(2).at(3)},
        {0, 0, incval.at(3).at(2), incval.at(3).at(3)},
    };
}

namespace nancyempty{
    // Variables and constants used for the empty-target run 111Nb->110Zr setting
    // Cut Isotopes for higherorder correction, F3-7, and F8-11
    const std::vector<std::vector<double>> cutval{ // for corrections we use 85Ge
            {2.6429,       // center x
             41.91,        // center y
             0.008,        // radius x
             0.6     },    // radius y
            {2.6460,       // center x
             40.343,       // center y
             0.008,        // radius x
             0.5     }     // radius y
    };

    const calibpar hoparame(2.644,            // F7absF5X
                            -1.453E-5,        // F7linF5X
                            -0.000'252'8,     // F7linF5A
                            0.000'103,        // F7linF3X
                            cutval[0][0],     // F7absF5X0
                            2.647,            // F11absF9X
                            -5.871E-5,        // F11linF9X
                            0.000'290'5,      // F11linF9A
                            0,                // F11linF11X
                            8.052E-5,         // F11linF11A
                            cutval[1][0]);    // F11absF9X0

    // For the PID-plot ratios are needed. Boundaries for inc and outg. defined
    const std::vector<std::vector<double>> incval{
            {2.706, 41.0, 0.009, 0.5},
            {0.,0.,0.,0.},
            {0.,0.,0.,0.},
            {0.,0.,0.,0.}
    };

    const std::vector<std::vector<double>> targetval{
            {2.750, 39.27, 0.011, 0.6},
            {0.,0.,0.,0.},
            {0.,0.,0.,0.},
            {0.,0.,0.,0.}
    };
}