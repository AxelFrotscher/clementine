//
// Created by afrotscher on 6/26/18.
//

#pragma once
#include "histogram_cuts.hh"

namespace runinfo{
    const std::vector<int> transsize = {395'267,
                                        356'843,
                                        0,
                                        1'026'943};    // Runs containing this number are transmission
    const std::vector<int> emptysize = {513'225,
                                        0,  // no empty run for this setting
                                        1'265'431,  // no empty run here either
                                        865'630};    // Empty-target measurement
    const std::vector<int> fulldata  = {36'004'149,
                                        8'051'865,
                                        38'803'327,
                                        58'865'437}; // // full physics run
    const std::string prefix = "/d/d02-1/ag_ob/SEASTAR2_DATA/root/";

    const std::vector<std::vector<std::vector<int>>> plasticrange{
            {{480,620},{700,920},{220,330},{270,1510}},
            {{337,520},{450,750},{220,330},{330,1200}},
            {{350,550},{450,800},{215,380},{350,1100}},
            {{375,575},{450,850},{215,400},{350,1200}}
    };

    const std::vector<std::vector<uint>> pidZrange{
        {200,36,46},{200,28,38},{200,29,39},{200,31,41}}; // y boundaries

    const std::vector<std::vector<std::string>> reactionmodes{
        {"111NbP2P","110NbP2P","110MoP3P","111MoP3P","112MoP3P",
         "113MoP3P","112TcP3P","113TcP3P","114TcP3P"},
        {"88AsP2P","89AsP2P","89SeP2P","90SeP2P","89AsP3P","89SeP3P", "90SeP3P"},
        {"93BrP2P","94BrP2P","95BrP2P","94KrP2P","95KrP2P","96KrP2P","97RbP2P",
         "93BrP3P","94BrP3P","95BrP3P","94KrP3P","95KrP3P","96KrP3P","97RbP3P"},
        {"99RbP2P","100RbP2P","100SrP2P","101SrP2P","102SrP2P","102YP2P","103YP2P",
         "99RbP3P","100RbP3P","100SrP3P","101SrP3P","102SrP3P","102YP3P","103YP3P"}};
}

namespace nancy{
    //Variables and constants used for 111Nb->110Zr setting
    // Cut Isotopes for higherorder correction, F3-7, and F8-11

    // center x, centery, radius x, radius y
    const std::vector<std::vector<std::vector<double>>> cutval{ // for corrections we use 85Ge
            {{2.6449,42.0,0.008,0.6}, {2.6437,41.86,0.008,0.5}},
            {{2.6509,34.0,0.006,0.6}, {2.6608,33.26,0.011,0.5}},
            {{2.6700,36.0,0.006,0.5}, {2.6748,35.46,0.010,0.5}}, // 96Kr
            {{2.7040,37.0,0.006,0.5}, {2.714, 36.43,0.010,0.6}}  // 100Rb
    };

    const std::vector<calibpar> hoparame{
        { 2.644,             // F7absF5X,     110Nb
         -1.453E-5,          // F7linF5X
          0.000'252'8,       // F7linF5A
         -0.000'103,         // F7linF3X
          cutval[0][0][0],   // F7absF5X0
          2.647,             // F11absF9X
         -5.871E-5,          // F11linF9X
          0.000'290'5,       // F11linF9A
          0,                 // F11linF11X
          8.052E-5,          // F11linF11A
          cutval[0][1][0]},  // F11absF9X0
        { 2.651,             // F7absF5X +     88Ge
         -1.999E-5,          // F7linF5X +
          7.127E-5,          // F7linF5A +
         -0.000'271'5,       // F7linF3X +
          cutval[1][0][0],   // F7absF5X0
          2.672,             // F11absF9X +
         -1.334E-4,          // F11linF9X +
          0.000'292'2,       // F11linF9A +
          0,                 // F11linF11X+
          1.079E-4,          // F11linF11A+
          cutval[1][1][0]},  // F11absF9X0
        { 2.67,             // F7absF5X        94Se
         -4.079E-5,          // F7linF5X +
          0.000'122'3,       // F7linF5A +
          0.000'102,         // F7linF3X +
          cutval[2][0][0],   // F7absF5X0
          2.683,             // F11absF9X
         -1.069E-4,          // F11linF9X +
          0.000'227'9,       // F11linF9A +
          0,                 // F11linF11X
          5.212E-5,          // F11linF11A +
          cutval[2][1][0]},  // F11absF9X0
        { 2.706,             // F7absF5X  +      100Kr
         -4.296E-5,          // F7linF5X  +
          8.216E-5,          // F7linF5A  +
         -0.000'122,         // F7linF3X  +
          cutval[3][0][0],   // F7absF5X0
          2.716,             // F11absF9X  +
         -8.227E-5,          // F11linF9X  +
         -0.000'107'1,       // F11linF9A  +
          0,                 // F11linF11X +
          2.5E-5,            // F11linF11A +
          cutval[3][1][0]},  // F11absF9X0
          };

    // For the PID-plot ratios are needed. Boundaries for inc and outg. defined
    const std::vector<double> incval111Nb{
            2.7125, // center x
            41.0,  // center y
            0.007, // radius x
            0.6    // radius y
    };

    const std::vector<double> targetval112Nb{2.729, 41.0, 0.01, 0.5};
    const std::vector<double> targetval111Nb{2.707, 40.9, 0.01, 0.5};
    const std::vector<double> targetval110Nb{2.683, 40.95, incval111Nb.at(2),
                                             incval111Nb.at(3)};
    const std::vector<double> targetval109Nb{2.66, 40.9, incval111Nb.at(2),
                                             incval111Nb.at(3)};
    const std::vector<double> targetval111Zr{2.775, 40, 0.01, 0.4};
    const std::vector<double> targetval110Zr{2.75, 40.0, 0.01, 0.5};
    const std::vector<double> targetval109Zr{2.725, 39.9, 0.01, 0.5};
    const std::vector<double> targetval108Zr{2.6992, 39.91, 0.01, 0.4};

    const std::vector<double> incval110Mo{2.6247, 41.9, 0.006, 0.6};
    const std::vector<double> incval110Nb{2.6882, 40.9, 0.006, 0.6};
    const std::vector<double> incval111Mo{2.6482, 41.9, 0.006, 0.6};
    const std::vector<double> incval112Mo{2.672, 41.9, 0.006, 0.6};
    const std::vector<double> incval113Mo{2.696, 41.92, 0.004, 0.6};
    const std::vector<double> incval113Tc{2.6334, 42.86, 0.004, 0.6};
    const std::vector<double> incval112Tc{2.6103, 42.9, 0.004, 0.6};
    const std::vector<double> incval114Tc{2.65676, 42.882, 0.004, 0.5};

    // Second setting
    const std::vector<double> incval90Se{2.651, 34.0, 0.006, 0.6};
    const std::vector<double> targetval89As{2.711, 32.3, 0.012, 0.5};
    const std::vector<double> targetval88Ge{2.759, 31.23, 0.012, 0.5};
    const std::vector<double> incval89Se{2.6225, 34.0, 0.005, 0.6};
    const std::vector<double> targetval88As{2.679, 32.27, 0.01, 0.5};
    const std::vector<double> targetval87Ge{2.7247, 31.244, 0.01, 0.5};
    const std::vector<double> incval88As{2.670, 33.0, 0.006, 0.6};
    const std::vector<double> incval89As{2.701, 33.0, 0.006, 0.6};
    const std::vector<double> targetval87Ga{2.8121, 30.348, 0.01, 0.5};

    // Third setting
    const std::vector<double> incval97Rb{2.6246, 37.0, 0.006, 0.6};
    const std::vector<double> incval94Kr{2.6116, 36.0, 0.006, 0.6};
    const std::vector<double> incval95Kr{2.640, 36.0, 0.006, 0.6};
    const std::vector<double> incval96Kr{2.670, 36.0, 0.006, 0.6};
    const std::vector<double> incval93Br{2.6578, 35.0, 0.007, 0.6};
    const std::vector<double> incval94Br{2.6878, 35.0, 0.007, 0.6};
    const std::vector<double> incval95Br{2.7172, 35.0, 0.008, 0.6};

    const std::vector<double> targetval96Kr{2.6756, 35.42, 0.009,0.6};
    const std::vector<double> targetval95Br{2.7215, 34.41, 0.012,0.6};
    const std::vector<double> targetval93Br{2.6627, 34.30, 0.009,0.6};
    const std::vector<double> targetval92Se{2.7091, 33.39, 0.011,0.6};
    const std::vector<double> targetval94Br{2.6932, 34.43, 0.010,0.6};
    const std::vector<double> targetval93Se{2.7419, 33.44, 0.011,0.6};
    const std::vector<double> targetval94Se{2.7711, 33.49, 0.011,0.6};
    const std::vector<double> targetval91As{2.7622, 32.42, 0.011,0.6};
    const std::vector<double> targetval92As{2.7933, 32.49, 0.011,0.6};
    const std::vector<double> targetval93As{2.8236, 32.54, 0.011,0.6};

    // Fourth Setting
    const std::vector<double> incval102Y {2.6172, 39.00, 0.006, 0.6};
    const std::vector<double> incval103Y {2.6426, 39.00, 0.006, 0.6};
    const std::vector<double> incval100Sr{2.6332, 38.00, 0.006, 0.6};
    const std::vector<double> incval101Sr{2.6587, 38.00, 0.006, 0.6};
    const std::vector<double> incval102Sr{2.6857, 38.00, 0.006, 0.6};
    const std::vector<double> incval99Rb {2.6770, 37.00, 0.006, 0.6};
    const std::vector<double> incval100Rb{2.7046, 37.00, 0.006, 0.6};

    const std::vector<double> targetval101Sr{2.6697,37.438, 0.01, 0.6};
    const std::vector<double> targetval100Rb{2.7145,36.412, 0.01, 0.6};
    const std::vector<double> targetval102Sr{2.6966,37.467, 0.01, 0.6};
    const std::vector<double> targetval101Rb{2.7410,36.499,0.01,0.6};
    const std::vector<double> targetval99Rb{2.6874,36.389,0.011,0.6};
    const std::vector<double> targetval98Kr{2.7333,35.422,0.011,0.6};
    const std::vector<double> targetval99Kr{2.7613,35.485,0.011,0.5};
    const std::vector<double> targetval100Kr{2.7891,35.506,0.011,0.6};
    const std::vector<double> targetval97Br{2.7829,34.503,0.011,0.6};
    const std::vector<double> targetval98Br{2.8115,34.553,0.011,0.6};
}

namespace nancytrans{
    // Variables and constants used for the empty-target runs
    // Cut Isotopes for higherorder correction, F3-7, and F8-11
    const std::vector<std::vector<std::vector<double>>> cutval{
        {{2.6429, 41.91, 0.008, 0.6}, {2.6460, 40.343, 0.010, 0.5}}, //85Ge
        {{2.6509, 34.00, 0.007, 0.6}, {2.6629, 33.354, 0.012, 0.6}}, //90Se
        {{0, 0, 0, 0}, {0, 0, 0, 0}},
        {{2.6876, 38.00, 0.007, 0.6}, {2.6937, 37.554, 0.013, 0.6}} // 102Sr
    };

    const std::vector<calibpar> hoparame{
        { 2.644,            // F7absF5X
         -1.453E-5,         // F7linF5X
          0.000'252'8,      // F7linF5A
         -0.000'103,        // F7linF3X
          cutval[0][0][0],  // F7absF5X0
          2.647,            // F11absF9X
         -5.871E-5,         // F11linF9X
          0.000'290'5,      // F11linF9A
          0,                // F11linF11X
          8.052E-5,         // F11linF11A
          cutval[0][1][0]}, // F11absF9X0
        { 2.651,            // F7absF5X   +// 88Ge
         -3.146E-5,         // F7linF5X   +
          5.905E-5,         // F7linF5A   +
         -0.000'207,        // F7linF3X   +
          cutval[1][0][0],  // F7absF5X0
          2.665,            // F11absF9X  +
         -9.226E-5,         // F11linF9X  +
          7.7079E-5,        // F11linF9A  +
          0,                // F11linF11X +
          6.194E-5,         // F11linF11A +
          cutval[1][1][0]}, // F11absF9X0
        { 2.644,            // F7absF5X
          -1.453E-5,         // F7linF5X
          0.000'252'8,      // F7linF5A
          -0.000'103,        // F7linF3X
          cutval[2][0][0],  // F7absF5X0
          2.647,            // F11absF9X
          -5.871E-5,         // F11linF9X
          0.000'290'5,      // F11linF9A
          0,                // F11linF11X
          8.052E-5,         // F11linF11A
          cutval[2][1][0]}, // F11absF9X0
        { 2.688,            // F7absF5X   + // 100Kr
         -4.472E-5,         // F7linF5X   +
          9.028E-5,         // F7linF5A   +
         -0.000'212'6,      // F7linF3X   +
          cutval[3][0][0],  // F7absF5X0
          2.695,            // F11absF9X  +
         -8.14E-5,          // F11linF9X  +
         -0.000'160'2,      // F11linF9A  +
          0,                // F11linF11X +
          3.16E-6,          // F11linF11A +
          cutval[3][1][0]}, // F11absF9X0
         };

    // For the PID-plot ratios are needed. Boundaries for inc and outg. defined
    // center x centery radius x radius y
    const std::vector<std::vector<double>> incval{
        {2.681, 41.0, 0.009, 0.6},  //110Nb
        {2.6509, 34.00, 0.007, 0.6},
        {0.,0.,0.,0.},   // cut value should match centered beam nuclide
        {2.7333,37.00,0.006,0.6}  // 101Rb
    };

    const std::vector<std::vector<double>> targetval{
        {2.6799, 40.95, 0.012, incval.at(0).at(3)}, //110Nb
        {2.6629, 33.354, 0.012, 0.65},
        {0.,0.,0.,0.},
        {2.7383, 36.571, 0.012, incval.at(3).at(3)}, //101Rb
    };
}

namespace nancyempty{
    // Variables and constants used for the empty-target run 111Nb->110Zr setting
    // Cut Isotopes for higherorder correction, F3-7, and F8-11
    const std::vector<std::vector<std::vector<double>>> cutval{ // for corrections we use 85Ge
        {{2.6429, 41.91, 0.008, 0.6}, {2.6460, 40.343, 0.008, 0.5}},
        {{2.6509, 34.00, 0.008, 0.6}, {2.661, 33.302, 0.012, 0.65}},
        {{2.6706, 36.00, 0.007, 0.6}, {2.7131, 34.18, 0.013, 0.6}},  //96Kr
        {{2.6876, 38.00, 0.007, 0.6}, {2.7319, 36.123, 0.012, 0.6}}, // 102Sr
    };

    const std::vector<calibpar> hoparame{
            { 2.644,            // F7absF5X
             -1.453E-5,         // F7linF5X
              0.000'252'8,      // F7linF5A
             -0.000'103,        // F7linF3X
              cutval[0][0][0],  // F7absF5X0
              2.647,            // F11absF9X
             -5.871E-5,         // F11linF9X
              0.000'290'5,      // F11linF9A
              0,                // F11linF11X
              8.052E-5,         // F11linF11A
              cutval[0][1][0]}, // F11absF9X0
            { 2.651,            // F7absF5X +
             -3.083E-5,         // F7linF5X +
              8.187E-5,        // F7linF5A
             -0.000'275,        // F7linF3X
              cutval[1][0][0],  // F7absF5X0
              2.672,            // F11absF9X +
             -1.375E-4,        // F11linF9X +
              3.195E-4,         // F11linF9A +
              0,                // F11linF11X
              1.055E-4,         // F11linF11A
              cutval[1][1][0]}, // F11absF9X0
            { 2.67,            // F7absF5X +  94Se
             -2.087E-5,         // F7linF5X +
              0.000'106'7,      // F7linF5A +
             -0.000'207,        // F7linF3X +
              cutval[2][0][0],  // F7absF5X0
              2.716,            // F11absF9X +
             -1.03E-4,          // F11linF9X +
             -4.504E-6,         // F11linF9A +
              0,                // F11linF11X+
              6.056E-5,         // F11linF11A+
              cutval[2][1][0]}, // F11absF9X0
            { 2.688,            // F7absF5X  + // 100Kr (102Sr)
             -4.512E-5,         // F7linF5X  +
              8.898E-5,         // F7linF5A  +
             -0.000'208,        // F7linF3X  +
              cutval[3][0][0],  // F7absF5X0
              2.733,            // F11absF9X  +
             -1.043E-4,         // F11linF9X  +
             -0.000'244'5,      // F11linF9A  +
              0,                // F11linF11X
             -1.166E-5,         // F11linF11A +
              cutval[3][1][0]}, // F11absF9X0
    };

    // For the PID-plot ratios are needed. Boundaries for inc and outg. defined
    const std::vector<std::vector<double>> incval{
            {2.706, 41.0, 0.009, 0.5},
            {2.6510, 34.00, 0.008, 0.6},
            {2.717,35.00, 0.007,0.6},
            {2.7328, 37.00,0.006,0.6}
    };

    const std::vector<std::vector<double>> targetval{
            {2.750, 39.27, 0.011, 0.6},
            {2.663, 33.302, 0.012, 0.65},
            {2.76, 33.22, 0.012, incval.at(2).at(3)},
            {2.7759,35.212,0.012,0.6}
    };
}