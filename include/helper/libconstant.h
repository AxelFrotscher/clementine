//
// Created by afrotscher on 6/26/18.
//

#pragma once
#include "histogram_cuts.hh"

namespace runinfo{
    using std::array, std::vector;
    using a7 = array<int,7>;

    constexpr a7 transsize = {395'267,
                              356'843,
                              0,
                              1'026'943, // Runs containing this number are transmission
                              31'716,
                              0,
                              0};
    constexpr a7 emptysize = {513'225,
                              0,  // no empty run for this setting
                              1'265'431,  // no empty run here either
                              865'630,
                              1'045'377,
                              0,
                              1'598'422};          // Empty-target measurement
    constexpr a7 fulldata  = {36'004'149,
                              8'051'865,
                              38'803'327,
                              58'865'437,
                              40'608'365,
                              55'363'036,
                              160'490'664};//192'409'095}; full physics run

    const  vector<vector<vector<int>>> plasticrange{
            {{480, 620},  {700, 920}, {220, 330}, {270, 1510}},
            {{337, 520},  {450, 750}, {220, 330}, {330, 1200}},
            {{350, 550},  {450, 800}, {215, 380}, {350, 1100}},
            {{375, 575},  {450, 850}, {215, 400}, {350, 1200}},
            {{850, 1250}, {450, 800}, {320, 550}, {250, 1500}},
            {{900, 1500}, {450, 900}, {330, 630}, {250, 2250}},
            {{600, 1600}, {450, 800}, {400, 700}, {250, 2400}}
    };

    const std::vector<std::vector<uint>> pidZrange{
        {200,36,46},{200,28,38},{200,29,39},{200,31,41}, {200,21,31}, {200,23,33},
        {200,25,35}}; // y boundaries

    const std::vector<std::vector<std::string>> reactionmodes{
        { "112MoP2P", "112MoP3P"},
         // "110NbP0P", "111NbP2P","110NbP2P","110MoP2P",
         // "111MoP2P", "110MoP3P","111MoP3P", "113MoP3P", "112TcP3P", "113TcP3P",
         // "114TcP3P"
        {"89SeP2P","90SeP2P","89SeP3P","90SeP3P"}, //  "88GeP0P", "89AsP0P"
        {"94KrP2P","96KrP2P","97RbP2P","94KrP3P","96KrP3P","97RbP3P"},
        // "93BrP2P","94BrP2P","95BrP2P","95KrP2P", "93BrP3P","94BrP3P","95BrP3P",
        // "95KrP3P",
        {"100SrP2P","101SrP2P","102SrP2P","102YP2P","103YP2P",
         "100SrP3P","101SrP3P","102SrP3P","102YP3P","103YP3P"},
        {"68FeP2P", "68CoP2P", "69CoP2P", "70CoP2P", "70NiP2P", "71NiP2P",
         "68FeP3P", "68CoP3P", "69CoP3P", "70CoP3P", "70NiP3P", "71NiP3P"}, //"66MnP2P", "67MnP2P", "67FeP2P",
        {"75ZnP2P", "76ZnP2P", "74CuP2P", "75CuP2P", "72NiP2P", "74NiP2P",
         "75ZnP3P", "76ZnP3P", "74CuP3P", "75CuP3P", "72NiP3P",
         "74NiP3P"}, // "73NiP2P", "73NiP3P",
        {"82GeP2P", "83GeP2P", "80GaP2P", "81GaP2P", "82GaP2P", "78ZnP2P", "79ZnP2P",
         "80ZnP2P", "82GeP3P", "83GeP3P", "80GaP3P", "81GaP3P", "82GaP3P", "78ZnP3P",
         "79ZnP3P", "80ZnP3P"}}; // test
}

namespace nancy{
    //Variables and constants used for 111Nb->110Zr setting
    // Cut Isotopes for higherorder correction, F3-7, and F8-11

    // center x, centery, radius x, radius y
    const std::vector<std::vector<std::vector<double>>> cutval{ // for corrections we use 85Ge
            {{2.6449,42.0, .008, .6}, {2.6437,41.86, .008, .5}},
            {{2.6509,34.0, .006, .6}, {2.6608,33.26, .011, .5}},
            {{2.6700,36.0, .006, .5}, {2.6748,35.46, .010, .5}},  // 96Kr
            {{2.7040,37.0, .006, .5}, {2.714, 36.43, .010, .6}},  // 100Rb
            {{2.6374,25.1, .007, .6}, {2.6246,25.62, .020, .6}},  // 66Mn
            {{2.585, 29.0, .007, .6}, {2.585 ,29.21, .012, .6}},  // 75Cu
            {{2.6666,29.95,.007, .6}, {2.6666,30.39, .020, .6}},  // 80Zn
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
        { 2.638,             // F7absF5X   +     66Cr (66Mn)
         -2.432E-5,          // F7linF5X   +
         -2.687E-5,          // F7linF5A   +
         -7.297E-7,          // F7linF3X   +
          cutval[4][0][0],   // F7absF5X0
          2.624,             // F11absF9X  +
         -7.632E-5,          // F11linF9X  +
         -0.000'609'4,       // F11linF9A  +
          0,                 // F11linF11X +
          0*-1.092E-4,       // F11linF11A +
          cutval[4][1][0]},  // F11absF9X0
        { 2.585,             // F7absF5X   +    75Cu (70Fe)
         -2.455E-5,          // F7linF5X   +
          2.277E-5,          // F7linF5A   +
         -0*7.297E-7,          // F7linF3X +
          cutval[5][0][0],   // F7absF5X0
          2.589,             // F11absF9X  +
         -5.673E-5,          // F11linF9X  +
         -2.126E-4,          // F11linF9A  +
          0,                 // F11linF11X
          0*-1.092E-4,       // F11linF11A +
          cutval[5][1][0]},  // F11absF9X0
        { 2.667,             // F7absF5X,   +  78Ni (80Zn)
         -1.922E-5,          // F7linF5X    +
          7.875E-5,          // F7linF5A    +
          0.000'103,         // F7linF3X    +
          cutval[6][0][0],   // F7absF5X0
          2.672,             // F11absF9X   +
         -6.957E-5,          // F11linF9X   +
         -7.475E-4,          // F11linF9A   +
          0,                 // F11linF11X
          0*8.052E-5,        // F11linF11A  +
          cutval[6][1][0]},  // F11absF9X0
    };

    using a10 = const std::array<double ,10>;
    using a4 = const std::array<double ,4>;
    // For the PID-plot ratios are needed. Boundaries for inc and outg. defined
    a10 incval111Nb{
            2.7125, // center x
            41.0,   // center y
            0.009,  // radius x
            0.6,    // radius y
            .4772,  // transmission (p,2p)
            6.731,  // Brho lower acceptance cut (p,2p)
            6.892,  // Brho higher acceptance cut (p,2p)
            .4668,  // transmission (p,3p)
            6.69,   // Brho lower acceptance cut (p,3p)
            6.80    // Brho higher acceptance cut (p,3p)
    };

    // Fifth value is simulated total transmission from F7 through LH2-Target
    // for P,2P. 6th is for P,3P

    a4 targetval112Nb{2.729, 41.0, .01, .5,};
    a4 targetval111Nb{2.707, 40.9, .01, .5,};
    a4 targetval110Nb{2.683, 40.95,.01, .6,};
    a4 targetval109Nb{2.66, 40.9, .01, .5};
    a4 targetval111Zr{2.775, 40,   .01, .4,};
    a4 targetval110Zr{2.75,  40.0, .01, .5,};
    a4 targetval109Zr{2.725, 39.9, .01, .5,};
    a4 targetval108Zr{2.6992,39.91,.01, .4,};

  a10 incval110Mo{2.6247, 41.9, .006, .6, .467, 6.677, 6.844, .461, 6.637, 6.754};
  a10 incval110Nb{2.6882, 40.9, .006, .6, .477, 6.705, 6.868, .467, 6.665, 6.775};
  a10 incval111Mo{2.6482, 41.9, .006, .6, .469, 6.702, 6.867, .459, 6.661, 6.778};
  a10 incval112Mo{2.672,  41.9, .006, .6, .469, 6.727, 6.891, .458, 6.687, 6.800};
  a10 incval113Mo{2.696,  41.92,.004, .6, .467, 6.751, 6.914, .459, 6.711, 6.824};
  a10 incval112Tc{2.6103, 42.9, .004, .6, .461, 6.700, 6.868, .451, 6.661, 6.779};
  a10 incval113Tc{2.6334, 42.86,.004, .6, .461, 6.724, 6.890, .451, 6.685, 6.802};
  a10 incval114Tc{2.6568, 42.88,.004, .5, .458, 6.758, 6.913, .451, 6.709, 6.825};

    // Second setting
    a10 incval90Se{2.651, 34.0, .006, .6, .546, 6.817, 7.076, .543, 6.748, 6.954};
    a10 incval89Se{2.623, 34.0, .005, .6, .545, 6.795, 7.055, .543, 6.725, 6.933};
    a10 incval88As{2.670, 33.0, .006, .6, .552, 6.796, 7.052, .546, 6.724, 6.926};
    a10 incval89As{2.701, 33.0, .006, .6, .545, 6.818, 7.073, .546, 6.747, 6.947};
    a4 incval88Ge{2.7533,32.0, .008, .6};
    a4 targetval87Ga{2.812, 30.35, .01, .5,};
    a4 targetval87Ge{2.725, 31.24, .01, .5,};
    a4 targetval88Ge{2.762, 31.23, .012,.5,};
    a4 targetval88As{2.679, 32.27, .01, .5,};
    a4 targetval89As{2.711, 32.3,  .012,.5,};
    //a4 targetval90Se{2.662, 33.27, .01, .5};

    // Third setting
    a10 incval93Br{2.6578, 35.0, 0.007, 0.6, .536, 6.679, 6.912, .532, 6.622, 6.800};
    a10 incval94Br{2.6878, 35.0, 0.007, 0.6, .533, 6.704, 6.933, .533, 6.647, 6.821};
    a10 incval95Br{2.7172, 35.0, 0.008, 0.6, .531, 6.728, 6.955, .532, 0.000, 6.842};
    a10 incval94Kr{2.6116, 36.0, 0.006, 0.6, .532, 6.677, 6.914, .529, 6.622, 6.806};
    a10 incval95Kr{2.640,  36.0, 0.006, 0.6, .528, 6.701, 6.935, .528, 6.646, 6.827};
    a10 incval96Kr{2.670,  36.0, 0.006, 0.6, .526, 6.726, 6.956, .527, 6.671, 6.848};
    a10 incval97Rb{2.6246, 37.0, 0.006, 0.6, .521, 6.724, 6.959, .523, 6.671, 6.853};

    a4 targetval91As{2.7622, 32.42, .011, .6};
    a4 targetval92As{2.7933, 32.49, .011, .6};
    a4 targetval93As{2.8236, 32.54, .011, .6};
    a4 targetval92Se{2.7091, 33.39, .011, .6};
    a4 targetval93Se{2.7419, 33.44, .011, .6};
    a4 targetval94Se{2.7711, 33.49, .011, .6};
    a4 targetval93Br{2.6627, 34.30, .009, .6};
    a4 targetval94Br{2.6932, 34.43, .010, .6};
    a4 targetval95Br{2.7215, 34.41, .012, .6};
    a4 targetval96Kr{2.6756, 35.42, .009, .6};

    // Fourth Setting
    a10 incval99Rb {2.6770, 37.00, .006, .6, .513, 6.768, 6.985, .495, 6.714, 6.877};
    a10 incval100Rb{2.7046, 37.00, .006, .6, .513, 6.793, 7.007, .496, 6.739, 6.898};
    a10 incval100Sr{2.6332, 38.00, .006, .6, .508, 6.767, 6.987, .491, 6.714, 6.881};
    a10 incval101Sr{2.6587, 38.00, .006, .6, .508, 6.790, 7.008, .492, 6.738, 6.903};
    a10 incval102Sr{2.6857, 38.00, .006, .6, .506, 6.814, 7.029, .495, 6.762, 6.925};
    a10 incval102Y {2.6172, 39.00, .006, .6, .503, 6.789, 7.010, .488, 6.738, 6.908};
    a10 incval103Y {2.6426, 39.00, .006, .6, .501, 6.813, 7.031, .491, 6.762, 6.930};

    a4 targetval97Br {2.7829,34.503, .011, .6};
    a4 targetval98Br {2.8115,34.553, .011, .6};
    a4 targetval98Kr {2.7333,35.422, .011, .6};
    a4 targetval99Kr {2.7613,35.485, .011, .5};
    a4 targetval100Kr{2.7891,35.506, .011, .6};
    a4 targetval99Rb {2.6874,36.389, .011, .6};
    a4 targetval100Rb{2.7145,36.412, .01,  .6};
    a4 targetval101Rb{2.7410,36.499, .01,  .6};
    a4 targetval101Sr{2.6697,37.438, .01,  .6};
    a4 targetval102Sr{2.6966,37.467, .01,  .6};

    // Fifth Setting      , T(P,2P) Brhocut(P,2P),..
    a10 incval66Mn{2.6372, 25.09, .006, .6, .595, 6.629, 6.937, .596, 6.526, 6.772};
    a10 incval67Mn{2.6773, 25.07, .006, .6, .593, 6.654, 6.958, .593, 6.550, 6.793};
    a10 incval67Fe{2.5738, 26.05, .006, .6, .593, 6.632, 6.944, .593, 6.536, 6.786};
    a10 incval68Fe{2.6124, 26.07, .006, .6, .590, 6.656, 6.965, .591, 6.559, 6.807};
    a10 incval68Co{2.5151, 27.03, .006, .6, .589, 6.636, 6.953, .590, 6.544, 6.803};
    a10 incval69Co{2.5524, 27.00, .006, .6, .587, 6.660, 6.974, .587, 6.567, 6.824};
    a10 incval70Co{2.5893, 27.01, .006, .6, .584, 6.683, 6.995, .585, 6.590, 6.845};
    a10 incval70Ni{2.4957, 27.98, .006, .6, .584, 6.664, 6.981, .585, 6.576, 6.838};
    a10 incval71Ni{2.5317, 27.96, .006, .6, .582, 6.686, 7.002, .582, 6.599, 6.859};

    a4 targetval70Co{2.6020, 26.74, .015, .6};
    a4 targetval69Co{2.5652, 26.70, .015, .6};
    a4 targetval69Fe{2.6636, 25.69, .015, .6};
    a4 targetval68Fe{2.6255, 25.65, .015, .6};
    a4 targetval67Fe{2.5887, 25.57, .015, .6};
    a4 targetval68Mn{2.7291, 24.61, .015, .6};
    a4 targetval67Mn{2.6899, 24.55, .015, .6};
    a4 targetval66Mn{2.6515, 24.50, .015, .6};
    a4 targetval66Cr{2.7611, 23.51, .015, .6};
    a4 targetval65Cr{2.7222, 23.44, .015, .6};
    a4 targetval65V {2.8366, 22.45, .015, .6};
    a4 targetval64V {2.8007, 22.34, .015, .6};

    // Sixth setting (70Fe)
    a10 incval75Zn{2.4984, 29.91, .006, .6, .579, 6.597, 6.890, .579, 6.521, 6.763};
    a10 incval76Zn{2.5312, 29.91, .006, .6, .576, 6.621, 6.911, .577, 6.542, 6.782};
    a10 incval74Cu{2.5503, 28.95, .006, .6, .582, 6.594, 6.883, .582, 6.514, 6.750};
    a10 incval75Cu{2.5847, 28.98, .006, .6, .579, 6.619, 6.904, .580, 6.537, 6.771};
    a10 incval72Ni{2.5697, 27.94, .006, .6, .587, 6.567, 6.855, .587, 6.484, 6.715};
    a10 incval73Ni{2.6057, 27.98, .006, .6, .585, 6.592, 6.876, .585, 6.508, 6.736};
    a10 incval74Ni{2.6411, 28.00, .006, .6, .582, 6.617, 6.900, .582, 6.532, 6.760};

    a4 targetval74Cu{2.5643, 29.05, .014, .6};
    a4 targetval75Cu{2.6017, 29.09, .014, .6};
    a4 targetval73Ni{2.6193, 27.98, .014, .6};
    a4 targetval74Ni{2.6546, 28.02, .014, .6};
    a4 targetval72Ni{2.5848, 27.98, .014, .6};
    a4 targetval71Co{2.6436, 26.88, .014, .6};
    a4 targetval72Co{2.6790, 26.95, .014, .6};
    a4 targetval73Co{2.7146, 27.00, .014, .6};
    //a4 targetval70Co{2.6059, 26.87, .014, .6};
    //a4 targetval69Fe{2.6715, 25.76, .014, .6};
    a4 targetval70Fe{2.7074, 25.79, .014, .6};
    a4 targetval71Fe{2.7426, 25.88, .014, .6};
    a4 targetval72Fe{2.7784, 25.92, .014, .6};

    // Seventh Setting (78Ni)
    a10 incval82Ge{2.5612,31.85, .006, .6, .562, 6.723, 6.995, .561, 6.647, 6.870};
    a10 incval83Ge{2.593, 31.91, .006, .6, .559, 6.745, 7.016, .560, 6.671, 6.891};
    a10 incval80Ga{2.580, 30.91, .006, .6, .567, 6.698, 6.971, .565, 6.620, 6.840};
    a10 incval81Ga{2.6126,30.93, .006, .6, .565, 6.721, 6.992, .564, 6.643, 6.861};
    a10 incval82Ga{2.6446,30.93, .006, .6, .563, 6.746, 7.013, .563, 6.668, 6.882};
    a10 incval78Zn{2.5986,29.92, .006, .6, .573, 6.672, 6.945, .566, 6.589, 6.807};
    a10 incval79Zn{2.6329,29.95, .006, .6, .570, 6.696, 6.966, .566, 6.614, 6.828};
    a10 incval80Zn{2.6663,29.96, .006, .6, .569, 6.722, 6.987, .567, 6.640, 6.850};
    a10 incval81Zn{2.7007,29.96, .006, .6, .566, 6.748, 7.008, .565, 6.668, 6.874};
    a10 incval77Cu{2.6538,28.94, .006, .6, .576, 6.671, 6.936, .568, 6.586, 6.798};
    a10 incval78Cu{2.6890,28.97, .006, .6, .573, 6.697, 6.960, .569, 6.613, 6.822};
    a10 incval79Cu{2.7245,29.00, .006, .6, .571, 6.723, 6.982, .569, 6.640, 6.844};

    a4 targetval81Ga{2.6121, 31.46, .014, .6};
    a4 targetval82Ga{2.6454, 31.52, .014, .6};
    a4 targetval79Zn{2.6348, 30.32, .014, .6};
    a4 targetval80Zn{2.6668, 30.38, .014, .6};
    a4 targetval81Zn{2.7000, 30.41, .014, .6};
    a4 targetval77Cu{2.6577, 29.29, .014, .6};
    a4 targetval78Cu{2.6901, 29.30, .014, .6};
    a4 targetval79Cu{2.7237, 29.32, .014, .6};
    a4 targetval80Cu{2.7586, 29.40, .014, .6};
    a4 targetval76Ni{2.7168, 28.17, .014, .6};
    a4 targetval77Ni{2.7500, 28.22, .014, .6};
    a4 targetval78Ni{2.7857, 28.30, .014, .6};
    a4 targetval75Co{2.7777, 27.10, .014, .6};
    a4 targetval76Co{2.8148, 27.20, .014, .6};

    //Cutting the IC values to combat pile-up
    const std::vector<std::vector<double>> iclimit{
        {8000,10000,4000,5700}, //110Nb
        {5000,6500,2000,3000}, //88Ge
        {5800,7300,2500,3500}, //94Se
        {6400,8100,2800,4000}, //100Kr
        {1600,2800,2000,3400}, //66Cr
        {2000,3300,2700,4000}, //70Fe
        {2400,3800,3500,4800}  //78Ni
    };
}

namespace nancytrans{
    // Variables and constants used for the empty-target runs
    // Cut Isotopes for higherorder correction, F3-7, and F8-11
    const std::vector<std::vector<std::vector<double>>> cutval{
        {{2.6429, 41.91, 0.008, 0.6}, {2.6460, 40.343, 0.010, 0.5}}, //85Ge
        {{2.6509, 34.00, 0.007, 0.6}, {2.6629, 33.354, 0.012, 0.6}}, //90Se
        {{0, 0, 0, 0}, {0, 0, 0, 0}},
        {{2.6876, 38.00, 0.007, 0.6}, {2.6937, 37.554, 0.013, 0.6}}, // 102Sr
        {{0, 0, 0, 0}, {0, 0, 0, 0}},
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
        { 2.688,            // F7absF5X   + // 66Cr
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

    using a10 = std::array<double ,10>;
    using a4 =  std::array<double ,4>;
    // For the PID-plot ratios are needed. Boundaries for inc and outg. defined
    // center x centery radius x radius y
    const std::vector<a10> incval{
        {2.681, 41.0, 0.009, 0.6},  //110Nb
        {2.701, 33.00, 0.007, 0.6}, // 89As
        {0.,0.,0.,0.},   // cut value should match centered beam nuclide
        {2.7333,37.00,0.006,0.6},  // 101Rb
        {2.7333,37.00,0.006,0.6}  // 101Rb
    };

    const std::vector<a4> targetval{
        {2.6799, 40.95, 0.012, incval.at(0).at(3)}, //110Nb
        {2.7122, 32.328, 0.012, 0.65}, // 89As
        {0.,0.,0.,0.},
        {2.7383, 36.571, 0.012, incval.at(3).at(3)}, //101Rb
        {2.7383, 36.571, 0.012, incval.at(3).at(3)},
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

    using a10 = std::array<double ,10>;
    using a4 =  std::array<double ,4>;
    // For the PID-plot ratios are needed. Boundaries for inc and outg. defined
    const std::vector<a10> incval{
            {2.681,   41.0, .009, .5},
            {2.6510, 34.00, .008, .6},
            {2.717,  35.00, .007, .6},
            {2.7328, 37.00, .006, .6},
            {2.5981, 26.76, .006, .6},  // 68Fe(p,2p)
            {2.7328, 37.00, .006, .6},
            {2.6693, 30.46, .006, .6}   // 81Ga(p,2p)
    };

    const std::vector<a4> targetval{
            {2.7505, 39.27, 0.011, 0.6},
            {2.663, 33.302, 0.012, 0.65},
            {2.76, 33.22, 0.012, .6},
            {2.7759,35.212,0.012,0.6},
            {2.5425, 26.30, .016, .7},//{2.6000,25.202,0.012,0.6, 1},
            {2.7759,35.212,0.012,0.6},
            {2.5966, 30.86, .014, .62}//{2.6510,29.834,0.012,0.6, 1}
    };
}
