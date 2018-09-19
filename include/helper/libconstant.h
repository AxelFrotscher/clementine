//
// Created by afrotscher on 6/26/18.
//

#pragma once
#include "histogram_cuts.hh"

namespace runinfo{
    const int transsize = 395'267; // Runs containing this number are transmission
    const int emptysize = 513'225; // Empty-target measurement
    const std::string prefix = "/d/d02-1/ag_ob/SEASTAR2_DATA/root/";
    enum class Reaction{Nb111PPN,   // 111Nb(p,pn)110Nb
                        Nb111PP2N};
}

namespace nancy{
    //Variables and constants used for 111Nb->110Zr setting

    // Cut Isotopes for higherorder correction, F3-7, and F8-11
    const std::vector<std::vector<double>> cutval{ // for corrections we use 85Ge
        {2.6449,       // center x
         42.0,         // center y
         0.008,        // radius x
         0.6     },    // radius y
        {2.6437,       // center x
         41.86,        // center y
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
    const std::vector<double> incval111Nb{
            2.7125, // center x
            41.0,  // center y
            0.007, // radius x
            0.6    // radius y
    };

    const std::vector<double> targetval110Nb{
            2.683, //2.750, // center x
            40.95, //40.0, // center y
            incval111Nb.at(2),
            incval111Nb.at(3)
    };

    const std::vector<double> targetval109Nb{
            2.66, //2.750, // center x
            40.9, //40.0, // center y
            incval111Nb.at(2),
            incval111Nb.at(3)
    };

    const std::vector<double> targetval110Zr{
            2.75, //2.750, // center x
            40.0, //40.0, // center y
            0.01,
            0.5
    };

    const std::vector<double> incval110Nb{
            2.6882, // center x
            40.9,  // center y
            0.006, // radius x
            0.6    // radius y
    };

    const std::vector<double> targetval109Zr{
            2.725, //2.750, // center x
            39.9, //40.0, // center y
            0.01,
            0.5
    };

    const std::vector<double> incval111Mo{
            2.6482, // center x
            41.9,  // center y
            0.006, // radius x
            0.6    // radius y
    };

    const std::vector<double> incval112Mo{
            2.672, // center x
            41.9,  // center y
            0.006, // radius x
            0.6    // radius y
    };

    const std::vector<double> incval113Tc{
            2.6334, // center x
            42.86,  // center y
            0.004, // radius x
            0.6    // radius y
    };

    const std::vector<double> incval112Tc{
            2.6103, // center x
            42.9,  // center y
            0.004, // radius x
            0.6    // radius y
    };

    const std::vector<double> incval114Tc{
            2.65676, // center x
            42.882,  // center y
            0.004, // radius x
            0.5    // radius y
    };

    const std::vector<double> targetval111Nb{
            2.707, //2.750, // center x
            40.9, //40.0, // center y
            0.01,
            0.5
    };

    const std::vector<double> targetval112Nb{
            2.729, //2.750, // center x
            41.0, //40.0, // center y
            0.01,
            0.5
    };

    const std::vector<double> incval110Mo{
            2.6247, // center x
            41.9,  // center y
            0.006, // radius x
            0.6    // radius y
    };

    const std::vector<double> targetval108Zr{
            2.6992, //2.750, // center x
            39.91, //40.0, // center y
            0.01,
            0.4
    };

    const std::vector<double> incval113Mo{
            2.696, // center x
            41.92,  // center y
            0.004, // radius x
            0.6    // radius y
    };

    const std::vector<double> targetval111Zr{
            2.775, //2.750, // center x
            40, //40.0, // center y
            0.01,
            0.4
    };
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
    const std::vector<double> incval{
            2.706, // center x
            41.0,  // center y
            0.009, // radius x
            0.6    // radius y
    };

    const std::vector<double> targetval{
            2.706, // center x
            41.0, // center y
            incval.at(2),
            incval.at(3)
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
    const std::vector<double> incval{
            2.706, // center x
            41.0,  // center y
            0.009, // radius x
            0.5    // radius y
    };

    const std::vector<double> targetval{
            2.750, // center x
            39.27, // center y
            0.011,
            0.6
    };
}