//
// Created by afrotscher on 6/26/18.
//

#pragma once
#include "histogram_cuts.hh"

namespace runinfo{
    const int emptysize = 395267; // Runs containing this number are empty
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
                            -0.0002528,       // F7linF5A
                            0.000103,         // F7linF3X
                            cutval[0][0],     // F7absF5X0
                            2.647,            // F11absF9X
                            -5.871E-5,        // F11linF9X
                            0.0002905,        // F11linF9A
                            0,                // F11linF11X
                            8.052E-5,         // F11linF11A
                            cutval[1][0]);    // F11absF9X0

    // For the PID-plot ratios are needed. Boundaries for inc and outg. defined
    const std::vector<double> incval{
            2.7084, // center x
            41.0,  // center y
            0.009, // radius x
            0.5    // radius y
    };

    const std::vector<double> targetval{
            2.750, // center x
            40.0, // center y
            incval.at(2),
            incval.at(3)
    };
}

namespace nancyempty{
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
                            -0.0002528,       // F7linF5A
                            0.000103,         // F7linF3X
                            cutval[0][0],     // F7absF5X0
                            2.647,            // F11absF9X
                            -5.871E-5,        // F11linF9X
                            0.0002905,        // F11linF9A
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
            2.706, // center x
            41.0, // center y
            incval.at(2),
            incval.at(3)
    };
}