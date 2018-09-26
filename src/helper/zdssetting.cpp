//
// Created by afrotscher on 9/26/18.
//

#include <iostream>
#include "zdssetting.h"

using std::vector;

void setting::loadnumbers(int i) {
    // Loading the right numbers for th right setting
    std::cout << "Setting number " << i << std::endl;
    switch(i){
        case 0:{ // 110Nb
            analysedfile = 57; // 57 Index of analysed file (first, offset)
            goodruns = vector<uint>{1,3,4,5,7,8,9,10,11,13,14,19,26,29,32,33,35,
                                    37, 38,40,41,42,46,49,50,51,52,54,57,58,59};
            transmissionrun = vector<uint>{54};
            emptyrun = vector<uint>{52,53};
            break;
        }
        case 1:{ // 88Ge
            analysedfile = 119;
            goodruns = vector<uint>{3,8,9,10,11,12};
            transmissionrun = vector<uint>{118};
            emptyrun = vector<uint>{119};
            break;
        }
        case 2:{ // 94Se
            analysedfile = 131;
            goodruns = vector<uint>{4,5,8,19,21,22,23,24,25,26,27,28,29,30,31,
                                    32,33,34,35,36,39,40,45,46,47,48};
            transmissionrun = vector<uint>{132};
            emptyrun = vector<uint>{151};
            break;
        }
        case 3:{ // 100Kr
            analysedfile = 180;
            goodruns = vector<uint>{4,6,7,8,9,10,11,12,13,14,15,17,19,20,16,21,
                                    22,23,24,25,26,27,28,29,30,31,32,33,34,35,
                                    36,38};
            transmissionrun = vector<uint>{183};
            emptyrun = vector<uint>{182};
            break;
        }
        default: std::__throw_invalid_argument("Chosen setting is not implemented.\n");
    }
}