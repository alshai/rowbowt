#include <iostream>
#include <fstream>
#include "rle_window_array.hpp"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "please provide input and output file names\n" << std::endl;
        exit(1);
    }
    rle_window_arr mps(argv[1]);
    std::ofstream ofs(argv[2]);
    if (!ofs.good()) {
        std::cerr << "error opening file: " << argv[2] << std::endl;
        exit(1);
    }
    mps.serialize(ofs);
    ofs.close();
    return 0;
}
