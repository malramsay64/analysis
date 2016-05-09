
#include "neighbours.h"
#include "Frame.h"
#include "input.h"
//#include "analyse.h"
#include <string>
#include <sstream>
#include <algorithm>

#define BOX 0
#define MOL 1



bool time_structure=false, regio=false, movie=false, moved=false, m_orient=false;

// Key frame every <n> frames
int key_rate = 100;

int main(int argc, char *argv[]){
    // Read options file
    char opt;
    std::string fname;
    while ((opt = getopt(argc, argv, options.c_str()))!= EOF){
        switch (opt){
            case 'i':
                fname = std::string{optarg};
                break;
            default:
                useage();
                exit (1);
                break;
        }
    }
    // Read data frame
    std::vector<Frame> frames;
    std::ifstream infile{fname};
    bool read = true;
    while (read){
        try {
        Frame f = read_frame(infile);
        frames.push_back(f);
        }
        catch (...){
        read=false;
        break;
        }
    }
    // Compute data frame

    return 0;
}
