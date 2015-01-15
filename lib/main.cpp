#include "output.h"
#include "neighbours.h"
#include "frame.h"
#include "input.h"
#include <iomanip>
#include "parallel.h"
#include "unit_cell.h"
#include "analyse.h"
#include <string>
#include <sstream>

#define BOX 0
#define MOL 1

using namespace std;


// Key frame every <n> frames
int key_rate = 100;

int main(int argc, char *argv[]){
    string in_fname, out_fname;
    //int reference = BOX;
    
    /* Read arguments
     Valid arguments are:
     -i  input filename
     X -o  output filename
     X -g  type of ordering parameter
     X -x  reference frame
     */
    if (argc != 3 && argc != 4){
        cout << "Arguments incorrect\n";
        cout << "  -i <input filename>\n" << endl;
    }
    if (argc == 4){
        istringstream ss(argv[3]);
        if (!(ss >> key_rate)){
            cerr << "Invalid number " << argv[1] << '\n';
        }
    }
    if (strcmp(argv[1],"-i") == 0){
            in_fname = argv[2];
    }
    
    int num_frames = 0;
    int data = true;
    vector<Frame *> key_frames;
    Frame  *current_frame;
    
    // Input file
    ifstream myfile;
    myfile.open(in_fname.c_str());
    
    
    while (data){
        // Read frame
        try {
            current_frame = new Frame;
            read_data(&myfile, current_frame);
            if (num_frames % key_rate == 0){
                analyse(current_frame, key_frames, 1);
                key_frames.push_back(current_frame);

            }
            else {
                analyse(current_frame, key_frames);
                delete current_frame;
            }
            num_frames++;
        }
        catch (int e) {
            data = false;
        }
    }
    
    fprintf(stderr, "Frames Read: %i\n", num_frames);
    
    for (auto &frame: key_frames){
        delete frame;
    }
    
    return 0;
}
