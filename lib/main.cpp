
#include "neighbours.h"
#include "frame.h"
#include "input.h"
#include <iomanip>
#include "parallel.h"
#include "unit_cell.h"
#include "analyse.h"
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>

#define BOX 0
#define MOL 1


using namespace std;

// Key frame every <n> frames
int key_rate = 100;

int main(int argc, char *argv[]){
    string in_fname, out_fname;
    bool quench=false, fast=false;
    int step_size = 1;
    //int reference = BOX;
    
    /* Read arguments
     * Valid arguments are:
     * -i <file>  input filename
     * -q look at quenched configuration
     * -s <n> look at every n steps
     * -f fast, only first and last
     */
    string arguments = "Valid Arguments are:\n -i <file>\t Input file\n -q\t\t Include quenched\n -s <n>\t\t Every nth frame\n -f \t\t Fast, only first and last frame\n";
    if (argc > 2){
        for (int i = 1; i < argc; i++){
            if (strcmp(argv[i],"-i") == 0){
                in_fname = argv[i+1];
                i++;
            }
            else if (strcmp(argv[i],"-q") == 0){
                quench = true;
            }
            else if (strcmp(argv[i],"-f") == 0){
                fast = true;
            }
            else if (strcmp(argv[i],"-s") == 0){
                istringstream ss(argv[i+1]);
                if (!(ss >> step_size)){
                    cerr << "Invalid number " << argv[1] << '\n';
                }
                i++;
            }
            else {
                cout << "Arguments Incorrect!" << endl << arguments << endl;
                return 1;
            }
        }
    }
    else {
        cout << "Arguments Incorrect!" << endl << arguments << endl;
        return 1;
    }
    
    long filesize;
    long num_frames = 0;
    long frame_count = -1;
    int data = true;
    vector<Frame *> key_frames;
    Frame  *current_frame;
    
    // Input file
    
    ifstream myfile;
    myfile.open(in_fname.c_str());
    
    
    // File Size
    ifstream inFile(in_fname.c_str());
    filesize = count(istreambuf_iterator<char>(inFile),
                     istreambuf_iterator<char>(), '\n');
    
    while (data && frame_count < num_frames){
        // Read frame
        current_frame = new Frame;
        read_data(&myfile, current_frame);
        if (num_frames == 0){
            analyse(current_frame, key_frames,1);
            key_frames.push_back(current_frame);
            inFile.seekg(myfile.tellg());
            num_frames = filesize/(filesize - count(istreambuf_iterator<char>(inFile),
                                                    istreambuf_iterator<char>(), '\n'));
            // If not quenced reduce num_frames by 1
            if (!quench){
                num_frames--;
            }
            frame_count = 0;
        }
        else if (frame_count == num_frames-1){
            analyse(current_frame, key_frames,1);
            delete current_frame;
        }
        else {
            if (!fast && frame_count % step_size == 0){
                analyse(current_frame, key_frames);
            }
            delete current_frame;
        }
        frame_count++;
    }
    
    fprintf(stderr, "Frames Read: %ld\n", frame_count);
    
    for (auto &frame: key_frames){
        delete frame;
    }
    
    return 0;
}
