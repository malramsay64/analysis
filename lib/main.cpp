
#include "neighbours.h"
#include "frame.h"
#include "input.h"
//#include "analyse.h"
#include "modular.h"
#include <string>
#include <sstream>
#include <algorithm>

#define BOX 0
#define MOL 1


using namespace std;

bool time_structure=false, regio=false, movie=false, moved=false;

// Key frame every <n> frames
int key_rate = 100;

int main(int argc, char *argv[]){
    string in_fname, out_fname;
    bool quench=false, fast=false;
    int step_size = 1, print, key_rate = 100;
    //int reference = BOX;
    
    /* Read arguments
     * Valid arguments are:
     * -i <file>  input filename
     * -q look at quenched configuration
     * -s <n> look at every n steps
     * -f fast, only first and last
     * -m movie, output movie
     * -d distance moved in quench
     * -k <n> key frame rate
     * -r regio data
     * -t time depedent structure
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
                    cerr << "Invalid number " << argv[i+1] << '\n';
                }
                i++;
            }
            else if (strcmp(argv[i],"-m") == 0){
                movie = true;
            }
            else if (strcmp(argv[i],"-d") == 0){
                moved = true;
                quench = true;
            }
            else if (strcmp(argv[i],"-r") == 0){
                regio = 1;
            }
            else if (strcmp(argv[i],"-k") == 0){
                istringstream ss(argv[i+1]);
                if (!(ss >> key_rate)){
                    cerr << "Invalid number " << argv[i+1] << '\n';
                }
                i++;
            }
            else if (strcmp(argv[i],"-t") == 0){
                time_structure = true;
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
    long frames_read = 0;
    long frame_count = -1;
    int data = true;
    vector<Frame *> key_frames;
    Frame  *current_frame = new Frame;
    //Frame *previous_frame;
    
    // Input file
    
    ifstream myfile;
    myfile.open(in_fname.c_str());
    
    
    // File Size
    ifstream inFile(in_fname.c_str());
    filesize = count(istreambuf_iterator<char>(inFile),
                     istreambuf_iterator<char>(), '\n');
    
    while (data && frame_count < num_frames){
        if (num_frames == 0){
            // Read first frame
            read_data(&myfile, current_frame);
            frames_read++;
            
            // Analyse
            print = 1;
            mod_analyse(current_frame, key_frames, print);
            
            // Add key frame
            key_frames.push_back(new Frame(*current_frame));
            //current_frame = new Frame;
            // Fraction of file read
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
            // Read last frame
            update(&myfile, current_frame);
            frames_read++;
            
            // Analyse
            print = 1;
            mod_analyse(current_frame, key_frames, print);
            
            // Delete unneeded frame
        }
        else if (frame_count % key_rate == 0){
            // Read frame
            update(&myfile, current_frame);
            frames_read++;
            
            // Analyse
            print = 0;
            mod_analyse(current_frame, key_frames, print);
            
            // Add key frame
            key_frames.push_back(new Frame(*current_frame));
            //current_frame = new Frame;
        }
        else if (!fast && frame_count % step_size == 0){
            // Read frame
            update(&myfile, current_frame);
            frames_read++;
            
            print = 0;
            mod_analyse(current_frame, key_frames, print);
        }
        else {
            // Move to next frame
            skip_frame(&myfile);
        }
        frame_count++;
    }
    cerr << "Frames Read: " << frames_read << " of " << frame_count << ", " << key_frames.size() << " key frames" << endl;
    
    for (auto &frame: key_frames){
        delete frame;
    }
    
    return 0;
}
