#include "output.h"
#include "neighbours.h"
#include "frame.h"
#include "input.h"
#include <iomanip>
#include "parallel.h"
#include "unit_cell.h"
#include <string>

#define BOX 0
#define MOL 1

using namespace std;

int main(int argc, char *argv[]){
    string in_fname, out_fname;
    int reference = BOX;
    
    /* Read arguments
     Valid arguments are:
     -i  input filename
     X -o  output filename
     X -g  type of ordering parameter
     X -x  reference frame
     */
    if (argc != 3){
        cout << "Arguments incorrect\n";
        cout << "  -i <input filename>\n" << endl;
    }
    if (strcmp(argv[1],"-i") == 0){
            in_fname = argv[2];
    }
    int data = 1;
    int num_frames = 0;
    ifstream myfile;
    myfile.open(in_fname.c_str());
    Frame *init_frame, *current_frame;
    vector<Frame *> frames(0,0);
    // Read first frame
    frames.push_back(new Frame);
    read_data(&myfile, frames.at(num_frames));
    num_frames++;
    init_frame = frames.front();
    // data present is 0
    while (data){
        
        // Read next frame
        try {
            frames.push_back(new Frame);
            read_data(&myfile, frames.at(num_frames));
            num_frames++;
        }
        catch (...) {
            data = 0;
            delete frames.back();
            frames.pop_back();
        }
    }
    
    string order;
    order = "order";
    // Order parameter last 2 frames
    vector<Frame *>::reverse_iterator f;
    for (f = frames.rbegin(); f != frames.rend() && f-frames.rbegin() < 2; f++){
        order_parameter("order", reference, *f);
    }
    
    current_frame = frames.back();
    
    short_range_order(frames.back());
    MSD(frames);
    diffusion(frames);
    rotation(frames);
    cout << "Diffusion Constant: " << fixed << setprecision(4)  << par_diffusion_constant(init_frame, current_frame)*1000 << endl;
    cout << "Average bond fraction: " << fixed << setprecision(4) << bonded_frac(current_frame) << endl;
    cout << "Local Order: " << setprecision(4) << fixed << par_local_order(current_frame) << endl;
    cout << "Global Order: " << setprecision(4) << fixed << par_global_order(current_frame) << endl;
    cout << "Circle Order: " << setprecision(4) << fixed << par_circle_order(current_frame) << endl;

    stats(current_frame);
    
    print_gnuplot(init_frame);
    print_angles(init_frame);
    
    print_gnuplot(current_frame);
    print_angles(current_frame);
    
    angle_list angles;
    angles = like_me(current_frame);
    //ofstream cell;
    //cell.open("cell.dat");
    //angles.print(&cell);
    
    fprintf(stderr, "Frames Read: %i\n", num_frames);
    
    
    vector<Frame *>::iterator i;
    for (i = frames.begin(); i != frames.end(); i++){
        delete *i;
    }
    return 0;
}
