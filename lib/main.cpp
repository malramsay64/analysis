#include "output.h"
#include "neighbours.h"
#include "frame.h"
#include "input.h"
#include <iomanip>
#include "parallel.h"
//#include "mean.h"
#include <string>

#define BOX 0
#define MOL 1

//static int NUM_THREADS = 8;

using namespace std;

int main(int argc, char *argv[]){
    string in_fname, out_fname;
    int k;
    int reference = BOX;
 
    /* Read arguments
       Valid arguments are:
        -i  input filename
        -o  output filename
        -g  type of ordering parameter
        -x  reference frame
     */
    if (argc != 3){
        cout << "Number of arguments incorrect\n";
        cout << "  -i <input filename>\n" << endl;
        //-o <output filename>\n";
        //cout << "  -g <ordering parameter> (COM, particle, site)\n";
        cout << "  -x <reference frame> (box, mol)\n";
        //return 1;
        //in_fname = "1.lammpstrj";
        //out_fname = "out.lammpstrj";
        //reference = 1;
    }
    for (int i = 0; i < (argc-1)/2; i++){
        k = 2*i+1;
        if (strcmp(argv[k],"-i") == 0){
            in_fname = argv[2*i+2];
        }
        /*
        else if (strcmp(argv[k],"-o") == 0){
            out_fname = argv[2*i+2];
        }
        else if (strcmp(argv[k],"-g") == 0){
            if (strcmp(argv[k+1],"COM") == 0){
            }
            else if (strcmp(argv[k+1], "particle") == 0){
            }
            else if (strcmp(argv[k+1],"site") == 0){
            }
            else {
                //cout << "Incorrect order parameter, use one of (COM|particle|site)\n";
                return 1;
            }
        }
        else if (strcmp(argv[k], "-x") == 0){ 
            if (strcmp(argv[k+1],"box") == 0){
                reference = BOX;
                break;
            }
            else if (strcmp(argv[k+1],"mol") == 0){
                reference = MOL;
                break;    
            }
            else {
                //cout << "Incorrect argument for frame of reference use one of (box|mol)\n";
                return 1;
            }
        }
         */
    }
    int data = 1;
    int num_frames = 0;
    ifstream myfile;
    myfile.open(in_fname.c_str());
    Frame *init_frame, *current_frame;
    vector<Frame *> frames;
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
            //cout << e.what() << endl;
            data = 0;
            delete frames.back();
            frames.pop_back();
        }
    }
    //cout << "Find Neighbours" << endl;
    string order;
    order = "order";
    // Order parameter last 2 frames
    vector<Frame *>::reverse_iterator f;
    for (f = frames.rbegin(); f != frames.rend() && f-frames.rbegin() < 2; f++){
        //cout << "Frame begin " << i << endl;
        //par_neigh(*f);
        order_parameter("order", reference, *f);
        //cout << "Frame complete " << i << endl;
    }
    current_frame = frames.at(num_frames-1);
    // output order parameter for last and second last frames
    //cout << "Order Parameter" << endl;
    //cout << "Short Range Order" << endl;
    short_range_order(frames.back());
    MSD(num_frames, frames);
    diffusion(num_frames, frames);
    rotation(num_frames, frames);
    cout << "Diffusion Constant (parallel): " << par_diffusion_constant(init_frame, current_frame) << endl;
    cout << "Average rotational change (parallel): " << par_average_rotation(init_frame, current_frame) *180/PI << endl;
    cout << "Average dist moved (parallel): " << setprecision(4) << par_average_moved(init_frame, current_frame) << endl;
	cout << "Average bond fraction (parallel): " << setprecision(4) << bonded_frac(current_frame) << endl;
    if (reference == BOX){
        cout << "Reference Frame: Box" << endl;
    }
    else if (reference == MOL){
        cout << "Reference Frame: Mol" << endl;
    }
    //cout << "Size: " << init_frame->size() << endl;
    cout << "Threads: " << NUM_THREADS << endl;
	stats(current_frame);
   
    print_gnuplot(init_frame);
    print_angles(init_frame);
    
    print_gnuplot(current_frame); 
    print_angles(current_frame);
    
    angle_list angles;
    //angles = like_me(current_frame);
    //angles.print(); 
    
    fprintf(stderr, "Frames Read: %i\n", num_frames);
    //cout << current_frame->size()  << " " << current_frame->xmin() << endl;    
    /*
    ofstream file;
    file.open("trj_contact/out.lammpstrj", ios::out);
    dyn_queue q = dyn_queue(&current_frame->particles.front());
    particle *p = q.pop();
    while (p){
        set_colour(p, current_frame);
        //cout << p->id << " " << noboolalpha << p->get_traversed() << " " << get_colour(p, current_frame) << endl;
        print(current_frame, &file);
        p = q.pop();
    }
    file.close();
   */
    //print_gnuplot(current_frame); 
    /*
    cout << "Randomise" << endl;
    randomise_orientation(current_frame);
    cout << "order" << endl;
    short_range_order(current_frame);
    
    order="random";
    order_parameter(order, reference, current_frame);
    
    */
    // Free memory
    vector<Frame *>::iterator i;
    for (i = frames.begin(); i != frames.end(); i++){
        delete *i;
    }
    return 0;
}
