/*
kenneth distefano, Feb. 2024
FILE: main_rk_v11.cpp
PURPOSE: 
    - the main purpose of this version is to generalize the square lattice into a 
        rectangular lattice to investigate the effects of the stable region on a 
        varying size unstable region (non-homogeneous LV model)

    - UPDATE (Apr 24, 2025):
        - it is now possible to close off the two regions => system can no longer 
            be diffusively coupled
            - created const int variable OPEN_SYSTEM_TRIGGER
            - created 1d bool array w/n main*.cpp: open_system_arr[]
            - created function get_environment()
            ~ modified functions reaction_direction(), 

    ~ UPDATE (June 17, 2025):
        ~ include more command line arguments for easier parameter sweeps


    - simulate reaction-diffusion processes on a 2d lattice
        - specifically predator-prey lotka-volterra model
        - single and multiple-species annihilation model
        - includes routine to spatially vary the carrying-capacity at random and by 
            splitting lattice into two regions (in half)
UPDATES FROM PREVIOUS VERSION (v9):
    - 
    - new functions:
        - Apr 24, 2025: get_environment()
            - passes int OPEN_SYSTEM_TRIGGER to determine when to close the two 
                regions from each other for inhomogeneous system
            - passes array of type bool by reference so it can be modified
                - arr[MCS] = [true, true, ..., false, false,...]
            - main idea is to first let a transient predator-prey wavefront enter 
                vulnerable region after total extinction to see if vulnerable region 
                is sustained
            - July 4, 2025:
                - found a bug w/n for-loop bounds
                    - changed for-loop bounds from [0,NUM_STEPS) to [0,NUM_STEPS]
                    - on the final step, the interface because open when it was 
                        supposed to be closed causing particles to enter into the vulnerable region on the final time step
    - depreciated functions:
        - 
    - modified functions:
        - initialize_lattice():
            - now passes lat[][LCOL] instead of lat[][L]
            - got rid of passed argument "const int lat_length" b/c global consts
            - got rid of passed argument "const int max_cap" b/c wasn't used
            - updated initalizations L --> LROW, LCOL
            - got rid of "ABABAB" initialization scheme
            
            ~ make sure random initialization scheme works
        
        - get_avail_rand_sites()
            - now passes lat[][LCOL] instead of lat[][L]
            - got rid of "const in lat_length" b/c do need it (global const)
        
        - print_Site_lattice()
            - now passes lat[][LCOL] instead of lat[][L]
            - got rid of "const in lat_length"
        
        - reaction_hop()
            - now passes lat[][LCOL] instead of lat[][L]
            - converted site index to i,j
                - (site index)/ LCOL and (site index) % LCOL
                - and for neighboring site
        
        - get_num_density()
            - now passes lat[][LCOL] instead of lat[][L]
            - updated loop that traverses through lattice i<L-->i<LROW, j<L-->j<LCOL
        
        - check_num_density()
            - now passes lat[][LCOL] instead of lat[][L]
            - updated loop that traverses through lattice i<L-->i<LROW, j<L-->j<LCOL
        
        - get_total_lattice_pop()
            - now passes lat[][LCOL] instead of lat[][L]
            - updated loop that traverses through lattice i<L-->i<LROW, j<L-->j<LCOL
        
        - print_lattice_pop()
            - now passes lat[][LCOL] instead of lat[][L]
            - updated loop that traverses through lattice i<L-->i<LROW, j<L-->j<LCOL
        
        - output_lattice_config()
            - now passes lat[][LCOL] instead of lat[][L]
            - updated loop that traverses through lattice i<L-->i<LROW, j<L-->j<LCOL
            - this function originally copied and pasted code from the overload 
                function from Site.cpp. so i just rewrote it s.t. I use the 
                insertion operator as I originially intended
                    - cout << lattice[i][j]; ==> (A,B;K)
            - Apr 24, 2025: changed type ofstream to ostream to account for passed 
                argument to be cout
        
        - find_Site()
            - now passes lat[][LCOL] instead of lat[][L]
            - updated range loop that loops through available site indices and 
                converts indices to row/col pair
            ~ make sure the function does what it is suppose to do!
        
        - output_num_density()
            - now passes lat[][LCOL] instead of lat[][L]
            - Apr 24, 2025: changed type ofstream to ostream to account for passed 
                argument to be cout

        - get_ofname()
            - naming convention now considered rectangular lattice
            - Apr 25, 2025: now denotes OPEN_SYSTEM_TRIGGER w/n output file name
            - July 28, 2025: no longer includes seed within output file name (too 
                long and obnoxious)
                - removed passed argument runSeed within function call

        - get_nn()
            - now considers rectangular lattice: L --> LCOL 

        - print_initial_conditions()
            - updated with ternary operator if LROW==LCOL
            - April 28, 2025: now outputs OPEN_SYSTEM_TRIGGER range

        - reaction_direction()
            - Apr 24, 2025
            - includes two more passing arguments: randommly picked site index and 
                environment state (open or closed)
            - now checks if environment is closed and if site index is near interface
                - if so it randomly choses between up, down, horizontal
                    - then chooses left or right depending on site w.r.t. interface
                    - NOTE: probability appears to leak out because horizontal is 
                        being chosen more often than not
        
    -notes:
        - within reaction_hop(), reaction types SAM and MAM, may not reflect the 
            desired outcome if carrying-capacity is not uniform. See above for 
            specific modifcations of function
    
    - to change:
        ~ get rid of unnecessary passed args. to functions (passing global constants)
        ~ find similarities between REACTION_TYPEs to create a general reaction 
            functions for predation, birth, coagulation, etc.
            - subsequently change constants MU, SIGMA, LAMBDA to birth rate, death 
                rate, predation rate, coagulation rate, etc.
        ~ do i need print_lattice_pop()?
        ~ modify check_num_density() s.t. it calls get_num_density()
        ~ modify print_lattice_pop() s.t. it calls get_num_density()
    
    - where can I optimize my code?
        - get_avail_rand_sites()
            - this function erases all available lattice Sites, then runs through 
                all Sites again to determine which Sites are occupied
                    => POTENTIAL FIX: keep track of the Site(s) that were affected 
                        by the reaction and only check if those Site(s). All of the 
                        other Sites should not be affected, thus untouched and still 
                        occupied.

        - i perform a lot of for-loops
        - a lot of functions and there may be over lap on functionality
            - ex. get_total_lattice_pop() and get_num_density() both traverse 
                through the lattice and count up all of the particles, but 
                get_total_lattice_pop() returns an int, whereas get_num_density() 
                returns by reference an array containing populations of all species
            - check out setters and getters within my Site class


*/

#include <iostream>
#include <fstream>
#include <vector>
#include "header_rk_v11.h"
using namespace std;

/* -------------------------- declarations ---------------------- */
Site lattice[LROW][LCOL];           // array to hold all Sites
int si_index;                       // site chosen to during mcs
int vec_index;                      // vector index of avail_rand_sites
int num_density[NUM_SPECIES];       // ints to save num den per specie for 1 mcs
int hoppingSpecies;                 // which species to hop/react during mcs
int randP;                          // randomly generated particle
string direction;                   // direction of reaction
ofstream ave_out, lat_out;          // output recorded data
string numDenName, latConfigName;   // name for output data files
random_device rd;                   // non-deterministic generator
int sim_seed = rd();                // generated deterministic seed
mt19937 mt_gen;                     // URNG engine typedef
vector<int> avail_rand_sites;       // vector, indicies of occupied sites
bool open_system_arr[NUM_STEPS+1];  // 1d bool array; system is open/closed at time t

/* -------------------------- main ---------------------- */
int main(int argc, char *argv[])
{
    // make sure the correct number of arguments are passed
    check_command_line_inputs(argv, argc);

    // seed simulation
    mt_gen.seed(sim_seed);
    
    // print intial conditions to screen
    print_initial_conditions(sim_seed, argv[1]);

    // create output file name
    get_ofname(numDenName, latConfigName, argv[1]);
    
    // open output data files
    ave_out.open(numDenName);
    if(MAKE_MOVIE){lat_out.open(latConfigName);}

    // initize number of particles for each species
    for(int i=0; i<NUM_SPECIES; i++){
        num_density[i] = N_i[i];
    }
    
    // place each particle of each species on lattice
    initialize_lattice(lattice, N_i, INIT_TYPE, mt_gen);
    
    // save initial data
    if(MAKE_MOVIE){output_lattice_config(lat_out, lattice, 0);}
    output_num_density(ave_out, lattice, num_density, 0);

    // sanity check: print num den of each species to screen; compare to N_i
    check_num_density(lattice);

    // initialize array to keep track if environment is open or closed
    get_environment(open_system_arr, OPEN_SYSTEM_TRIGGER);
    
    // this loop begins at 1 to avoid miscalculating num_densities[0]
    for(int ti=1; ti<=NUM_STEPS; ti++)
    {            
        // 1 MCS = completetion of this for-loop, looped N_A+ N_B+... times
        for(int n=0; n < sum_1d_int_array(num_density, NUM_SPECIES); n++)
        {
            // determine which sites are occupied
            get_avail_rand_sites(lattice, avail_rand_sites, mt_gen);

            // choose random particle
            randP = gen_rand_int(1, get_total_lattice_pop(lattice), mt_gen);

            // determine which Site the randomly picked particle is at
            si_index = find_Site(lattice, avail_rand_sites, randP);
            
            // determine right or left reaction
            direction = reaction_direction(mt_gen, open_system_arr[ti], si_index);

            // react
            reaction_hop(lattice, si_index, direction, mt_gen);

        } // end of simulation for-loop

        // output current lattice configuration to create movie later
        if(MAKE_MOVIE){output_lattice_config(lat_out, lattice, ti);}

        // output current population of each species depending on CARRY_CAP_TYPE
        output_num_density(ave_out, lattice, num_density, ti);

        // output status of simulation to screen (should output every tenth timestep)
        // throws an error when NUM_STEPS is less than 10
        if(!(ti % (NUM_STEPS/10))){
            print_simulation_progress(ti, num_density);
        }

    
    } // end of timestep for-loop

    // close output data files
    ave_out.close();
    if(MAKE_MOVIE){lat_out.close();}
    
    
    // exit code
    return 0;
}
