// kenneth distefano, Feb. 2024
// FILE: header_rk_v11.h (copied from version 9)
// PURPOSE: the main purpose of this version is to generalize the square lattice 
//          into a rectangular lattice to investigate the effects of the stable 
//          region on a varying size unstable region (non-homogeneous LV model)

#ifndef HEADER_RK_H_V11
#define HEADER_RK_H_V11
#include <random>
#include "Site.h"
using namespace std;

/* -------------------------- global constants ---------------------- */
const int LROW = 200;							// total vertical length
const int LCOL_UNSTABLE = 100;				// horiz. len of unstable area
const int LCOL = LROW/2 + LCOL_UNSTABLE;        // total horizontal side length
const int SITES = LROW*LCOL;                    // total # of sites on 2-d lattice
const int K = 3;								// uniform carrying capacity K =/= -1
const int K_LOW = 1;
const int K_HIGH = 3;
const int N_A = 36000;                            // init num of particles for A
const int N_B = 36000;                            // init num of particles for B
// const int N_C = 25;                             // init num of particles for C
// const int N_D = 5;                              // init num of particles for D
const int N_i[NUM_SPECIES] = {N_A, N_B};
// const int N_i[NUM_SPECIES] = {N_A};
// const int N_i[NUM_SPECIES] = {N_A, N_B, N_C};
const int NUM_STEPS = 2000;                 // number of time steps to simulate
const string INIT_TYPE = "random";          // how to initialize chain
const string REACTION_TYPE = "LV";          // kind of reaction, "MAM", "LV", etc.
const string CARRY_CAP_TYPE = "uniform";    // kind of carrying capacity
const float MAM_PROB = .4;
const float MU = .1;                        // spontaneous predator death prob.
const float SIGMA = .2;                     // prey birth prob.
const float LAMBDA = 0.0;						// predation prob.
const bool MAKE_MOVIE = false;              // output to "lat_config_*.dat" data file

// the time at which the inhomogeneous system is no longer diffusively coupled
const int OPEN_SYSTEM_TRIGGER = CARRY_CAP_TYPE=="uniform" ? NUM_STEPS : 200;

/* -------------------------- prototypes ---------------------- */
// desc: generate a random integer between min and max
// pre: pass only integers, min < max
// post: returns single int, min <= int <= max
int gen_rand_int(const int min, const int max, mt19937 & mt);

// desc: Uniform continuous distribution for random numbers. A continuous 
//          random distribution on the range [min, max) with equal probability 
//          throughout the range. The URNG should be real-valued and deliver 
//          number in the range [0, 1).
// post: returns double within [0,1)
double gen_rand_double(mt19937 & mt);

// desc: initializes a 1d chain of integers randomly
// pre: int array passed by reference automatically
// post: 
void initialize_lattice(Site lat[][LCOL], const int max_pop[],
                        const string init_type, mt19937& mt);

// desc: 
// pre: 
// post: 
void get_avail_rand_sites(const Site lat[][LCOL],
                            vector<int> & rand_nums, mt19937 & mt);

// desc: 
// pre: 
// post:
void get_ofname(string & rho, string & lat, const string runID);

// desc: function for testing to print to screen contents of 1d Site array
// pre: Site array, arr_length= side length of lattice, string array name
// post: prints to screen the contents of 1d Site array
void print_Site_lattice(Site lat[][LCOL], const string arr_name);

// desc: print to screen the contents of a int vec
void print_int_vector(const vector<int> & vec, const string vecName);

// desc: 
// pre: 
// post: 
void print_loop_status(const string index_label, const int index,
                        const string num_tabs, const string loop_type);

// desc: testing function to print a line of dashes to screen (for organization)
void print_dashes();

// desc: print to screen the elements of a 1d int array
void print_1d_int_array(const int arr[],const int arrSize,const string arrName);

// desc: funct that moves the species to a neigboring site, then reaction occurs
// pre: REACTION_TYPES include: "MAM", "LV", "SAM"
// post: 
void reaction_hop(Site lat[][LCOL], const int site_index, const string direction, 
                    mt19937 & mt);

// desc: 
// pre: 
// post: 
string reaction_direction(mt19937 & mt, const bool env, const int siteIndex);

// desc: determine the index of nearest neighbor
// pre: siIndex is an int within the range [0, SITES-1)
// post: returns single int within the same range
int get_nn(const int siIndex, const string dir);

// desc: determine the total population for each species of entire lattice
// pre: lat[] must have length of SITES and speciesPop[] must have length of 
//          NUM_SPECIES
// post: speciesPop[] is pass by reference so the population of each species is 
//          stored within said array
void get_num_density(const Site lat[][LCOL], const int desiredK, int speciesPop[]);

// desc: obtain the sum of all the elements within a 1d int array
// post: returns single int
int sum_1d_int_array(const int arr[], const int arrSize);

// desc: return a random int from within the passed vector array
// pre: "index" is passed by reference so it can be used after function call
// post: single int returned
int get_int_from_vec(const vector<int> & vec, int & index, mt19937 & mt);

// desc: pick a species with nonzero population from a single Site. The probability 
//          of choosing a species depends on its local Site population w.r.t. the 
//          total local Site population
// post: returns a single int
int pick_reacting_species(const Site s, mt19937 & mt);

// desc: print to screen the initial conditions for the simulation
void print_initial_conditions(const int currentSeed, const string runID);

// desc: function to exit simulation and to uniformly print to screen
void exiting_via_exit1();

// desc: check to make sure the number of particles placed on lattice Sites 
//          equals the desired number
void check_num_density(const Site lat[][LCOL]);

// desc: print to screen an update relating to the progress of the simulation
void print_simulation_progress(const int timeStep, const int popArr[]);

// desc: print the current system time
void print_sys_time();

// desc: check command-line inputs
// command line example:[user@joe dir]$ ./a.out fileID
void check_command_line_inputs(char *arr[], const int arr_size);

// desc: check to see if any particles are present
int get_total_lattice_pop(const Site lat[][LCOL]);

// desc: print to screen the total population of each species
void print_lattice_pop(const Site lat[][LCOL]);

// desc: function to output the current lattice configuration a specific/uniform way
void output_lattice_config(ostream & out, Site lat[][LCOL], const int mcs);

// desc: function to determine which Site the randomly picked particle is located
int find_Site(const Site lat[][LCOL],const vector<int> & vec,const int randParticle);

// desc: function to output the current number density a specific/uniform way. this 
//          helps with a non-uniform, spatially varying carrying-capacity.
// post: the array, pops[], returns with the values of the total population of each 
//          species regardless of CARRY_CAP_TYPE because the number of attempted 
//          reactions depend of the total number of particles w/n the system.
void output_num_density(ostream & out, Site lat[][LCOL], int pops[], const int mcs);

// desc: modifies boolArr[] of type bool and of length NUM_STEPS+1 b/c passed by 
//      reference. Index of boolArr[] corresponds to a particular time t. Each 
//      element of boolArr will be either true or false denoting if the system is 
//      either open or closed at that particular time. 
void get_environment(bool boolArr[], const int trigger);





#endif
