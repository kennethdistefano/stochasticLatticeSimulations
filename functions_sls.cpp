// kenneth distefano, Feb. 2024
// FILE: functions_rk_v11.cpp

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>        // shuffle(), find()
#include <ctime>            // ctime()
#include <fstream>          // outputting lat. config.
#include "header_rk_v11.h"
using namespace std;

// function to return a random int
int gen_rand_int(const int min, const int max, mt19937 & mt)
{
    // local vars
    uniform_int_distribution<int> dist(min, max);

    return dist(mt);
}

// function to return a random float
double gen_rand_double(mt19937 & mt){
    // local vars
    uniform_real_distribution<double> dist(0, 1);

    return dist(mt);
}

// function to initialize the lattice
void initialize_lattice(Site lat[][LCOL], const int max_pop[],
                        const string init_type, mt19937& mt)
{
    // local function vars
    int speciesPop[NUM_SPECIES];    //modifiable int arr for species max pop

    // print to screen simulation update
    cout << "initializing lattice ";
    print_sys_time();

    // make sure carrying-capacity is not less than zero
    if(K == -1){
        cout << "\n\tWARNING! K = -1 is reserved!" << endl;
        exiting_via_exit1();
    }

    // initialize maximum population for each species
    for(int i=0; i<NUM_SPECIES; i++){
        speciesPop[i] = max_pop[i];
    }

    // if a uniform carrying capacity is desired
    if(CARRY_CAP_TYPE == "uniform"){
        // change the carrying capacity if K =/= DEFAULT_CAPACITY
        if(K != DEFAULT_CAPACITY){
            for(int i=0; i<LROW; i++){
                for(int j=0; j<LCOL; j++){
                    lat[i][j].setCapacity(K);
                }
            }
        }

        // check to see if max. num. of parts. has been exceeded N_{max}= KL^d
        if(sum_1d_int_array(max_pop, NUM_SPECIES) > K * SITES){
            cout <<"\n\tLOOK OUT! Double check desired initial population for"
                    <<" each species." <<endl;
            cout << "\t\tN_A + N_B + ... + N_m > KL^d\n" << endl;
            cout << "\tThe number of mcs will not reflect the correct number of " 
                    << "particles within simulation." << endl;
            exiting_via_exit1();
        }
    }
    
    // spatially varying carrying capcity: each Site has randomly dist. K
    else if(CARRY_CAP_TYPE == "random"){
        // traverse down lattice and randomly assign K for that Site
        for(int i=0; i<LROW; i++){
            for(int j=0; j<LCOL; j++){
                lat[i][j].setCapacity(gen_rand_int(K_LOW, K_HIGH, mt));
            }
        }

        // check to see if max. num. of parts. has been exceeded N_{max}= K_HIGH*L^d
        if(sum_1d_int_array(max_pop, NUM_SPECIES) > K_HIGH * SITES){
            cout <<"\n\tLOOK OUT! Double check desired initial population for"
                    <<" each species. " << CARRY_CAP_TYPE 
                    << "carrying capacity chosen." <<endl;
            cout << "\t\tN_A + N_B + ... + N_m > K_HIGH * L^d\n" << endl;
            cout << "\tThe number of mcs will not reflect the correct number of " 
                    << "particles within simulation." << endl;
            exiting_via_exit1();
        }
    }
    
    // split lattice with low and high carrying capacity
    else if(CARRY_CAP_TYPE == "half"){
        // traverse through lattice
        for(int i=0; i<LROW; i++){
            for(int j=0; j<LCOL; j++){
                if(j < LROW/2){
                    lat[i][j].setCapacity(K_LOW);
                }
                else{
                    lat[i][j].setCapacity(K_HIGH);
                }
            }
            
        }

        // make sure correct number of preds/prey will be placed on lattice
        if(sum_1d_int_array(max_pop, NUM_SPECIES) > 
            (SITES*K_HIGH - (LROW*LROW*(K_HIGH-K_LOW))/2)){
            cout <<"\n\tLOOK OUT! Double check desired initial population for"
                    <<" each species." <<endl;
            cout<<"\t\tN_A + N_B + ... + N_m > "
                <<"SITES*K_HIGH - LROW*LROW*(K_HIGH-K_LOW)/2\n" << endl;
            cout << "\tThe number of mcs will not reflect the correct number of " 
                    << "particles within simulation." << endl;
            exiting_via_exit1();
        }

    }

    // incorrect carrying capacity type
    else{
        // print warning to screen and exit safely
        cout<<"WARNING! INCORRECTLY defined \"CARRY_CAP_TYPE\" WITHIN "
                <<"\"initialize_chain()\""<< endl;
        cout<<"The incorrect input: "<< CARRY_CAP_TYPE <<endl;
        cout<<"Correct input: uniform, random, or half" << endl;
        exiting_via_exit1();
    }

    // if initialization type is "random"
    if(init_type == "random")
    {
        // local if-statment vars
        int availSiteIndex, availSpeciesIndex;
        int chosenSpecies, chosenSiteIndex, chosenSiteRow, chosenSiteCol;
        vector<int> availSpecies, availSites;
        bool AreSitesAvail = true;
        bool AreSpeciesAvail = true;

        // all sites on lattice are initially unoccupied
        for(int i=0; i<SITES; i++)
        {
            availSites.push_back(i);
        }

        // init the species vector: max species pop has not been reached yet
        for(int i=0; i<NUM_SPECIES; i++){
            // check if the max population of any species is zero
            if(max_pop[i]){
                // initialize available species
                availSpecies.push_back(i);
            }
        }

        // place particles on lattice
        while(AreSitesAvail && AreSpeciesAvail)
        {
            // shuffle indicies of vectors: lattice Sites and avail species
            shuffle(availSites.begin(), availSites.end(), mt);
            shuffle(availSpecies.begin(), availSpecies.end(), mt);

            // generate random int inbetween [0, vector.size() - 1]
            availSiteIndex = gen_rand_int(0, availSites.size() - 1, mt);
            availSpeciesIndex = gen_rand_int(0, availSpecies.size() - 1, mt);

            // define temporary chosen lattice Site
            chosenSiteIndex = availSites.at(availSiteIndex);
            chosenSiteRow = chosenSiteIndex / LCOL;
            chosenSiteCol = chosenSiteIndex % LCOL;
            chosenSpecies = availSpecies.at(availSpeciesIndex);

            // decrease the number of times the chosen species can be selected
            speciesPop[chosenSpecies] -= 1;

            // place chosen species on chosen site
            lat[chosenSiteRow][chosenSiteCol].incrSingSpecPop(1, chosenSpecies);

            // now check to see if Site has reached max carrying capacity
            if(lat[chosenSiteRow][chosenSiteCol].getTotSitePop()
                    == lat[chosenSiteRow][chosenSiteCol].getCapacity())
            {   
                availSites.erase(availSites.begin() + availSiteIndex);
            }

            // now check to see if the max species population has been reached
            if(speciesPop[chosenSpecies] == 0){
                // remove species from vector availSpecies
                availSpecies.erase(remove(availSpecies.begin(), 
                                            availSpecies.end(),
                                            chosenSpecies), availSpecies.end());
            }

            // check to see if there are no more available Sites or species
            if(!availSites.size()){
                AreSitesAvail = false;
            }
            else if(!availSpecies.size()){
                AreSpeciesAvail = false;
            }

        }
    }

    /*
    if initialization type is "standard"
        - standard means that each Site will have at least one species present
        - DEFAULT_CAPACITY must equal NUM_SPECIES
    */
    else if(init_type == "standard")
    {
        // make sure DEFAULT_CAPACITY == NUM_SPECIES
        if(DEFAULT_CAPACITY != NUM_SPECIES){
            cout << "WARNING!! \"standard\" initialization was chosen and "
                    << "DEFAULT_CAPACITY does not equal NUM_SPECIES." << endl;
            exiting_via_exit1();
        }

        // loop through each Site within lattice
        for(int i=0; i<LROW; i++){
            for(int j=0; j<LCOL; j++){
                // loop through each species at each site
                for(int k=0; k<NUM_SPECIES; k++){
                    lat[i][j].setSingleSpeciesPopulation(1, k);
                }
            }
        }
    }

    // improper string for init_type
    else
    {
        // print warning to screen and exit safely
        cout<<"WARNING! INCORRECTLY PASSED \"init_type\" TO "
                <<"\"initialize_chain()\""<< endl;
        cout<<"The incorrect input: "<< init_type <<endl;
        cout<<"Correct input: random or standard" << endl;
        exiting_via_exit1();
    }

    // let user know initialization is complete
    cout << "initialization was completed ";
    print_sys_time();


    return;
}


// function to determine which lattice sites are occupied
void get_avail_rand_sites(const Site lat[][LCOL],
                            vector<int> & rand_nums, mt19937 & mt)
{
	// make sure vector is empty before adding occupied lattice site indicies
    rand_nums.clear();
    
    // determine which sites are occupied
	for(int i=0; i<LROW; i++){
		for(int j=0; j<LCOL; j++){
            // if Site population does not equal 0, then true
            if(lat[i][j].getTotSitePop()){
                // record site index
                rand_nums.push_back(i*LCOL + j);
            }
        }
	}

	return;
}

// function to return the name of the output file for the average number density
void get_ofname(string & rho,string & lat,const int runSeed,const string runID){
    // local vars
    string run_info, Ktype;

    if(CARRY_CAP_TYPE == "uniform"){
        Ktype = "_K" + to_string(K);
    }
    else if(CARRY_CAP_TYPE == "half"){
        Ktype = "_halfK" + to_string(K_LOW) + "-" + to_string(K_HIGH);
    }
    else if(CARRY_CAP_TYPE == "random"){
        Ktype = "_randK" + to_string(K_LOW) + "-" + to_string(K_HIGH);
    }
    // incorrect carrying capacity type
    else{
        // print warning to screen and exit safely
        cout<<"WARNING! INCORRECTLY defined \"CARRY_CAP_TYPE\" WITHIN "
                <<"\"get_ofname()\""<< endl;
        cout<<"The incorrect input: "<< CARRY_CAP_TYPE <<endl;
        cout<<"Correct input: uniform, random, or half" << endl;
        exiting_via_exit1();
    }

    // create string with simulation details
    run_info = "species"+to_string(NUM_SPECIES)+"_L"
                +(LROW==LCOL ? to_string(LROW) : to_string(LROW)+"-"+to_string(LCOL))
                +"_steps"+to_string(NUM_STEPS)+"_Na"+to_string(N_A)
                + Ktype+"_seed"+to_string(runSeed)+"_id"+runID+".dat";
    
    // create output file name for the number density of each run
    rho = "num_density_" + run_info;

    // create output file name for the lattice configuration of each attempt
    lat = "lat_config_" + run_info;
    
    return;
}

// function to print contents of an array to screen
void print_Site_lattice(Site lat[][LCOL], const string arr_name){
    
    cout << "printing "<< arr_name <<"[], arr_size= "<< LROW*LCOL 
            <<endl;
    cout << arr_name <<"[] = "<<endl;
    for(int i=0; i<LROW; i++){
        for(int j=0; j<LCOL; j++){
            cout << lat[i][j] << " ";
        }
        cout<<endl;
    }
    cout << endl;

    return;
}

// function to print the contents of the vector to screen
void print_int_vector(const vector<int> & vec, const string vecName){
    
    cout << vecName<<".size()= "<<vec.size()<< ", ";
    cout << vecName<<"[]= ";
    for(int i : vec){
        cout << i <<" ";
    }
    cout << endl;
    return;
}

// function to print current status of the loop
void print_loop_status(const string index_label, const int index,
                        const string num_tabs, const string loop_type)
{
    // local vars
    cout << num_tabs << "inside "<<loop_type<<", "
            << index_label << "= " << index << endl;

    return;
}

// test function to print line of dashes
void print_dashes(){
    cout << "-----------------------------------------------------------"<<endl;
    return;
}

// print the elements of a 1d int
void print_1d_int_array(const int arr[],const int arrSize,const string arrName){
    
    cout << arrName << "[]= ";
    for(int j=0; j<arrSize; j++){
        cout << arr[j] << " ";
    }
    cout << endl;
}

// function to perform reaction
void reaction_hop(Site lat[][LCOL], const int site_index, const string direction, 
                    mt19937 & mt)
{
    // local variables
    int nn_index, nnRowIndex, nnColIndex;
    int siRowIndex, siColIndex;
    int A_i;                    // rand chosen hopping species
    int A_j;                    // rand chosen hopped species
    vector<int> nnSpecies, onSiteSpecies;// how to know wat species are present

    // convert site_index of type int [0-SITES-1), into row and col of type ints
    siRowIndex = site_index / LCOL;
    siColIndex = site_index % LCOL;

    // check boundaries for 1d chain and get nearest neighbor index
    nn_index = get_nn(site_index, direction);

    // convert nn_index of type int [0-SITES-1), into row and col of type ints
    nnRowIndex = nn_index / LCOL;
    nnColIndex = nn_index % LCOL;

    // choose species at chosen Site w.r.t. its local pop. that will hop/react
    A_i = pick_reacting_species(lat[siRowIndex][siColIndex], mt);

    // determine the species types on current Site
    lat[siRowIndex][siColIndex].getAllSpeciesAtSite(onSiteSpecies);

    // determine which species are present on neighboring site
    lat[nnRowIndex][nnColIndex].getAllSpeciesAtSite(nnSpecies);
    
    // single-species annihilation model (carr cap & num species must equal 1)
    if(REACTION_TYPE == "SAM"){
        /*
            i will have to check, but the "SAM" reaction type may be 
            unnecessary because it may be a limiting case of "MAM" reaction type
        */
        // check if number of species is equal to 1
        if(NUM_SPECIES != 1){
            cout << "WARNING!! A single-species annihilation reaction cannot be"
                    << " performed if NUM_SPECIES= " << NUM_SPECIES << ". "
                    << "Double-check REACTION_TYPE initialization." << endl;
            exiting_via_exit1();
        }
        else if(K != 1){
            print_dashes();
            cout << "WARNING!! K =/= 1 and SAM reaction was chosen."<<endl;
            print_dashes();
            exiting_via_exit1();
        }

        // SAM, have single-species at current Site, hop to neighbor Site
        lat[nnRowIndex][nnColIndex].incrSingSpecPop(1, 0);
        lat[siRowIndex][siColIndex].decrSingSpecPop(1, 0);

        // check to see if Site pop has exceeded Site's max capacity 
        if(lat[nnRowIndex][nnColIndex].getTotSitePop() 
                > lat[nnRowIndex][nnColIndex].getCapacity()){
            lat[nnRowIndex][nnColIndex].setSingleSpeciesPopulation(0, 0);
        }
    }
    
    /*
    multispecies annihilation model
        - currently, multiple species can exit on the lattice at the same time 
            if random initialization place multiple species there and if K is
            large enough to account for them.
            - therefore, "MAM" reaction checks if an annihilation can 
                occur on the current Site (no hopping) first, then chosen 
                hopping particle can hop only if Site is occupied by one species
        - During hopping, annihilation only occurs between the hopping species
            and randomly generated different species on the nearest neighbor 
            site. A_i + A_j –> 0, where i =/= j
    */
    else if(REACTION_TYPE == "MAM"){
        // local else-if vars
        double randDouble = gen_rand_double(mt);

        // check if there are multiple species on same Site
        if(onSiteSpecies.size() > 1){
            // check if MAM occurs
            if(randDouble < MAM_PROB){
            
                // determine hopped species
                do{
                    // pick a random species until A_i =/= A_j
                    A_j = pick_reacting_species(lat[siRowIndex][siColIndex], mt);
                } while(A_i == A_j);

                // perform annihilation reaction on current Site
                lat[siRowIndex][siColIndex].decrSingSpecPop(1, A_i);
                lat[siRowIndex][siColIndex].decrSingSpecPop(1, A_j);
            }
        }
        
        // if multiple species do not exist on current Site, look at nn
        else{
            // no hop, nor annihil occurs if max Site capacity has been reached
            if(lat[nnRowIndex][nnColIndex].getTotSitePop() < K){
                // check if no neighbors are present
                if(nnSpecies.size() == 0){
                    // only hopping occurs, have A_i hop
                    lat[siRowIndex][siColIndex].decrSingSpecPop(1, A_i);
                    lat[nnRowIndex][nnColIndex].incrSingSpecPop(1, A_i);
                }
                
                // annihilation could occur
                else{
                    // choose nn species
                    A_j = pick_reacting_species(lat[nnRowIndex][nnColIndex],mt);

                    // if hopped species == hopping species, just hop
                    if(A_i == A_j){
                        lat[siRowIndex][siColIndex].decrSingSpecPop(1, A_i);
                        lat[nnRowIndex][nnColIndex].incrSingSpecPop(1, A_j);
                    }

                    // else, MAM could occur
                    else{
                        // attempt annihilation reaction
                        if(randDouble < MAM_PROB){
                            lat[siRowIndex][siColIndex].decrSingSpecPop(1, A_i);
                            lat[nnRowIndex][nnColIndex].decrSingSpecPop(1, A_j);
                        }
                        else{
                            // hopping occurs, but not MAM reaction
                            lat[nnRowIndex][nnColIndex].incrSingSpecPop(1, A_i);
                            lat[siRowIndex][siColIndex].decrSingSpecPop(1, A_i);
                        }

                    }
                }
                    
            } // end of nn Site pop < K if-statement
        
        } // end of nearest neighbor MAM reaction attempt
        
    } // end of "MAM" reaction type

    /*
        predator-prey reaction ––> predator = A (0), prey = B (1)
        spontaneous predator death: A –> 0 with prob mu
        prey birth: B –> B+B with prob sigma
        predation: A+B –> A+A with prob lambda
    */
    else if(REACTION_TYPE == "LV"){
        // local LV else-if statment variables
        double randVar = gen_rand_double(mt);
        double deathRandVar = gen_rand_double(mt);
        int attemptPredation = gen_rand_int(0, 1, mt);

        // first check if there's more than one species on current site
        if(onSiteSpecies.size() > 1){
            // attempt either predation/predator death reactions on current Site
            // make sure nnSite can accept a pred/prey after predation/birth
            if(A_i == 0){
                // attempt an on-Site predation reaction
                if(attemptPredation){
                    if(randVar < LAMBDA && 
                        lat[nnRowIndex][nnColIndex].getTotSitePop() 
                            < lat[nnRowIndex][nnColIndex].getCapacity())
                    {
                        // remove prey on-Site and place predator on nnSite
                        lat[siRowIndex][siColIndex].decrSingSpecPop(1, 1);
                        lat[nnRowIndex][nnColIndex].incrSingSpecPop(1, 0);
                    }
                }

                // or attempt a spontaneous predator death reaction
                else{
                    if(deathRandVar < MU){
                        // remove a predator from Site
                        lat[siRowIndex][siColIndex].decrSingSpecPop(1, 0);
                    }
                }
            }
            
            // attempt prey birth on neighboring Site
            else if(A_i == 1 && 
                lat[nnRowIndex][nnColIndex].getTotSitePop() 
                    < lat[nnRowIndex][nnColIndex].getCapacity() &&
                randVar < SIGMA)
            {
                lat[nnRowIndex][nnColIndex].incrSingSpecPop(1, 1);
            }
        }
        
        // attempt reactions on neighboring Site
        else{
            // determine the species composition of neighboring Site
            if(nnSpecies.size() == 0){
                // predator only attempts death for empty neighboring Site
                if(A_i == 0 && randVar < MU){
                    lat[siRowIndex][siColIndex].decrSingSpecPop(1, 0);
                }

                // prey attempts birth like normal
                else if(A_i == 1 && randVar < SIGMA){
                    lat[nnRowIndex][nnColIndex].incrSingSpecPop(1, 1);
                }

            }

            // neigboring Site is not empty
            else{
                // determine if predator was chosen for reaction attempt
                if(A_i == 0){
                    // attempt either a predation or pred. death reaction
                    if(attemptPredation){
                        // now look at neighbor and determine if prey is present
                        if(find(nnSpecies.begin(), nnSpecies.end(), 1) != 
                            nnSpecies.end())
                        {
                            // pred on Site and prey is on nnSite; predation
                            // nnSite cap doesnt have tobe checked bc nnSite pop is const
                            if(randVar < LAMBDA){
                                lat[nnRowIndex][nnColIndex].decrSingSpecPop(1, 1);
                                lat[nnRowIndex][nnColIndex].incrSingSpecPop(1, 0);
                            }
                        }
                    } 

                    // attempt predator death
                    else{
                        if(deathRandVar < MU){
                            lat[siRowIndex][siColIndex].decrSingSpecPop(1, 0);
                        }
                    }
                }

                // prey was chosen to attempt reactions
                else{
                    if(randVar < SIGMA && 
                        lat[nnRowIndex][nnColIndex].getTotSitePop() 
                            < lat[nnRowIndex][nnColIndex].getCapacity())
                    {
                        // increase prey population on neighboring Site
                        lat[nnRowIndex][nnColIndex].incrSingSpecPop(1, 1);
                    }
                }
            }
        
        } // end of nearest neighbor reaction attempt
        
    } // end of LV reaction type else-if statement

    /*
        may-leonard reaction ––> species Bi, i=1,2,3
        predation without birth: Bi + Bi+1 –> Bi with rate mu
        birth: Bi –> Bi + Bi with rate sigma
        coagulation: Bi + Bi –> Bi with rate lambda
    */
    else if(REACTION_TYPE == "ML"){
        // stuff
    }
    
    // make sure correct reaction type was initially defined
    else{
        cout << "WARNING!! Not a vaild REACTION_TYPE, choose another. "<<endl;
        exiting_via_exit1();
    }

    return;
}

// function to return with direction particle will hop and react
string reaction_direction(mt19937 & mt)
{
    // local vars
    string dir;
    int rand_int = gen_rand_int(0, 99, mt);

    // react to the left
    if(rand_int < 25){
        dir = "left";
    }
    // react to the right
    else if(rand_int >= 25 && rand_int < 50){
        dir = "right";
    }
    // react up
    else if(rand_int >= 50 && rand_int < 75){
        dir = "up";
    }
    // react down
    else{
        dir = "down";
    }

    return dir;
}

// function to return the index of the nearest neighbor
int get_nn(const int siIndex, const string dir)
{
    // local vars
    int nnIndex;

    // react to the left
    if(dir == "left"){
        // check edge case
        if(siIndex % LCOL == 0){
            nnIndex = siIndex + (LCOL - 1);
        }
        else{
            nnIndex = siIndex - 1;
        }
    }

    // react to the right
    else if(dir == "right"){
        // check edge case
        if(siIndex % LCOL == (LCOL - 1)){
            nnIndex = siIndex - (LCOL - 1);
        }
        else{
            nnIndex = siIndex + 1;
        }
    }

    // react up
    else if(dir == "up"){
        // check edge case
        if(siIndex < LCOL){
            nnIndex = siIndex + (SITES - LCOL); // same as adding by (LROW-1)*LCOL
        }
        else{
            nnIndex = siIndex - LCOL;
        }
    }

    // react down
    else{
        // check edge case
        if(siIndex >= SITES-LCOL){
            nnIndex = siIndex % LCOL;
        }
        else{
            nnIndex = siIndex + LCOL;
        }
    }

    return nnIndex;
}

// function to calculate the number density of each species
void get_num_density(const Site lat[][LCOL], const int desiredK, int speciesPop[])
{
    /*
    the passed argument, desiredK, is used for "half" lattice configurations and if 
    desiredK == -1, then if statement always executes, thus counting all species pops
    */

    // make sure the passed by reference int array is all zeros
    for(int i=0; i<NUM_SPECIES; i++){
        speciesPop[i] = 0;
    }

    // walk down lattice
    for(int i=0; i<LROW; i++){
        for(int j=0; j<LCOL; j++){
            // count each particle present for each species
            for(int k=0; k<NUM_SPECIES; k++){
                // make sure desired carrying-capacity is recorded
                if(desiredK == lat[i][j].getCapacity() || desiredK == -1){
                    speciesPop[k] += lat[i][j].getSingleSpeciesPopulation(k);
                }
            }
        }
    }

    return;
}

// add all elements within a 1d int array
int sum_1d_int_array(const int arr[], const int arrSize){
    // local vars
    int sum = 0;
    
    // traverse down array and add each element
    for(int i=0; i<arrSize; i++){
        sum += arr[i];
    }

    return sum;
}

int get_int_from_vec(const vector<int> & vec, int & index, mt19937 & mt){
    // local vars
    int chosenElement;

    // vector must not be empty!
    if(!vec.size()){
        cout << "\n\tWARNING!! vector passed to get_int_from_vec() was empty!"
                << endl;
        print_int_vector(vec, "vec");
        exiting_via_exit1();
    }

    // choose an index at random
    index = gen_rand_int(0, vec.size() - 1, mt);

    // get vector element from randomly chosen index
    chosenElement = vec.at(index);

    return chosenElement;
}

int pick_reacting_species(const Site s, mt19937 & mt){
    // local vars
    vector<int> availSpecies;
    int chosenSpecies, species;

    // loop through each species
    for(int i=0; i<NUM_SPECIES; i++){
        
        // insert species into vector depending on it's pop. on this Site
        for(int j=0; j<s.getSingleSpeciesPopulation(i); j++){
            availSpecies.push_back(i);
        }
    }

    // randomly choose int from vector
    chosenSpecies = get_int_from_vec(availSpecies, species, mt);

    return chosenSpecies;
}

// function to print to screen the initial conditions
void print_initial_conditions(const int currentSeed, const string runID){
    // make sure only two species for LV simulation with A=pred, B=prey
    if(NUM_SPECIES != 2 && REACTION_TYPE == "LV"){
        cout<<"WARNING!! LV reaction type was chosen, but NUM_SPECIES != 2"
                << endl;
        exiting_via_exit1();
    }
    
    cout << "seed= " << currentSeed << ", "
            << "ID= " << runID << ", "
            << NUM_SPECIES << " species, "
            << ( (CARRY_CAP_TYPE == "uniform")?(CARRY_CAP_TYPE+" K= "+to_string(K)):
                (CARRY_CAP_TYPE + ": K_LOW= " + to_string(K_LOW)
                + " & K_HIGH= " + to_string(K_HIGH)) )
            << ((LROW == LCOL) ? ", L= "+to_string(LROW)+", " : 
                (", LROW= " + to_string(LROW)+ " & LCOL= " + to_string(LCOL))) + ", "
            << NUM_STEPS << " mcs, "
            << INIT_TYPE << " intitalization, "
            << REACTION_TYPE << " reaction, "
            << "mu= " << MU << ", sigma= " << SIGMA << ", lambda= " << LAMBDA<< endl;
    
    print_1d_int_array(N_i, NUM_SPECIES, "N_i");
    return;
}

// function to exit via exit(1)
void exiting_via_exit1(){
    cout << "Exiting via exit(1)..." << endl;
    exit(1);

    return;
}

// funct to make sure the desired num of particles have been placed on the lat
void check_num_density(const Site lat[][LCOL]){
    // local vars
    int sums[NUM_SPECIES] = {0};
    int desired_density = sum_1d_int_array(N_i, NUM_SPECIES);
    int calculated_density;

    // loop through number of Sites
    // i think this nested loop is redoing what get_num_density() does
    for(int i=0; i<LROW; i++){
        for(int j=0; j<LCOL; j++)
            // loop through the number of species
            for(int k=0; k<NUM_SPECIES; k++){
                sums[k] += lat[i][j].getSingleSpeciesPopulation(k);
            }
    }
    calculated_density = sum_1d_int_array(sums, NUM_SPECIES);

    // exit simulation if the desired number of particles do not equal the calc
    if(desired_density != calculated_density){
        cout << "\n\tWARNING!! The number of particles placed on lattice Sites "
                <<"do not equal the number of total desired particles!" << endl;
        print_1d_int_array(N_i, NUM_SPECIES, "N_i");
        print_1d_int_array(sums, NUM_SPECIES, "sums");
        exiting_via_exit1();
    }

    return;
}

// funct to print to screen the current progress of the simulation
void print_simulation_progress(const int timeStep, const int popArr[]){
    
    // print simulation progress to screen
    cout << "finished " << timeStep << "/" << NUM_STEPS << " steps ";
    print_sys_time();
    print_1d_int_array(popArr, NUM_SPECIES, "num_density");
    cout << endl;

    return;
}

void print_sys_time(){
    // local vars
    time_t currentTime;

    // get current time
    time(& currentTime);

    // print update to screen
    cout <<  "@ " << ctime(& currentTime);

    return;
}

// function to make sure command line agrument is submitted
void check_command_line_inputs(char *arr[], const int arr_size){

    // make sure one command line arguments were given
    if(arr_size != 2){
        cout << "WARNING!! Include runID as single command line argument" << endl;
        
        cout << "argv[]= [";
        for(int i=0; i<arr_size; i++){
            cout << arr[i] << (i == arr_size - 1 ? "]" : ", ");
        }
        cout << endl;
        
        exiting_via_exit1();
    }

    return;
}

int get_total_lattice_pop(const Site lat[][LCOL]){
    // local vars
    int totPop=0;
    
    for(int i=0; i<LROW; i++){
        for(int j=0; j<LCOL; j++){
            totPop += lat[i][j].getTotSitePop();
        }
    }

    return totPop;
}

// funct to cout the tot pop of each species w/n the lattice
void print_lattice_pop(const Site lat[][LCOL]){
    // local vars
    int popSums[NUM_SPECIES] = {0};

    // walk through lattice
    // i think this nested loop is redoing what get_num_density() does
    for(int i=0; i<LROW; i++){
        for(int j=0; j<LCOL; j++){
            // loop through the number of species
            for(int k=0; k<NUM_SPECIES; k++){
                // sum up species pop
                popSums[k] += lat[i][j].getSingleSpeciesPopulation(k);
            }
        }
    }

    // print the total population of each species to screen
    cout << "total population of each species= ";
    for(int i=0; i<NUM_SPECIES; i++){
        cout << popSums[i] << " ";
    }
    cout << endl;

    return;
}

// function to print current lattice configuration to screen
void output_lattice_config(ofstream & out, Site lat[][LCOL], const int mcs){
    // output current run and current Monte Carlo timestep
    out << "# mcs= " << mcs << endl;

    // walk down lattice and print to screen
    for(int i=0; i<LROW; i++){
        for(int j=0; j<LCOL; j++){
            out << lat[i][j] << " ";

        }
        out << endl;
    }

    return;
}

int find_Site(const Site lat[][LCOL],const vector<int> & vec,const int randParticle){
    // local vars
    int Site_index = -1;
    int row, col;
    int sum = 0;

    // traverse through the available Sites
    for(int i : vec){
        // convert Site index to coordinates (row, col)
        row = i / LCOL;
        col = i % LCOL;

        // add up particles until desired particle is reached
        sum += lat[row][col].getTotSitePop();

        // check if desired particle is reached
        if(sum >= randParticle){
            Site_index = i;
            break;
        }
    }

    return Site_index;
}

void output_num_density(ofstream & out, Site lat[][LCOL], int pops[], const int mcs){
    
    // output time step to file
    out << mcs << "\t";
    
    // output scheme for spatially uniform carrying-capcity
    if(CARRY_CAP_TYPE == "uniform"){
        // get populations of each species
        get_num_density(lat, K, pops);

        // output to file
        for(int i=0; i<NUM_SPECIES; i++){
            out << pops[i] << "\t";
        }
    }

    // output scheme for spatially heterogeneous carrying-capcity (2 regions)
    else if(CARRY_CAP_TYPE == "half"){
        // get populations of each speices for each region
        for(int i : {K_LOW, K_HIGH}){
            // get populations for each species depending on desired K
            get_num_density(lat, i, pops);

            // output to file
            for(int j=0; j<NUM_SPECIES; j++){
                out << pops[j] << "\t";
            }
        }

        /* the array, pops[], must return with the the total population of each 
        species regardless of CARRY_CAP_TYPE because the number of attempted 
        reactions depend of the total number of particles w/n the system. */
        get_num_density(lat, -1, pops);
    }

    // any other CARRY_CAP_TYPE is not supported at the moment
    else{
        cout << "\n\tWARNING!! Incorrect CARRY_CAP_TYPE! Cannot output species "
                << "populations correctly." << endl;
        exiting_via_exit1();
    }

    // row finished, begin new line
    out << endl;

    return;
}