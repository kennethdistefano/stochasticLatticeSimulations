// kenneth distefano, summer 2023
// FILE: Site.cpp
// PURPOSE: v9 Site.cpp contains the function definitions for the Site class
/* LAST UPDATE: on Feb. 21, 2024 
    --> overloaded operator used cout, not os!! within output_lattice_config() 
    within functions_rk_v*.cpp, function definition for operator <<() was just 
    copied and pasted...
*/

#include <iostream>
#include "Site.h"
#include "vector"

using namespace std;

Site :: Site(const int species[], const int capacity)
{
    setAllSpeciesPopulation(species);
    setCapacity(capacity);
    
}

Site :: Site()
{
    // set default values
    for(int i=0; i<NUM_SPECIES; i++){
        m_speciesArr[i] = DEFAULT_SPECIES_ARR[i];
    }
    m_capacity = DEFAULT_CAPACITY;
}

void Site :: setAllSpeciesPopulation(const int species[])
{
    //traverse down passed array to equate corresponding indices to member array
    for(int i=0; i<NUM_SPECIES; i++){
        m_speciesArr[i] = species[i];
    }
    return;
}

void Site :: setSingleSpeciesPopulation(const int pop, const int speciesLoc){
    m_speciesArr[speciesLoc] = pop;
    return;
}

void Site :: setCapacity(const int cap)
{
    m_capacity = cap;
    return;
}

int Site :: getSingleSpeciesPopulation(const int speciesIndex) const{
    return m_speciesArr[speciesIndex];
}

char Site :: getSingleSpeciesType(const int speciesIndex) const{
    return char(speciesIndex + A_ASCII_INT_VALUE); 
}

ostream & operator << (ostream & os, Site & site){

    os << "(";
    for(int i=0; i<NUM_SPECIES; i++){
        os << site.getSingleSpeciesPopulation(i) 
                << (i == NUM_SPECIES - 1 ? ";" : ",");
    }
    os << site.getCapacity() << ")";
    
    return os;
}

void Site :: incrSingSpecPop(const int dPop, const int speciesLoc){
    m_speciesArr[speciesLoc] += dPop;
    return;
}

void Site :: decrSingSpecPop(const int dPop, const int speciesLoc){
    m_speciesArr[speciesLoc] -= dPop;
    return;
}

int Site :: getTotSitePop() const{
    // local vars
    int popSum = 0;

    // loop through member species array
    for(int i=0; i<NUM_SPECIES; i++){
        popSum += m_speciesArr[i];
    }
    return popSum;
}

int Site :: getNumSpeciesAtSite() const{
    // local vars
    int counter = 0;
    
    // traverse down member species array
    for(int i=0; i<NUM_SPECIES; i++){
        // if species population is nonzero, then true
        if(m_speciesArr[i]){
            counter++;
        }
    }

    return counter;

}

void Site :: getAllSpeciesAtSite(vector<int> & vec){
    // walk down species array
    for(int i=0; i<NUM_SPECIES; i++){
        // check if specific species is present
        if(m_speciesArr[i]){
            vec.push_back(i);
        }
    }

    return;
}
