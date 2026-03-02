/*
kenneth distefano, summer 2023
FILE: Site.h
PURPOSE: v11

UPDATES FROM PREVIOUS VERSION (v8):
    - 
*/

#ifndef SITE_H
#define SITE_H
#include <vector>

using namespace std;

/* -------------------------- global constants ---------------------- */
const int A_ASCII_INT_VALUE = 65;
const int NUM_SPECIES = 2;
const int DEFAULT_SPECIES_ARR[NUM_SPECIES] = {0};
const int DEFAULT_CAPACITY = NUM_SPECIES;

/* -------------------------- type declarations ---------------------- */
class Site
{
    public:
        // desc: Site() is a constructor which passes all member variables 
        Site(const int species[], const int capacity);

        // desc: Site() is the default constructor 
        Site();

        // desc: setEntireSpeciesArr() accesses the member array of the object
        //          to give each index a value
        void setAllSpeciesPopulation(const int species[]);

        // set the population of a specific species
        void setSingleSpeciesPopulation(const int pop, const int speciesLoc);

        // desc: setCarryCapacity() accesses the member variable of the object 
        //          to give it a value
        void setCapacity(const int cap);

        // desc: getCapacity() allows for access of private member variable 
        //          m_capacity
        int getCapacity() const {return m_capacity;}

        // get the current population of a specific species
        int getSingleSpeciesPopulation(const int speciesIndex) const;

        // get species type (convert int into char)
        char getSingleSpeciesType(const int speciesIndex) const;

        // desc: operator overload
        // post: now << can be used to output to screen the type Site. Print to
        //          screen information about one site
        // note: output becomes weird if number of speices exceeds 26 because
        //          the ASCII values for A, B, C, ..., are used
        friend ostream & operator << (ostream & os, Site & site);

        // increase population of specific species
        void incrSingSpecPop(const int dPop, const int speciesLoc);

        // decrease population of specific species
        void decrSingSpecPop(const int dPop, const int speciesLoc);

        // get the total population of the Site (sum of all species population)
        int getTotSitePop() const;

        // desc: determine how many speices are present at Site
        // post: returns single int
        int getNumSpeciesAtSite() const;

        // desc: determine which specific species are present at Site
        void getAllSpeciesAtSite(vector<int> & vec);

    private:
        int m_speciesArr[NUM_SPECIES];
        int m_capacity;
};

#endif
