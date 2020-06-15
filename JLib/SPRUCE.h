#ifndef SPRUCE_h
#define SPRUCE_h

#include<iostream>
#include<string>
#include<fstream>
#include<stdint.h>
#include<vector>
#include<LasReadWrite.h>
#include<Misc.h>
#include<PCProc.h>
#include<iomanip>



/* takes a sting of a SPRUCE .las file and parses by "_" and "." for using year, plot etc...
	inputs: SIN = a string of the .las file name
	outputs: OP = vector of strings that are separated by "_" and "."
*/
std::vector<std::string> ParseSPRUCELas(std::string SIN);


// Takes boardwalk UTMs and returns the mean distance between opposing sides (2D)
double BWML(std::vector<double> BWC);


// Takes boardwalk UTMs and and returns UTMs of the 8 vertices
std::vector<double> BWVerts(std::vector<double> BWC);


/* Finds the shortest distance from a point to a SPRUCE boardwalk
        inputs: BWV = a pointer to a vector of boardwalk vectices.. from BWVerts 
*/
double BWD (double *BWV, double *p);


/* Calculates the distance to the boardwalk for all points in a LasFile
        inputs: LF = pointer to LasFile to be operated on
                BMV = pointer to vector or array of boardwalk vertices in the order x1,y1,x2,y2,.. etc
        output: OPV = A vector of boardwalk distance for each point in LF in order

*/
std::vector<double> LasBWD(LasFile *LF, double *BWV);


/* structure for SphagClosesRefPts to output
*/
struct SpClRefPtsO{
   LasFile LasF;
   std::vector<double> Error;
   std::vector<double> RefEle;
   std::vector<double> ReConEle;
};


// Legacy..
struct GenO{
        LasFile LasF;
        std::vector<double> OPV;
};


/* Finds the closest(XY) sphagnum points to reference points inputs should be scale converted first
        inputs: *SPG = a pointer to the LasFile containing sphagnum surface points
                *RP = a pointer to the LasFile containing reference points
        output: A LasFile containing only the points from SPG that were closest to each reference point (XY distance).. copies SPG header
        
*/
SpClRefPtsO SphagClosestRefPts(LasFile *SPG, LasFile *RP);




/* Performs the accuracy assessment for sphagnum points for spring 2017.. Takes a surface file, full plot point cloud, and reference file

*/
void SphagAcc(LasFile * SPG, LasFile * RF, LasFile * PC, std::string ACCFile, std::string BWFile, int32_t plot, int32_t year, double);


/* iteratively calculates the percent of each chamber classified as hollow as a function of tolerance and mean water table level.. pretty sure this is redundant and better function below...
*/
int PctHolWT();


// calculates all elevation profiles for 2017 sphagnum files and returns as 2D vector of doubles (UTM)
std::vector<std::vector<double> > ElevBatch();


// calculates elevation quantiles for all spruce plots 
void EleQuantBatch();


// iteratively calculates the percent of each chamber classified as hollow as a function of tolerance and quantile
void EleThresh();


// calculates percent hollow based on threshold and quantile... iteratively
void HolFromWT();

// calulates percent of plot classified as hollow based on a 40 percent of max hummock height.. The original definition used at SPRUCE.. min + (Max - Min)*.4 
std::vector<std::vector<double> > PctHolPctMaxEle(std::vector<double> *Zin, bool sort);

// structure output by SPRUCE_Hol containing all steps of processing and parameters used
struct SPRUCE_HolO{
     int32_t DL; // Deliniation length 
     int32_t Gsig; // sigma used in gaussian smoothing for gradient
     int32_t Lsig; // sigma used in gaussian smoothing for laplacian
     double Gthresh; // gradient threshold used
     double Lthresh; // laplacian threshold used
     double GridSize; // grid size 
     std::vector<double> Raw; // raw elevation data
     std::vector<double> GGB; // smoothed elevation used for gradient
     std::vector<double> LGB; // smoothed elevation used for laplacian
     std::vector<double> Grad; // gradient after smoothing
     std::vector<double> Lap; // laplacian after smoothing
     std::vector<double> GradT; // gradient threshold, displays elevation
     std::vector<double> LapT; // laplacian threshold, displays elevation
     std::vector<double> HumHol; // classification of hollow, displays hollows raw elevation
};





/* Performs classification of hollows from microtopography point cloud using first and second derivatives
     inputs:
          Lasf = pointer to LasFile to be classified 
          Gsig = sigma for gaussian smoothing for gradient
          Lsig = sigma for gaussian smoothing for laplacian
          Gthresh = gradient threshold
          Lthresh = laplacian threshold
          GS = grid cell size
     outputs:
          SPRUCE_HolO = Object containing vectors of the orignial data and every processing step until final classification
*/
SPRUCE_HolO SPRUCE_Hol(LasFile *Lasf, int32_t Gsig, int32_t Lsig, double Gthresh, double Lthresh, double GS);


#endif
