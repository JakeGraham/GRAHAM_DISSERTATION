#ifndef Misc_h
#define Misc_h

#include<vector>
#include<iostream>
#include<chrono>
#include<stdint.h>
#include<sstream>
#include<istream>
#include<cmath>
#include<cstring>
#include<algorithm>
#include<iterator>
#include<fstream>
#include<iomanip>


// reads txt file into vector<double> 
std::vector<double> ReadTxt(std::string infilename);


/* takes a string and a delimeter and splits the string by the delimeter and returns a vector of string for each substring afte the split
	inputs:

	outputs:
*/
std::vector<std::string> SplitString(std::string SIN, std::string delim);


// given two points XY returns distance between the points... inputs = doubles
double DP2P(double x1, double y1, double x2, double y2);


// given two points XY returns distance between the points... inputs = 32 bit signed ints
double DP2PSI(int32_t x1, int32_t y1, int32_t x2, int32_t y2);


// returns intersection point (XY) of two lines given their XY coordinates.. first four inputs are for line 1, last four are for 2
std::vector<double> IntrsctPt(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);


// takes pointers to x components of two points in 2 dimensions, the next memory location should be the y components respectively..needed for D2L.. distance between two points squared
double dist2(double *p1, double *p2);


/* Measures the distance between a point (2D) and a line (2D)
        inputs: LS = a pointer to the x component of the starting point of the line, the next memory location should be the y component
                LE = same as LS but for the end of the line
                P =  a pointer to the x component of the point to be measured to the line, next memory location should be the y component
        output: the distance between point P and the line between LS and LE
*/ 
double D2L(double *LS, double *LE, double *P);


/* function that gives the 1D grid position of the value
        input:  gridstart = location the grid starts from
                gridcell = the size of the grid cell
                value = the value being evaluated 
        output: gridpos = the grid cell "value" should be placed in

		!!!!!!!!!!  BUG! if the (value - gridstart) term is a larger negative than the data structure can hold, converts to a positive

*/
int64_t gridpos(int64_t gridstart, int64_t gridcell, int64_t value);


/* function that gives the 1D voxel given xyz gridcell
        inputs: x = x grid cell (from gridpos)
                y = y grid cell (from gridpos)
                z = z grid cell (from gridpos)
                xc = # of cells in the x grid
                yc = # of cells in the y gird
        output: v = the voxel position in 1D... 1D projection hierachry = z > y > x 
*/
uint32_t vxl(uint32_t x,uint32_t y,uint32_t z,uint32_t xc,uint32_t yc);


/* function that gives the 1D position for a reporojected 2D grid
        inputs: xpos = the x cell from gridpos 
                ypos = the y cell from gridpos
                xc = # of cells in the x grid
        output: gp = the 2D grid position projected in 1D
*/
uint32_t GridPos2D(int32_t xpos, int32_t ypos, uint32_t xc);


/* takes a 2D vector projected into 1D and outputs a 2D matrix as a text file
	inputs:	VIN = a pointer to the vector to be output
		DL = Delineation length, the number of cells in the x grid
		OFName = the name for the output file
	outputs: n/a
*/
void WriteVecToMatrix(std::vector<double> *VIN, uint32_t DL, char OFName[150]);

// same as WriteVecToMatrix() but handles vector<bool> instead of vector<double>
void WriteBVecToMatrix(std::vector<bool> *VIN, uint32_t DL, char OFName[150]);


/* Takes a vector and returns the sorted indices for the input vector.. for example [3 1 8 5] as input would retur
	a vector of [2 1 4 3]
	input: v = the vector to be sorted
	output: idx = a vector of sorted indices
*/
template <typename T> 
std::vector<size_t> sort_indexes(const std::vector<T> &v);

/* Checks if a given point is inside or outside a polygon give its vertices
     inputs:
          nvert = number of polygon vertices
          vertx = a pointer to an array of x coordinates of polygon vertices
          verty = a pointer to an array on y coordinates of polygon vertices
          testx = x coordinate of point to be evaluated
          testy = y coordiante of point to be evalutted
 
     outputs:
          OP = indicates whether the test point is inside (1) or outside (0) the polygon
*/

int32_t PointInPolygon(int32_t nvert, double *vertx, double *verty, double testx, double testy);

/* calculates the mean of a vector
      inputs =
           VIN a pointer to the vector to calculate the mean from 
      
      outputs:
           OP = the mean of VIN ... NOTE!! this is templated so the mean of integers will return integers not floats
*/
template <typename T>
T VecMean(std::vector<T> *VIN);


/*
  calculates population standard deviation
     inputs: 
          VIN = a pointer to the vector of type T to calculate the standard deviation from
 
     outputs: 
          OP = standard deviation of VIN, returns a double regardless of type T
*/
template <typename T>
double VecSD(std::vector<T> *VIN);

/* Creates a logistic curve with paramaters x0, L, k, and trans, and applies that relationship to the input vector... From the equation y = L / (1 + exp(-k(xi - x0))), trans is an optional paramater used to translate output values (adds this value to all output values)
     inputs: 
          VIN = pointer to the vector weights are to be generated for
          X0 = the centerpoint of the logistic curve (horizontal midpoint)
          L = the top vertical asymptote
          k = the "steepness" of the curve
          trans = value that is added to the output vector, useful for modifying possible weighted values (vertical midpoint)
     outputs:
          OP = a vector containing the weights from the logistic curve at x locations corresponding to VIN values
*/
std::vector<double> LogisticWeights(std::vector<double> *VIN, double x0, double L, double k, double trans, bool invert);


#endif                            
