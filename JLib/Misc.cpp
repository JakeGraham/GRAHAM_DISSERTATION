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
#include<Misc.h>

// reads txt file into vector<double> 
std::vector<double> ReadTxt(std::string infilename){
        std::ifstream infile;
        infile.open(infilename);
        std::vector<double> v1;
        std::string temp;
        // Imports data by line into vector<double>
        while (std::getline(infile, temp)) {
                std::istringstream buffer(temp);
                std::vector<double> line((std::istream_iterator<double>(buffer)),std::istream_iterator<double>());
                for(size_t i = 0; i < line.size(); i++){
                        v1.push_back(line[i]);
                }
        }
        infile.close();
        return(v1);
}


/* takes a string and a delimeter and splits the string by the delimeter and returns a vector of string for each substring afte the split
	inputs:

	outputs:
*/
std::vector<std::string> SplitString(std::string SIN, std::string delim){
        std::vector<std::string> parsed;
	std::string token;
	size_t pos = 0;
        while((pos = SIN.find(delim)) != std::string::npos){
                token = SIN.substr(0,pos);
                parsed.push_back(token);
                SIN.erase(0, pos + delim.length());
        }
        parsed.push_back(SIN);
	return(parsed);
}



// given two points XY returns distance between the points... inputs = doubles
double DP2P(double x1, double y1, double x2, double y2){
        double Dist = pow( ( pow( (x2 - x1), 2) + pow( (y2 - y1), 2) ), .5 );
        return(Dist);
}

// given two points XY returns distance between the points... inputs = 32 bit signed ints
double DP2PSI(int32_t x1, int32_t y1, int32_t x2, int32_t y2){
        int32_t Dist = pow( ( pow( (x2 - x1), 2) + pow( (y2 - y1), 2) ), .5 );
        return(Dist);
}


// returns intersection point (XY) of two lines given their XY coordinates.. first four inputs are for line 1, last four are for 2
std::vector<double> IntrsctPt(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
        std::vector<double> IP(2);
        IP[0] = ( (x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4) ) / ( (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4) );
        IP[1] = ( (x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4) ) / ( (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4) );
        return(IP);
}


// takes pointers to x components of two points in 2 dimensions, the next memory location should be the y components respectively..needed for D2L.. distance between two points squared
double dist2(double *p1, double *p2){
        double d2;
        d2 = pow(fabs(*p1 - *p2), 2) + pow(*(p1+1) - *(p2+1), 2);
        return(d2);
}


/* Measures the distance between a point (2D) and a line (2D)
        inputs: LS = a pointer to the x component of the starting point of the line, the next memory location should be the y component
                LE = same as LS but for the end of the line
                P =  a pointer to the x component of the point to be measured to the line, next memory location should be the y component
        output: the distance between point P and the line between LS and LE
*/ 
double D2L(double *LS, double *LE, double *P){
        double l2 = dist2(LS, LE);
        if(l2 == 0){
                return(dist2(P, LS));
        } else {
                double t = ((*P - *LS) * (*LE - *LS) + (*(P+1) - *(LS+1)) * (*(LE+1) - *(LS+1))) / l2;
                double tmp[2];
                tmp[0] = *LS + t * (*LE - *LS); 
                tmp[1] = *(LS+1) + t * (*(LE+1) - *(LS+1));
                return(pow(dist2(P, &tmp[0]), .5));
        }
}


/* function that gives the 1D grid position of the value
        input:  gridstart = location the grid starts from
                gridcell = the size of the grid cell
                value = the value being evaluated 
        output: gridpos = the grid cell "value" should be placed in

		!!!!!!!!!!  BUG! if the (value - gridstart) term is a larger negative than the data structure can hold, converts to a positive

*/
int64_t gridpos(int64_t gridstart, int64_t gridcell, int64_t value) {
        int64_t gridpos = floor((value - gridstart) / gridcell);
        return (gridpos);
} 


/* function that gives the 1D voxel given xyz gridcell
        inputs: x = x grid cell (from gridpos)
                y = y grid cell (from gridpos)
                z = z grid cell (from gridpos)
                xc = # of cells in the x grid
                yc = # of cells in the y gird
        output: v = the voxel position in 1D... 1D projection hierachry = z > y > x 
*/
uint32_t vxl(uint32_t x,uint32_t y,uint32_t z,uint32_t xc,uint32_t yc){
        uint32_t v = (z*(xc*yc))+(y*(xc))+x;
        return(v);
}


/* function that gives the 1D position for a reporojected 2D grid
        inputs: xpos = the x cell from gridpos 
                ypos = the y cell from gridpos
                xc = # of cells in the x grid
        output: gp = the 2D grid position projected in 1D
*/
uint32_t GridPos2D(int32_t xpos, int32_t ypos, uint32_t xc){
        uint32_t gp = (ypos*xc) + xpos;
        return(gp);
}



/* takes a 2D vector projected into 1D and outputs a 2D matrix as a text file
	inputs:	VIN = a pointer to the vector to be output
		DL = Delineation length, the number of cells in the x grid
		OFName = the name for the output file
	outputs: n/a
*/
void WriteVecToMatrix(std::vector<double> *VIN, uint32_t DL, char OFName[150]){
        std::ofstream ofile;
        ofile.open(OFName); 
        for(uint32_t i = 0; i < VIN->size() / DL; i++){
                for(uint32_t j = 0; j < DL; j++){ 
                        ofile << std::setprecision(12) << VIN->at((DL * i) + j) << " ";
                }
        ofile << std::endl;
        }
}

void WriteBVecToMatrix(std::vector<bool> *VIN, uint32_t DL, char OFName[150]){
	std::ofstream ofile;
        ofile.open(OFName);
	for(uint32_t i = 0; i < VIN->size() / DL; i++){
		for(uint32_t j = 0; j < DL; j++){
			if(VIN->at((DL * i) + j)){
				ofile << 1 << " ";
			} else {
				ofile << 0 << " ";
			}
		}
	ofile << std::endl;
	}	
}


/* Takes a vector and returns the sorted indices for the input vector.. for example [3 1 8 5] as input would retur
	a vector of [2 1 4 3]
	input: v = the vector to be sorted
	output: idx = a vector of sorted indices
*/
template <typename T> 
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

        // initialize original index locations
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        std::sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

        return idx;
}


/* Takes vector and fills in missing values with a specified value
       inputs:
                VIN = pointer to vector to be filled
                Eval = the value that is considered "empty" (usually 0)
                impute = the value that "empty" cells will be filled with
       outputs: NA
*/
template <typename T>
void GridImpute(std::vector<T> *VIN, T Eval, T impute){
     for(size_t i = 0; i < VIN->size(); i++){
          if(VIN->at(i) == Eval){
               VIN->at(i) = impute;
          }
     }
}


/*  adds (or subtracts) a value from all elements in a vector
     inputs:
          VIN = a pointer to the grid (vector) to be translated
          subtract = bool determining whether to add or subtract Trans
          Trans = value to be added (or subtracted)
     outputs: NA
*/
template <typename T>
void GridTranslate(std::vector<T> *VIN, bool subtract, T Trans){
     if(subtract){
          for(size_t i = 0; i < VIN->size(); i++){
               VIN->at(i) = VIN->at(i) - Trans;
          }
     } else {
          for(size_t i = 0; i < VIN->size(); i++){
               VIN->at(i) = VIN->at(i) + Trans;
          }
     }
}





/* Returns the value for the specified quantile of a vector.. Generic to work with any type
     intputs:
          VIN = pointer to the vector the quantile is being calculated for
          Quant = quantile to be calculated (e.g., if Quant = 0.05 the 5th quantile of VIN will be returned
     outputs:
          QV = the value of the quantile extracted from VIN
*/
template <typename T>
T VecQuantile(std::vector<T> *VIN, double Quant){
     std::sort(VIN->begin(), VIN->end());
     size_t QE = round(VIN->size() * Quant);
     T QV = VIN->at(QE);
     return(QV);
}



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
std::vector<double> LogisticWeights(std::vector<double> *VIN, double x0, double L, double k, double trans, bool invert){
     std::vector<double> OP(VIN->size());
     double tmp;
     if(invert){
          for(size_t i = 0; i < OP.size(); i++){
               OP[i] = (( L / (1 + exp(-k*(VIN->at(i)-x0))) ) * -1) + trans;
          }
     } else {
          for(size_t i = 0; i < OP.size(); i++){
               OP[i] = ( L / (1 + exp(-k*(VIN->at(i)-x0))) ) + trans;
          }
     }
     return(OP);
}


/* calculates the mean of a vector
      inputs =
           VIN a pointer to the vector to calculate the mean from 
      
      outputs:
           OP = the mean of VIN ... NOTE!! this is templated so the mean of integers will return integers not floats
*/
template <typename T> 
T VecMean(std::vector<T> *VIN){
     T OP = std::accumulate(VIN->begin(), VIN->end(), 0.0) / VIN->size();
     return(OP);
} 


/*
  calculates population standard deviation
     inputs: 
          VIN = a pointer to the vector of type T to calculate the standard deviation from
 
     outputs: 
          OP = standard deviation of VIN, returns a double regardless of type T
*/
template <typename T>
double VecSD(std::vector<T> *VIN){
     std::vector<double> tmp(VIN->size());
     double vmean = VecMean(VIN);
     for(size_t i = 0; i < VIN->size(); i++){
          tmp[i] = pow((vmean - VIN->at(i)), 2);
     }
     double OP = pow(std::accumulate(tmp.begin(), tmp.end(), 0.0) / tmp.size(), .5);
     return(OP);
}


/*  makes a binary classification based on a threshold for the input vector, can specify whether true is above or below the threshold
     inputs:
          VIN = a pointer to the vector the classification is being made on ( type T)
          T = threshold the classification is determined by ( type T)
          below = bool specifying whether above or below the threshold is TRUE

     outputs:
          OP = bool vector with binary classifictions   

*/
template <typename T> 
std::vector<bool> GridClassThresholdBool(std::vector<T> * VIN, T thresh, bool below){
     bool yaynay;
     if(below){
          yaynay = true;
     } else {
          yaynay = false;
     }

     std::vector<bool> OP(VIN->size());
     for(size_t i = 0; i < OP.size(); i++){
          if(VIN->at(i) > thresh){ 
               OP[i] = !yaynay;
          } else {
               OP[i] = yaynay;
          }
     }
     return(OP);
}



/* converts a vector<bool> to vector<double> .. useful for outputing binary vector (see GridClassThresholdBool) to txt file
     inputs:
          VIN = a pointer to the binary vector (classification vector)

     outputs:
          OP = vector<double> of either 1 or 0 specifying TRUE and FALSE, respectively
*/
std::vector<double> BoolVtoDblV(std::vector<bool> * VIN){
     std::vector<double> OP(VIN->size());
     for(size_t i = 0; i < OP.size(); i++){
          if(VIN->at(i)){
               OP[i] = 1;
          } else {
               OP[i] = 0; 
          }
     }
     return(OP);
}



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

int32_t PointInPolygon(int32_t nvert, double *vertx, double *verty, double testx, double testy){
     int32_t i, j, OP = 0;
     for (i = 0, j = nvert-1; i < nvert; j = i++) {
          if ( ((verty[i]>testy) != (verty[j]>testy)) && (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) ){
               OP = !OP;
          }
     }
     return(OP);
}


/* Converts a string to a char array
     inputs:
          SIN = a pointer to the string to be converted
     outputs:
          OP = a pointer to the char array

*/ 
char * StringToCharArray(std::string *Sin){
     char * OP;
     OP = (char*) malloc((Sin->length() + 1) * sizeof(char));
     strcpy(OP, Sin->c_str());
     return(OP);
}



/* converts a double to a string, rounded at the specified level
      inputs:
         value = value to be rounded and returned
         decimals = the number of decimals to round to 
      outputs:
         s = string with "value" rounded to "decimals"
 */

std::string toStrMaxDecimals(double value, int decimals){
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(decimals) << value;
    std::string s = ss.str();
    if(decimals > 0 && s[s.find_last_not_of('0')] == '.') {
        s.erase(s.size() - decimals + 1);
    }
    return s;
}


/* Writes a 2D vector to a text file 




*/
template <typename T> 
void Write2DVecToTxt(std::vector<std::vector<T> > VIN, std::string FName){
   std::ofstream Ofile;
   Ofile.open(FName);    
   for(size_t i = 0; i < VIN[0].size() ; i++){
      for(size_t j = 0; j < VIN.size(); j++){
         Ofile << std::setprecision(10) << VIN[j][i] << " ";
      }
      Ofile << std::endl;
   }
   Ofile.close();
}



template <typename T> 
std::vector<T> VecRemZeros(std::vector<T> *VIN){
   size_t ctr = 0;
   std::vector<T> OP(VIN->size());
   for(size_t i = 0; i < VIN->size(); i++){
      if(VIN->at(i) != 0){
         OP[ctr] = VIN->at(i);
         ctr++; 
      }
   }
   ctr = ctr - 1;
   OP.resize(ctr);
   return(OP);

}













