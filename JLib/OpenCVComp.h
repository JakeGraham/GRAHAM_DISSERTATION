#ifndef OpenCVComp_h
#define OpenCVComp_h

#include<iostream>
#include<iomanip>
#include<stdint.h>
#include<Misc.h>
#include<PCProc.h>
#include<LasReadWrite.h>
#include<SPRUCE.h>
#include<stdio.h>
#include<opencv2/imgproc.hpp>
#include<opencv2/imgcodecs.hpp>
#include<opencv2/highgui.hpp>
using namespace cv;


/* Converts a 2d grid of doubles (projected in 1D) to a openCV Mat object (done through memcpy). Also has the option to translate the data (i.e., subtract a value) which helps with displaying in openCV.
     inputs:  
              Grid = a pointer to a 2D grid (projected in 1D) to be converted to cv::Mat object 
              DL = delineation length (# of x grid cells)
              translate = bool that determines if Grid values will be shifted by subtracting transV from each grid value 
              transV = value to be used in translation
     outputs: 
              OP = the cv::Mat object that now contains the "pixel" values from Grid

                   !!! NOTE!! if translate = FALSE transv is inconsequential !!!
*/
Mat GridToCVMat(std::vector<double> *Grid, int32_t DL, bool translate , double transV);  


/* Structure to be output by CVMatToGrid. 
     Grid = the vector representing the 2D grid
     DL = the delineation length (i.e., # of x grid cells)
*/
struct MTG{
     std::vector<double> Grid;
     int32_t DL;
};


/* Converts a cv::Mat object to a 2D grid projected in 1D. Currently only works with a single band 64 bit cv::Mat image, will likely modify to accept multiband images and images with data types other than CV_64F (double).
     input: 
               MT = a pointer to the cv::Mat object to be converted
     outputs:
               OP = MTG object containing the vector (grid) and DL (delineation length i.e., the number of x grid cells)
*/
MTG CVMatToGrid(Mat *MT);


/* Creates a new Mat object that only contains values from the "IM" Mat where the "Mask" Mat is true (i.e., != 0)
     inputs: 
              IM =  pointer to Mat object from which values will be extracted for output
              Mask =  pointer to Mat object that will be used to "mask" IM
     outputs: 
              OP = Mat object containing values from IM where Mask != 0
              OPF = empty Mat object that is returned if dimensions of IM and Mask don't match

 !!!!!!!!!!!!!!!! NOTE ~~~~~~~~~~ Only works with Mat objects storing doubles.. may adjust to be generic.. not sure how ATM
*/
Mat MaskMat(Mat *IM, Mat *Mask);


/* Calculates the image gradient using sobel operators in OpenCV to a std::vector<double>
     inputs:
              M = pointer to Mat object to calculate the gradient of 
              Mat = if Mat == true, the output is an object of type CV::Mat, if false output is std::vector<double>
     outputs:
              OP = 2d grid (projected in 1D) containing absolute gradient from M
*/
std::vector<double> CVImGradGrid(Mat *M);


/*  Scales the values of one grid (or vector) based on the value from another grid 
     inputs:
          VIN = pointer to the vector being scaled
          Scalar = pointer to the vector of scaling weights
     outputs:
          OP = vector containing the values of VIN scaled by values in Scalar
*/
std::vector<double> GridScaleByGrid(std::vector<double> *VIN, std::vector<double> *Scalar);


/*  Takes a grid (or vector) and normalizes by the range. So, after the operation the min value = 0 and the max value = 1
     inputs:
          VIN = pointer to the grid or vector to be normalized
     outputs: NA
*/
void NormalizeGrid(std::vector<double> *VIN);


/* Structure to be output by gridDerivatives, contains original data, 1st & 2nd derivatives and sigmas for gaussian blurring used
     members:
          Raw = Original grid being operated on (unsmoothed)
          S1 = sigma to be used in gaussian smoothing prior to taking the 1st derivative.. if S1 == 0, no smoothing is done prior to calculating gradient
          Gradient = first derivative of Raw after performing gaussian smoothing with sigma = S1.. Calculated using the function sobel in OpenCV
          S2 = sigma to be used prior to taking the laplacian of the origianl grid.. if S2 == 0, no smoothing is done prior to calculating the laplacian
          Laplacian = 2nd derivative of Raw after performing guassian smoothing with sigma = S2.. Calculated using the laplacian function in OpenCV
          GS = size of grid cells for all grids
          DL = delineation length (number of x grid cells (i.e., number of colums))
*/
struct GridDerivatives{
     std::vector<double> Raw;
     std::vector<double> Gradient;
     std::vector<double> Laplacian;
     std::vector<double> GradGB;
     std::vector<double> LapGB;
     double GS;
     uint32_t DL;
     uint32_t S1;
     uint32_t S2;
};


/* Takes an input grid (projected in 1D) and calculates 1st and 2nd derivatives... optional whether to perform gaussian smoothing prior to derivative calculation
     inputs:
          VIN = a pointer to the vector representing the grid to be operated on
          DL = delineation length (number of x grid cells... number of columns..)
          GS = size of grid cells in VIN
          Translate = optional translation value for grid cell values prior to derivative calculation (i.e., this value is subtracted from each grid cell)
          Gsig = sigma to be used for gaussian smoothing prior to calculating gradient.. if Gsig == 0, no smoothing occurs for gradient
          Lsig = sigma to be used for gaussian smoothing prior to calculating laplacian .. if Lsig == 0, no smoothing occurs for laplacian

     outputs:
          OP = GridDerivatives object containing the original data (after optional translation), 1st & 2nd derivatives, and meta data about the grid and smoothing used prior to derivative calculation
*/
GridDerivatives gridDerivatives(std::vector<double> *VIN, uint32_t DL, double GS, double Translate, uint32_t Gsig, uint32_t Lsig);


#endif
