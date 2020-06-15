#ifndef ImProc_H
#define ImProc_H

#include<iostream> 
#include<vector>
#include<stdint.h>
#include<algorithm>
#include<cmath>


/* Takes matrix (projected into 1D) and performs a gaussian convolution 
        inputs: Grid = a pointer to the vector (matrix)
                Sigma = the standard deviation for the guassian kernel 
                DL = Deliniation length, the number of cells in the x grid
                thresh = if this value IS NOT 0, all points above the thresh hold will recieve one value and all under will recieve another
                impute = the value to be given to cells without a value 
        outputs: OP = convolved grid values (2D projected to 1D)
*/
std::vector<double> Gaus2D(std::vector<double> *Grid, double sigma, uint32_t DL, double thresh, double impute);


/* structure to be output by ImGrad
        X = magnitude of the gradient in the X direction
        Y = magnitude of the gradient in the Y deriection
        Mag = total gradient
        Dir = direction of the gradient
        DDir = discretized dircetion of the gradient, for non-maxima supression

*/
struct ImGradOP{
        std::vector<double> X;
        std::vector<double> Y;
        std::vector<double> Mag;
        std::vector<double> Dir;
        std::vector<uint16_t> DDir;
        int32_t DL;
};


/* takes pointer to a vector representing an image (or raster) (projected in 1D) and ruturns several values related to the gradient of the image (in 1D)
        inputs: VIN = pointer to the vector (matrix) to be processed
                DL = delineation length, the number of cells in the x grid
                thresh = if thresh DOES NOT equal 0, all elements above this value will be asigned one value and all under assigned another value
        output: OP = ImGradOP object containing members of gradient magnitude in the x direction, y direction, total magnitude, direction of gradient, and discretized direction

*/
ImGradOP ImGrad(std::vector<double> *VIN, int32_t DL, bool th, double thresh);

/* takes a matrix of image gradient that has been thesholded but retains original gradients in the retained areas


*/
std::vector<bool> NMS(ImGradOP *TG);



/* Takes matrix (projected into 1D) and performs a laplacian of gaussain (LoG) convolution 
        inputs: Lasf = a pointer to the LasFile from which ZM came
                Gsize = the grid cell size used to produce ZM
                clean = bool for if cleaning operation should be perfomred
                tol = the tolerance to be used to clean (meters)
        outputs: OP = convolved grid values (2D projected to 1D)
*/
//&& (  floor((((yg*DL) + (j * DL)) + (xg + k) )/DL) == yg )
std::vector<double> LoG2D(std::vector<double> *Grid, double sigma, uint32_t DL, bool th, double thresh, double impute);


/* Given a standard deviation for the gaussian, creates a LoG convolution kernel that is reweighted to sum to (nearly) zero
	input: sigma = the standard deviation of the gaussian
	output: kernel = 2D vector (projected in 1D) of weights for a LoG kernel for the given sigma, dimensions are 4*SD x 4*SD
*/
std::vector<double> LoGRW(int32_t sigma);

/*structure to be output by Conv2D 
	members: grid = the result of the convolution projected in 1D
		DL = the delineation length of the 1D projected output grid
*/
struct Conv2DOP{
	std::vector<double> grid;
	uint32_t DL;
};

/* convolves an image (grid) projected in 1D with a given kernel projected in 1D from pointers
	inputs: kernel = a pointer to the kernel to convolve grid with, projected in 1D, needs to be a square
		grid = the grid (or image) projected in 1D to be convolved
		DL = delineation length of the 1D projected grid
		impute = the value used to impute grid cells without a value or areas outside the grid extent
	outputs: OP = structure containing a vector of doubles that represent the result of the convolution,(2D grid projected in 1D) and an uint32_t which designates the delineation length of the output grid
*/
Conv2DOP Conv2D(std::vector<double> *grid, uint32_t DL,  std::vector<double> kernel, double impute);



/* Takes matrix (projected into 1D) and performs a laplacian of gaussain (LoG) convolution 
        inputs: Lasf = a pointer to the LasFile from which ZM came
                Gsize = the grid cell size used to produce ZM
                clean = bool for if cleaning operation should be perfomred
                tol = the tolerance to be used to clean (meters)
        outputs: OP = convolved grid values (2D projected to 1D)
*/
//&& (  floor((((yg*DL) + (j * DL)) + (xg + k) )/DL) == yg )
/*
std::vector<double> LoG2DRW(std::vector<double> *Grid, uint32_t DL, std::vector<double>  bool th, double thresh, double impute){
        std::vector<double> OP(Grid->size());
        int32_t yg, xg;
        for(size_t i = 0; i < Grid->size(); i++){
                yg = floor(i / DL);
                xg = i - yg*DL;i
		uint32_t corr = ((DLK - 1) / 2 ) + 1;
                for(size_t j = -DLK; j < DLK + 1; j++){
                        for(size_t k = -DLK; k < DLK + 1; k++){
                                if( ( (((yg*DL) + (j * DL)) + (xg + k) ) < OP.size()) && ( (((yg*DL) + (j * DL)) + (xg + k) ) > 0)  &&  ( floor((((yg*DL) + (j * DL)) + (xg + k) )/DL) == (yg + j) ) ){
                                        if(Grid->at(((yg*DL) + (j * DL)) + (xg + k)) != 0){
                                                OP[i] = OP[i] + Grid->at(((yg*DL) + (j * DL)) + (xg + k)) * kernel->at(((j + corr) * DLK) + (k + corr));
                                        } else {
                                        OP[i] = OP[i] + (impute) * kernel->at(((j + corr) * DLK) + (k + corr));
                                        }
                                } else {
					OP[i] = OP[i] + (impute) * kernel->at(((j + corr) * DLK) + (k + corr));
				}
                        }
                }
        }
	if(th){
	        if(thresh){
        	        for(uint32_t i = 0; i < OP.size(); i++){
                	        if(OP[i] < thresh){
                        	        OP[i] = 100;
	                        } else {
        	                        OP[i] = 0;
                	        }
               		} 
		} else {
			for(uint32_t i = 0; i < OP.size(); i++){
                	        if(OP[i] > thresh){
                        	        OP[i] = 0;
	                        }
               		}
		}
        }
        return(OP);
}
*/

#endif
