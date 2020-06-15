#include<iostream> 
#include<vector>
#include<stdint.h>
#include<algorithm>
#include<cmath>
#include<PCProc.h>


/*###################################################

Everything below is outdated and replaced with funcionality in OpenCV

###################################################*/

/* Takes matrix (projected into 1D) and performs a gaussian convolution 
        inputs: Grid = a pointer to the vector (matrix)
                Sigma = the standard deviation for the guassian kernel 
                DL = Deliniation length, the number of cells in the x grid
                thresh = if this value IS NOT 0, all points above the thresh hold will recieve one value and all under will recieve another
                impute = the value to be given to cells without a value 
        outputs: OP = convolved grid values (2D projected to 1D)
*/
//&& (  floor((((yg*DL) + (j * DL)) + (xg + k) )/DL) == yg ) 
std::vector<double> Gaus2D(std::vector<double> *Grid, double sigma, uint32_t DL, double thresh, double impute){
        std::vector<double> OP(Grid->size());
        int32_t yg, xg;
        for(uint32_t i = 0; i < Grid->size(); i++){
                yg = floor(i / DL);
                xg = i - yg*DL;
                for(int32_t j = -4*sigma; j < 4*sigma + 1; j++){
                        for(int32_t k = -4*sigma; k < 4*sigma + 1; k++){
                                if( ( (((yg*DL) + (j * DL)) + (xg + k) ) < OP.size()) && ( (((yg*DL) + (j * DL)) + (xg + k) ) > 0) &&  ( floor((((yg*DL) + (j * DL)) + (xg + k) )/DL) == (yg + j) )   ){
                                        if(Grid->at(((yg*DL) + (j * DL)) + (xg + k)) != 0){
                                                OP[i] = OP[i] + ((Grid->at(((yg*DL) + (j * DL)) + (xg + k))) * (1/(2*M_PI*(sigma*sigma))) * exp(-((k*k + j*j) / (2*sigma*sigma))));
                                        } else {
                                                OP[i] = OP[i] + impute * (1/(2*M_PI*(sigma*sigma))) * exp(-((k*k + j*j) / (2*sigma*sigma)));
                                        } 					
				} else {
					OP[i] = OP[i] + impute * (1/(2*M_PI*(sigma*sigma))) * exp(-((k*k + j*j) / (2*sigma*sigma)));
                                }
                        }
                }
        }
        if(thresh){
                for(uint32_t i = 0; i < OP.size(); i++){
                        if(OP[i] < thresh){
                                OP[i] = 100;
                        } else {
                                OP[i] = 0;
                        }
                }
        }
        return(OP);
}


/* takes pointer to a vector representing an image (or raster) (projected in 1D) and ruturns several values related to the gradient of the image (in 1D)
        inputs: VIN = pointer to the vector (matrix) to be processed
                DL = delineation length, the number of cells in the x grid
                thresh = if thresh DOES NOT equal 0, all elements above this value will be asigned one value and all under assigned another value
        output: OP = ImGradOP object containing members of gradient magnitude in the x direction, y direction, total magnitude, direction of gradient, and discretized direction

*/
ImGradOP ImGrad(std::vector<double> *VIN, int32_t DL, bool th, double thresh){
        ImGradOP OP;
	OP.DL = DL;
        OP.X.resize((VIN->size()));
        OP.Y.resize((VIN->size()));
        OP.Mag.resize(VIN->size());
        OP.Dir.resize(VIN->size());
        OP.DDir.resize(VIN->size());
        int32_t tmp = VIN->size();
        double PQ = 2*M_PI/16;
        for(int32_t i = 1; i < tmp - 1; i++){
                if( ((i + DL + 1) <  tmp) && ((i - DL - 1) > 0 ) && (floor((i + DL + 1) / DL) == floor((i + DL) / DL )) && (floor((i - DL - 1) / DL) == floor((i - DL) / DL ))){
			OP.X[i] = VIN->at(i + 1)*2 - VIN->at(i - 1)*2 + VIN->at(i + DL + 1) - VIN->at(i + DL - 1) + VIN->at(i - DL + 1) - VIN->at(i - DL - 1);
                        OP.Y[i] = VIN->at(i + DL)*2  + VIN->at(i + DL + 1) + VIN->at(i + DL - 1) - VIN->at(i - DL)*2 - VIN->at(i - DL + 1) - VIN->at(i - DL - 1);
                }
                OP.Mag[i] = pow(OP.X[i]*OP.X[i] + OP.Y[i]* OP.Y[i], .5);
                OP.Dir[i] = atan2((OP.Y[i]),(OP.X[i]));
                if(OP.Y[i] < 0){
                        OP.Dir[i] = OP.Dir[i] + 2*M_PI;
                }
                if(((OP.Dir[i] >= 0) && (OP.Dir[i] < PQ)) || (OP.Dir[i] >= 15*PQ) || (OP.Dir[i] >= 7*PQ && (OP.Dir[i] < 9*PQ))){
                        OP.DDir[i] = 1;
                } else if(((OP.Dir[i] >= PQ) && (OP.Dir[i] < 3*PQ)) || ((OP.Dir[i] >= 9*PQ) && (OP.Dir[i] < 11*PQ))) {
                        OP.DDir[i] = 2;
                } else if((((OP.Dir[i]) >= 3*PQ) && ((OP.Dir[i]) < 5*PQ)) || ((OP.Dir[i] >= 11*PQ) && (OP.Dir[i] < 13*PQ))){
                        OP.DDir[i] = 3; 
                } else if(((OP.Dir[i] >= 5*PQ) && (OP.Dir[i] < 7*PQ)) || ((OP.Dir[i] >= 13*PQ) && (OP.Dir[i] < 16*PQ))){
                        OP.DDir[i] = 4;
                }
        }
	if(thresh){
	        if(th){
        	        for(uint32_t i = 0; i < OP.Mag.size(); i++){
                	        if(OP.Mag[i] > thresh){
                        	        OP.Mag[i] = 100;
	                        } else {
        	                        OP.Mag[i] = 0;
                	        }
               		} 
		} else {
			for(uint32_t i = 0; i < OP.Mag.size(); i++){
                	        if(OP.Mag[i] < thresh){
                        	        OP.Mag[i] = 0;
	                        }
               		}
		}
        }
        return(OP);
}


/* takes a matrix of image gradient that has thesholded but retains original gradients in the retained areas


*/
std::vector<bool> NMS(ImGradOP *TG){
	std::vector<bool> OP(TG->Mag.size());
	for(uint32_t i = 0; i < OP.size(); i++){
		//std::cout << "iteration # " << i << std::endl;
		if(TG->Mag[i]){
			switch(TG->DDir[i]){
				case 1:
					if( (TG->Mag[i] < TG->Mag[i + 1] ) || (TG->Mag[i] < TG->Mag[i + 1]) ){
						TG->Mag[i] = 0;
					} else {
						OP[i] = true;
			//			std::cout << i <<  " Max dir = " <<  1 << std::endl;
					}
					break;
				case 2:
					if( (TG->Mag[i] < TG->Mag[i + TG->DL + 1]) || (TG->Mag[i] < TG->Mag[i - TG->DL - 1]) ){
						TG->Mag[i] = 0;
					} else {
						OP[i] = true;
			//			std::cout << i << " Max dir = " <<  2 << std::endl;
					}
					break;
				case 3:
					if( (TG->Mag[i] < TG->Mag[i + TG->DL]) || (TG->Mag[i] < TG->Mag[i -TG-> DL]) ){
						TG->Mag[i] = 0;
					} else {
						OP[i] = true;	
			//			std::cout << i <<  " Max dir = " <<  3 << std::endl;
					}
					break;
				case 4:
					if( (TG->Mag[i] < TG->Mag[i + TG->DL - 1]) || (TG->Mag[i] < TG->Mag[i - TG->DL + 1]) ){
						TG->Mag[i] = 0;
					} else {
						OP[i] = true;
			//			std::cout << i <<  " Max dir = " <<  4 << std::endl;
					}
					break;
				default:
					std::cout << "!! ERROR DIRECTION NOT VALID !!" << std::endl;
					break;

			}
		}
	} 
	return(OP);
}



/* Takes matrix (projected into 1D) and performs a laplacian of gaussain (LoG) convolution 
        inputs: Lasf = a pointer to the LasFile from which ZM came
                Gsize = the grid cell size used to produce ZM
                clean = bool for if cleaning operation should be perfomred
                tol = the tolerance to be used to clean (meters)
        outputs: OP = convolved grid values (2D projected to 1D)
*/
//&& (  floor((((yg*DL) + (j * DL)) + (xg + k) )/DL) == yg )
std::vector<double> LoG2D(std::vector<double> *Grid, double sigma, uint32_t DL, bool th, double thresh, double impute){
        std::vector<double> OP(Grid->size());
        int32_t yg, xg;
        for(uint32_t i = 0; i < Grid->size(); i++){
                yg = floor(i / DL);
                xg = i - yg*DL;
                for(int32_t j = -4*sigma; j < 4*sigma + 1; j++){
                        for(int32_t k = -4*sigma; k < 4*sigma + 1; k++){
                                if( ( (((yg*DL) + (j * DL)) + (xg + k) ) < OP.size()) && ( (((yg*DL) + (j * DL)) + (xg + k) ) > 0)  &&  ( floor((((yg*DL) + (j * DL)) + (xg + k) )/DL) == (yg + j) ) ){
                                        if(Grid->at(((yg*DL) + (j * DL)) + (xg + k)) != 0){
                                                OP[i] = OP[i] + ((Grid->at(((yg*DL) + (j * DL)) + (xg + k)))) * ((-1/(M_PI*pow(sigma,4))) * (1 - ((k*k + j*j) / (2*sigma*sigma))) * exp(-((k*k + j*j) / (2*sigma*sigma))));
                                        } else {
                                        OP[i] = OP[i] + (impute) * ((-1/(M_PI*pow(sigma,4))) * (1 - ((k*k + j*j) / (2*sigma*sigma))) * exp(-((k*k + j*j) / (2*sigma*sigma))));
                                        }
                                } else {
					OP[i] = OP[i] + (impute) * ((-1/(M_PI*pow(sigma,4))) * (1 - ((k*k + j*j) / (2*sigma*sigma))) * exp(-((k*k + j*j) / (2*sigma*sigma))));
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


/* Given a standard deviation for the gaussian, creates a LoG convolution kernel that is reweighted to sum to (nearly) zero
	input: sigma = the standard deviation of the gaussian
	output: kernel = 2D vector (projected in 1D) of weights for a LoG kernel for the given sigma, dimensions are 4*SD x 4*SD
*/
	

std::vector<double> LoGRW(int32_t sigma){
        double pos,val,sum;
	sum = 0; pos = 0;
        size_t ctr = 0;
        for(int32_t j = -4*sigma; j < 4*sigma + 1; j++){
                for(int32_t k = -4*sigma; k < 4*sigma + 1; k++){
                        val = ((-1/(M_PI*pow(sigma,4))) * (1 - ((k*k + j*j) / (2*sigma*sigma))) * exp(-((k*k + j*j) / (2*sigma*sigma))));
                        sum = sum + val;
                        if(val > 0){
                                pos++;
                        }
                        ctr++;
                }
        }
        double cor = sum/(pos);     
        std::vector<double> kernel(ctr);
        ctr = 0;
        for(int32_t j = -4*sigma; j < 4*sigma + 1; j++){
                for(int32_t k = -4*sigma; k < 4*sigma + 1; k++){
                        val = ((-1/(M_PI*pow(sigma,4))) * (1 - ((k*k + j*j) / (2*sigma*sigma))) * exp(-((k*k + j*j) / (2*sigma*sigma))));
                        if(val > 0){
                                sum = sum + val - cor;
                                kernel[(j+4*sigma)*8*sigma + (k+4*sigma) + ctr] = val - cor;
                        } else {
                                kernel[(j+4*sigma)*8*sigma + (k+4*sigma) + ctr] = val;
                                sum = sum + val;
                        }
                }
                ctr++;
        }
	return(kernel);
}


/* convolves an image (grid) projected in 1D with a given kernel projected in 1D from pointers
	inputs: kernel = a pointer to the kernel to convolve grid with, projected in 1D, needs to be a square
		grid = the grid (or image) projected in 1D to be convolved
		DL = delineation length of the 1D projected grid
		impute = the value used to impute grid cells without a value or areas outside the grid extent
	outputs: OP = structure containing a vector of doubles that represent the result of the convolution,(2D grid projected in 1D) and an uint32_t which designates the delineation length of the output grid
*/
Conv2DOP Conv2D(std::vector<double> *grid, uint32_t DL,  std::vector<double> kernel, double impute){
	Conv2DOP OP;
	OP.grid.resize(grid->size());
	OP.DL = DL;
	int32_t EZI = (grid->size()/2)-1;
	uint32_t xg, yg;
	for(size_t i = 0; i < grid->size(); i++){
		yg = floor(i / DL);
                xg = i - yg*DL;
		for(int32_t j = -EZI; j < EZI; j++){
			for(int32_t k = -EZI; k < EZI; k++){
				
					
				//			HERE!!!!!!!!!!!!!!!!!1
			}
		} 
	}

	return(OP);
}






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


