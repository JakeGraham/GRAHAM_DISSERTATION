#include<iostream>
#include<iomanip>
#include<stdint.h>
#include<math.h>
#include<Misc.h>
#include<PCProc.h>
#include<LasReadWrite.h>
#include<SPRUCE.h>
#include<stdio.h>
#include<opencv2/imgproc.hpp>
#include<opencv2/imgcodecs.hpp>
#include<opencv2/highgui.hpp>
#include<OpenCVComp.h> 
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
Mat GridToCVMat(std::vector<double> *Grid, int32_t DL, bool translate , double transV){  
     if(translate){
          for(size_t i = 0; i < Grid->size(); i++){ 
               if(Grid->at(i) != 0){ 
                    Grid->at(i) = (Grid->at(i) - transV); 
               } 
          }
     } 
     size_t nrow = DL;  // rows and colums are swapped in openCV... I think, still not clear
     size_t ncol = Grid->size()/DL;  // rows and columns are swapped in openCV
     Mat OP = Mat::zeros(ncol, nrow, CV_64F); // rows and columns are swapped in openCV
     memcpy(OP.data, &Grid->at(0), Grid->size()*sizeof(double)); 
     return(OP);
}


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
Mat GridToCVMat_3C(std::vector<double> *Grid, int32_t DL, bool translate , double transV){  
     if(translate){
          for(size_t i = 0; i < Grid->size(); i++){ 
               if(Grid->at(i) != 0){ 
                    Grid->at(i) = (Grid->at(i) - transV); 
               } 
          }
     } 
     size_t nrow = DL;  // rows and colums are swapped in openCV... I think, still not clear
     size_t ncol = Grid->size()/DL;  // rows and columns are swapped in openCV
     Mat OP = Mat::zeros(ncol, nrow, CV_64FC(3)); // rows and columns are swapped in openCV
     memcpy(OP.data, &Grid->at(0), Grid->size()*sizeof(double)); 
     memcpy(&OP.data[Grid->size()*sizeof(double)], &Grid->at(0), Grid->size()*sizeof(double)); 
     memcpy(&OP.data[2*(Grid->size()*sizeof(double))], &Grid->at(0), Grid->size()*sizeof(double)); 
     return(OP);
}


/* Converts a cv::Mat object to a 2D grid projected in 1D. Currently only works with a single band 64 bit cv::Mat image, will likely modify to accept multiband images and images with data types other than CV_64F (double).
     inputs: 
               MT = a pointer to the cv::Mat object to be converted
     outputs:
               OP = MTG object containing the vector (grid) and DL (delineation length i.e., the number of x grid cells)
*/
MTG CVMatToGrid(Mat *MT){
     MTG OP;
     OP.Grid.resize(MT->cols*MT->rows);
     OP.DL = MT->cols; 
     memcpy(&OP.Grid[0], MT->data, OP.Grid.size()*sizeof(double));
     return(OP);
}


/* Creates a new Mat object that only contains values from the "IM" Mat where the "Mask" Mat is true (i.e., != 0)
     inputs: 
              IM =  pointer to Mat object from which values will be extracted for output
              Mask =  pointer to Mat object that will be used to "mask" IM
     outputs: 
              OP = Mat object containing values from IM where Mask != 0
              OPF = empty Mat object that is returned if dimensions of IM and Mask don't match

 !!!!!!!!!!!!!!!! NOTE ~~~~~~~~~~ Only works with Mat objects storing doubles.. may adjust to be generic.. not sure how ATM
*/
Mat MaskMat(Mat *IM, Mat *Mask){
     if(IM->cols != Mask->cols | IM->rows != Mask->rows){
          std::cout << "Warning! Image Masking failed!!!! Input image and Mask do not have the same dimensions!!!!" << std::endl;
          Mat OPF;
          return(OPF);
     } 
     Mat OP = Mat::zeros(IM->rows, IM->cols, IM->type());
     for(size_t i = 0; i < IM->rows; i++){
          for(size_t j = 0; j < IM->cols; j++){
               if(Mask->at<double>(i,j)){
                    OP.at<double>(i,j) = IM->at<double>(i,j);
               }
          }
     }
     return(OP);
}


/* Calculates the image gradient using sobel operators in OpenCV to a std::vector<double>
     inputs:
              M = pointer to Mat object to calculate the gradient of 
              Mat = if Mat == true, the output is an object of type CV::Mat, if false output is std::vector<double>
     outputs:
              OP = 2d grid (projected in 1D) containing absolute gradient from M
*/
std::vector<double> CVImGradGrid(Mat *M){
     std::vector<double> OP(M->cols*M->rows);
     Mat sobx, soby;
     Sobel(*M, sobx, 6, 1, 0, 3);
     Sobel(*M, soby, 6, 0, 1, 3);
     MTG Xgrad = CVMatToGrid(&sobx);
     MTG Ygrad = CVMatToGrid(&soby);
     for(size_t i = 0; i < Xgrad.Grid.size(); i++){
          OP[i] = pow(pow(Xgrad.Grid[i], 2) + pow(Ygrad.Grid[i], 2), .5);
     }
     return(OP);
}


/*  Scales the values of one grid (or vector) based on the value from another grid 
     inputs:
          VIN = pointer to the vector being scaled
          Scalar = pointer to the vector of scaling weights
     outputs:
          OP = vector containing the values of VIN scaled by values in Scalar
*/
std::vector<double> GridScaleByGrid(std::vector<double> *VIN, std::vector<double> *Scalar){
     if(VIN->size() != Scalar->size()){
          std::cout << "!!!!!! ~~~     WARNING! Grid scaling failed! input vectors differ in length ~~~~ !!!!!!!" << std::endl;
          std::vector<double> OPF;
          return(OPF);
     }
     std::vector<double> OP(VIN->size());
     for(size_t i = 0; i < OP.size(); i++){
          OP[i] = VIN->at(i) * (Scalar->at(i)); // TESTTING HERE
     }
     return(OP);
}

/*  Takes a grid (or vector) and normalizes by the range. So, after the operation the min value = 0 and the max value = 1
     inputs:
          VIN = pointer to the grid or vector to be normalized
     outputs: NA
*/
void ScaleGridZeroOne(std::vector<double> *VIN){
     double Max = *max_element(std::begin(*VIN), std::end(*VIN));
     double Min = *min_element(std::begin(*VIN), std::end(*VIN));
     double range = Max - Min;
     for(size_t i = 0; i < VIN->size(); i++){
          VIN->at(i) = range / ((VIN->at(i) - Min)) ;
     }
}


/*  Takes a grid (or vector) and normalizes by the range. So, after the operation the min value = 0 and the max value = 1
     inputs:
          VIN = pointer to the grid or vector to be normalized
     outputs: NA
*/
void NormalizeGrid(std::vector<double> *VIN){
     double Max = *max_element(std::begin(*VIN), std::end(*VIN));
     double Min = *min_element(std::begin(*VIN), std::end(*VIN));
     double mean = accumulate( VIN->begin(), VIN->end(), 0.0)/VIN->size(); 
     double range = Max - Min;
     std::vector<double> diff(VIN->size());
     std::transform(VIN->begin(), VIN->end(), diff.begin(), [mean](double x) { return x - mean; });
     double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
     double stdev = std::sqrt(sq_sum / VIN->size());
     for(size_t i = 0; i < VIN->size(); i++){
          VIN->at(i) = (VIN->at(i) - mean) / stdev;
     }
}


/* Takes an input grid (projected in 1D) and calculates 1st and 2nd derivatives... optional whether to perform gaussian smoothing prior to derivative calculation
     inputs:
          VIN = a pointer to the vector representing the grid to be operated on
          DL = delineation length (number of x grid cells... number of columns..)
          GS = size of grid cells in VIN
          Translate = optional translation value for grid cell values prior to derivative calculation (i.e., this value is subtracted from each grid cell)
          Gsig = sigma to be used for gaussian smoothing prior to calculating gradient.. if Gsig == 0, no smoothing occurs for gradient
          Lsig = sigma to be used for gaussian smoothing prior to calculating laplacian .. if Lsig == 0, no smoothing occurs for laplacian

     outputs:
          OP = GridDerivatives object containing the original data (after optional translation), 1st & 2nd derivatives, and meta data about the grid and smoothing u
sed prior to derivative calculation
*/
GridDerivatives gridDerivatives(std::vector<double> *VIN, uint32_t DL, double GS, double Translate, uint32_t Gsig, uint32_t Lsig){
     GridDerivatives OP; // create output

     // put elevation grid into OpenCV Mat object
     Mat MT = GridToCVMat(VIN, DL, true, Translate);
     std::vector<double> Grid = CVMatToGrid(&MT).Grid;
     Mat LGB, GGB;
     int16_t WSSD = 9;

     // Make sure sigmas are odd numbers and perform smoothing
     // gaussian for gradient
     if(Gsig % 2 == 0){
          GaussianBlur(MT, GGB, Size(Gsig*WSSD+1, Gsig*WSSD+1), Gsig, Gsig);
     } else {
          GaussianBlur(MT, GGB, Size(Gsig*WSSD, Gsig*WSSD), Gsig, Gsig);
     }
     // gaussian for laplacian
     if(Lsig % 2 == 0){
          GaussianBlur(MT, LGB, Size(Lsig*WSSD+1, Lsig*WSSD+1), Lsig, Lsig);
     } else {
          GaussianBlur(MT, LGB, Size(Lsig*WSSD, Lsig*WSSD), Lsig, Lsig);
     }

     // calculate gradient 
     std::vector<double> GV = CVImGradGrid(&GGB);

     // calculate laplacian
     Mat lap;
     Laplacian(LGB, lap, LGB.depth(), 3);

     OP.DL = DL;
     OP.GS = GS;
     OP.Raw = *VIN;
     OP.S1 = Gsig;
     OP.Gradient = GV;
     OP.S2 = Lsig;
     OP.Laplacian = CVMatToGrid(&lap).Grid;
     OP.GradGB = CVMatToGrid(&GGB).Grid;
     OP.LapGB = CVMatToGrid(&LGB).Grid;
     return(OP);
}


/* converts gradient in length/length to degrees
     inputs:
          VIN = pointer to the vector containing gradient values
          GS = the length of a grid cell in VIN
     outputs: NA
*/
void GridGradientToDegrees(std::vector<double> *VIN, double GS){
     for(size_t i = 0; i < VIN->size(); i++){
          VIN->at(i) = ( atan(GS / VIN->at(i)) * 180 ) / M_PI; 
     }  
}



