#include<iostream>
#include<string>
#include<fstream>
#include<stdint.h>
#include<vector>
#include<LasReadWrite.h>
#include<Misc.h>
#include<PCProc.h>
#include<iomanip>
#include<SPRUCE.h>
#include<OpenCVComp.h>
#include<opencv2/highgui.hpp>




float SPRUCE_TempT(int32_t PlotID){
   float TempT;
   if(PlotID == 7 | PlotID == 21){
      TempT = -1;
   }
   else if(PlotID == 6 || PlotID == 19){
      TempT = 0;
   }
   else if(PlotID == 11 || PlotID == 20){
      TempT = 2.25;
   }
   else if(PlotID == 4 || PlotID == 13){
      TempT = 4.5;
   }
   else if(PlotID == 8 || PlotID == 16){
      TempT = 6.75;
   }
   else if(PlotID == 10 || PlotID == 17){
      TempT = 9;
   }
   return(TempT);
}




int32_t SPRUCE_CO2T(int32_t PlotID){
   int32_t CO2T;
   if(PlotID == 13 || PlotID == 17 || PlotID == 6 || PlotID == 8 || PlotID == 20){
      CO2T = 0;
   } else if(PlotID == 16 || PlotID == 19 || PlotID == 4 || PlotID == 11 || PlotID == 10){
      CO2T = 1;
   } else if(PlotID == 7 || PlotID == 21){
      CO2T = -1;
   }
   return(CO2T);
}











/* takes a sting of a SPRUCE .las file and parses by "_" and "." for using year, plot etc...
   inputs: SIN = a string of the .las file name
   outputs: OP = vector of strings that are separated by "_" and "."
*/
std::vector<std::string> ParseSPRUCELas(std::string SIN){
   std::vector<std::string> OP;
//        std::cout << SIN << std::endl;
        std::vector<std::string> parsed = SplitString(SIN, "_");
        std::vector<std::string> par2 = SplitString(parsed.back(), ".");
   for(uint32_t i = 0; i < parsed.size() - 1; i++){
      OP.push_back(parsed[i]);
   }
   for(uint32_t i = 0; i < par2.size(); i++){
      OP.push_back(par2[i]);
   }
   return(OP);
}

// Takes boardwalk UTMs and returns the mean distance between opposing sides (2D)
double BWML(std::vector<double> BWC){
/* Order of BW coords      X1       X2     Y1       Y2     of opposing sides
                        (BWC[0], BWC[12], BWC[1], BWC[13]);
                        (BWC[3], BWC[15], BWC[4], BWC[16]);
                        (BWC[6], BWC[18], BWC[7], BWC[19]);
                        (BWC[9], BWC[21], BWC[10], BWC[22]);
*/

        std::vector<double> LS( (BWC.size() / 3) / 2);
        for (uint16_t i = 0; i < BWC.size() / 2; i+= 3){
                // line below is pretty hard to understand.. consider breaking into components, possibly make function to calc hypotenuse
                LS[i / 3]  = pow( (pow( (BWC[i] - BWC[i + (BWC.size() / 2)]),2) + pow( (BWC[i + 1] - BWC[i + 1 + (BWC.size() / 2)]),2) ), .5);
        }
        double ML = 0;
        for (uint16_t i = 0; i < LS.size(); i++){
                ML += LS[i];
        }
        ML = ML / LS.size();
        std::cout << "Mean L = " << ML << std::endl;    // currently give mean distance from opposite sides of boardwalk
        return(ML);
}


// Takes boardwalk UTMs and and returns UTMs of the 8 vertices
std::vector<double> BWVerts(std::vector<double> BWC){
        std::vector<double> Vout(4);
        Vout[0] = BWC[0] - BWC[12]; // X component of line between two opposing sides
        Vout[1] = BWC[1] - BWC[13]; // Y component of line between two opposing sides
        Vout[2] = BWC[3] - BWC[15]; // X component of line between other two opposing sides
        Vout[3] = BWC[4] - BWC[16]; // Y component of line between other two opposing sides
        double CX =  BWC[12] + (Vout[0] / 2);
        double CY = BWC[13] + (Vout[1] / 2);
        std::vector<double> Corners(16);

        Corners[0] = BWC[0] + (Vout[1] / 2);
        Corners[1] = BWC[1] - (Vout[0] / 2); // 0
        Corners[2] = BWC[0] - (Vout[1] / 2);
        Corners[3] = BWC[1] + (Vout[0] / 2); // 1
        Corners[4] = BWC[12] + (Vout[1] / 2);
        Corners[5] = BWC[13] - (Vout[0] / 2); // 2
        Corners[6] = BWC[12] - (Vout[1] / 2);
        Corners[7] = BWC[13] + (Vout[0] / 2); // 3
        Corners[8] = BWC[3] + (Vout[3] / 2);
        Corners[9] = BWC[4] - (Vout[2] / 2); // 4
        Corners[10] = BWC[3] - (Vout[3] / 2);
        Corners[11] = BWC[4] + (Vout[2] / 2); // 5
        Corners[12] = BWC[15] + (Vout[3] / 2);
        Corners[13] = BWC[16] - (Vout[2] / 2); // 6 
        Corners[14] = BWC[15] - (Vout[3] / 2);
        Corners[15] = BWC[16] + (Vout[2] / 2); // 7

        std::vector<double> IP;
        double dist2cent;
        std::vector<double> Verts;
        int x1, x2, x3, x4;
        for(uint16_t j = 0; j < Corners.size() / 2; j += 2){
                x1 = j;
                for (uint16_t j2 = j + 2; j2 < Corners.size() / 2; j2 += 2){
                        x2 = j2;
                        for (uint16_t i = Corners.size() / 2; i < Corners.size(); i += 2){
                                x3 = i;
                                for (uint16_t i2 = i + 2; i2 < Corners.size(); i2 += 2){
                                        x4 = i2;
                                        IP = IntrsctPt(Corners[x1],Corners[x1 + 1],Corners[x2],Corners[x2 + 1],Corners[x3],Corners[x3 + 1],
                                                Corners[x4],Corners[x4 + 1]);
                                        dist2cent = DP2P(CX, CY, IP[0], IP[1]);
                                        if(dist2cent < 5 && dist2cent > 4.5){
                                                Verts.push_back(IP[0]); Verts.push_back(IP[1]);
                                        }
                                }
                        }
                }
        }

        int inds[] = {0,1,4,5,7,6,3,2};
        std::vector<double> Cverts(Verts.size());

/*      std::ofstream ofile;
        ofile.open("testcloud.txt");
*/
        for(uint32_t i = 0; i < Verts.size(); i+=2){
                //ofile << std::setprecision(9) << i/2 << " " << Verts[inds[i/2]*2] << " " << Verts[inds[(i/2)]*2 + 1] <<  " "  << 420 << std::endl;
                Cverts[i] = Verts[inds[i/2]*2]; 
      Cverts[i+1] = Verts[inds[i/2]*2 + 1];
        }

//      ofile.close();
        return(Cverts);
}



/* Finds the shortest distance from a point to a SPRUCE boardwalk
        inputs: BWV = a pointer to a vector of boardwalk vectices.. from BWVerts 
*/
double BWD (double *BWV, double *p){
   double mindist = 25;//  D2L((BWV + 15), BWV, p);
   for(int16_t i = 0; i < 15; i += 2){
      double dist = D2L((BWV + i), (BWV + i + 2), p);
      if(mindist > dist) {
         mindist = dist;
      }
   }
   return(mindist);
}


/* Calculates the distance to the boardwalk for all points in a LasFile
        inputs: LF = pointer to LasFile to be operated on
                BMV = pointer to vector or array of boardwalk vertices in the order x1,y1,x2,y2,.. etc
        output: OPV = A vector of boardwalk distance for each point in LF in order

*/
std::vector<double> LasBWD(LasFile *LF, double *BWV){
   std::vector<double> OPV;
   OPV.resize(*LF->NumberOfPointRecords);
   double P[2];
   size_t EZI = 0;
   for(uint32_t i = *LF->OffsetToPoints; i < LF->Data.size(); i += *LF->PointDataRecordLength){
      P[0] = CoordLongToDbl(*(int32_t*) &LF->Data[i], LF, 'x');
      P[1] = CoordLongToDbl(*(int32_t*) &LF->Data[i + sizeof(int32_t)], LF, 'y');
      OPV[EZI] = BWD(BWV, &P[0]);
      EZI++;
   }
   return(OPV);
}


/* Finds the closest(XY) sphagnum points to reference points inputs should be scale converted first
        inputs: *SPG = a pointer to the LasFile containing sphagnum surface points
                *RP = a pointer to the LasFile containing reference points
        output: A LasFile containing only the points from SPG that were closest to each reference point (XY distance).. copies SPG header
        
*/
SpClRefPtsO SphagClosestRefPts(LasFile *SPG, LasFile *RP){
   SpClRefPtsO GO;
   LasFile LasO;
   LasO.Data.resize(RP->Data.size());
   memcpy(&LasO.Data[0], &SPG->Data[0], *SPG->OffsetToPoints);
   LasFormat(&LasO);
   memcpy(&LasO.Data[107], &RP->Data[107], sizeof(uint32_t));
   double mindist, dist, closestz;
   std::vector<double> SCoords(3);
   std::vector<double> RCoords(3);
   std::vector<double> error(*RP->NumberOfPointRecords);
   std::vector<double> RefEle(*RP->NumberOfPointRecords);
   std::vector<double> ReConEle(*RP->NumberOfPointRecords);
   int32_t *rp, *sp, *spm;
   size_t MinI, ctr;
   MinI = -1;
   ctr = 0;
   for(uint32_t i = *RP->OffsetToPoints; i < RP->Data.size(); i += *RP->PointDataRecordLength){
      rp = (int32_t*) &RP->Data[i];
      mindist = 0;
      RCoords[0] = CoordLongToDbl(*rp,  RP, 'x');
      RCoords[1] = CoordLongToDbl(*(rp+1), RP,'y');
      RCoords[2] = CoordLongToDbl(*(rp+2), RP,'z');
      for(uint32_t j = *SPG->OffsetToPoints; j < SPG->Data.size(); j += *SPG->PointDataRecordLength){
         sp = (int32_t*) &SPG->Data[j];
         SCoords[0] = CoordLongToDbl(*sp, SPG, 'x');
         SCoords[1] = CoordLongToDbl(*(sp + 1), SPG, 'y');
         SCoords[2] = CoordLongToDbl(*(sp + 2), SPG, 'z');
         dist = DP2P(RCoords[0], RCoords[1], SCoords[0], SCoords[1]);
         if(mindist){
            if(dist < mindist){
               mindist = dist;
               MinI = j;
               spm = (int32_t*) &SPG->Data[j];
               closestz = SCoords[2]; 
            }
         } else {
            mindist = dist;
            MinI = j;
            spm = (int32_t*) &SPG->Data[j];
            closestz = SCoords[2];
         }
      }
      memcpy(&LasO.Data[i], &SPG->Data[MinI],*SPG->PointDataRecordLength);
      error[ctr] = closestz - (RCoords[2] - 0.5);
      RefEle[ctr] = RCoords[2] - 0.5;
      ReConEle[ctr] = closestz;
      ctr++; 
   }
   GO.RefEle = RefEle;
   GO.ReConEle = ReConEle; 
   GO.LasF = LasO;
   GO.Error = error;
   return(GO);
}




/* Performs the accuracy assessment for sphagnum points for spring 2017.. Takes a surface file, full plot point cloud, and reference file

*/
void SphagAcc(LasFile * SPG, LasFile * RF, LasFile * PC, std::string ACCFile, std::string BWFile, int32_t plot, int32_t year, std::string Dfile, double rad, double VS){
   std::string VoxSize = "0.01";
   double VoxS = 0.01;
//   LasScaleConv(PC, SPG, true, true);
//   LasScaleConv(PC, RF, true, true);
   SpClRefPtsO ACC = SphagClosestRefPts(SPG, RF);
   LasFormat(&ACC.LasF);
   std::ofstream accfile;
   accfile.open(ACCFile, std::ios::out | std::ios::binary);
   accfile.write((char*)&ACC.LasF.Data[0], ACC.LasF.Data.size());
//   LasScaleConv(PC, &ACC.LasF, true, false);
   std::vector<double> BWV = ReadTxt(BWFile);
   std::vector<uint32_t> GC = GridCounts(PC, VoxS);
   std::vector<double> BWDs = LasBWD(&ACC.LasF, &BWV[0]);
   int32_t * CP;
   std::vector<double> RPCoords(((RF->Data.size() - *RF->OffsetToPoints) / *RF->PointDataRecordLength)*2);
   std::vector<double> tmpV;
   size_t ctr = 0;
   for(size_t i = *RF->OffsetToPoints; i < RF->Data.size(); i += *RF->PointDataRecordLength){
      CP = (int32_t*) &RF->Data[i];
      tmpV = PointLongToDbl(*CP, *(CP+1), *(CP+2), RF);
      RPCoords[ctr] = tmpV[0];
      RPCoords[ctr + 1] = tmpV[1];
      ctr += 2; 
   }  
   LasFile ShrubVox = VoxelizeLas2Las(PC, VS);
   std::vector<uint32_t> SVols = PointVoxVolRadius(&ShrubVox, &RPCoords, rad);
   std::string OPFile = Dfile;
   std::ofstream ofile;
   ofile.open(OPFile, std::ios_base::app | std::ios_base::out);
   if((SVols.size() == BWDs.size()) & (SVols.size() == ACC.Error.size())){
      for(uint32_t i = 0; i < ACC.Error.size(); i++){
         std::cout << year << " " << plot << " "<< BWDs[i] << " " << ACC.RefEle[i] << " " << ACC.ReConEle[i] << " " << SVols[i] << " " << ACC.Error[i] << std::endl;
         ofile << year << " " << plot << " " << BWDs[i] << " " << ACC.RefEle[i] << " " << ACC.ReConEle[i] << " " << SVols[i] << " " << ACC.Error[i] << std::endl; 
      }
   } else {
      std::cout << "!!!! ~~~ !!!! WARNING OUTPUT VECTORS DON'T MATCH IN LENGTH !!!! ~~! !!!!" << std::endl;
   }
   std::cout << std::endl;
   std::cout << std::endl;
   ofile.close();
}






/* Calculates the percent of a vector below a tolerance (m) above the specified quantile value.. used for calculating percent hollow in SPRUCE plots
   inputs:
      quant = the quantile of VIN to be used
      tol = the tolerance above the quantile value.. (m) for SPRUCE plots
      VIN = a pointer to the vector representing SPRUCE plot elevations
   outputs:
      OP = the percent of the plot classified as hollow given the quantile and tolerance
*/
double PctHolThreshQuantTol(double quant, double tol, std::vector<double> *VIN){
   double i, ind, thresh, OP;
   std::sort(VIN->begin(),VIN->end());
   ind = VIN->size() / (1/quant);
   thresh = VIN->at(ind) + tol; 
   i = 0;
   while(VIN->at(i) < thresh){
      i++;
   }
   OP = (i+1) / VIN->size(); 
   return(OP);
}



// calulates percent of plot classified as hollow based on a 40 percent of max hummock height.. The original definition used at SPRUCE.. min + (Max - Min)*.4
std::vector<std::vector<double> > PctHolPctMaxEle(std::vector<double> *Zin, bool sort){
   if(sort){
      std::sort(Zin->begin(), Zin->end());
   }
   std::vector<std::vector<double> > OP;
   double thresh = ((Zin->back() - Zin->front()))*0.40 + Zin->front();
   double pcthol;
   double ctr = 0;
   for(size_t i = 0; i < Zin->size(); i++){
   /*   if(Zin->at(i) > thresh){
         ctr = i;
         break;
      }*/
      if(Zin->at(i) < thresh){
         ctr++;
      }   
   }
   std::cout << ctr << " " << Zin->size() << std::endl;
   pcthol = ctr/Zin->size();
   std::cout << "diff = " << Zin->back() - Zin->front() << " " << (Zin->back() - Zin->front()) * 0.4 << std::endl;
   std::cout << "min = " << Zin->front() << " max = " << Zin->back() << " thresh = " << thresh << " pct hol = " << pcthol << std::endl;
   std::cout << std::endl;
   std::cout << std::endl;
   return(OP);
}


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
SPRUCE_HolO SPRUCE_Hol(LasFile *Lasf, int32_t Gsig, int32_t Lsig, double Gthresh, double Lthresh, double GS){
     // add padding and grid the point cloud.. extracting Z values
     LasGridPad(Lasf, GS, true, 10, 10, 10);
     std::vector<uint32_t> GCS = GridCounts(Lasf, GS);
     std::vector<double> Grid = Las2DZMax(Lasf, GS, false, GS);
     double impute = *Lasf->MinZ + (*Lasf->MaxZ - *Lasf->MinZ)/2;

     // put elevation grid into OpenCV Mat object
     Mat MT = GridToCVMat(&Grid, GCS[0], true, *Lasf->MinZ + (*Lasf->MaxZ - *Lasf->MinZ)/2 );
     Mat LGB;
     Mat GGB;

     // Make sure sigmas are odd numbers and perform smoothing

     // gaussian for gradient
     if(Gsig % 2 == 0){
          GaussianBlur(MT, GGB, Size(Gsig*7+1, Gsig*7+1), Gsig, Gsig);
     } else {
          GaussianBlur(MT, GGB, Size(Gsig*7, Gsig*7), Gsig, Gsig);
     }

     // gaussian for laplacian
     if(Lsig % 2 == 0){
          GaussianBlur(MT, LGB, Size(Lsig*7+1, Lsig*7+1), Lsig, Lsig);
     } else {
          GaussianBlur(MT, LGB, Size(Lsig*7, Lsig*7), Lsig, Lsig);
     }

     // create gradient and perform thresholding 
     std::vector<double> GV = CVImGradGrid(&GGB);
     Mat grad = GridToCVMat(&GV, GCS[0], false, 0.0);
     Mat GT;
     threshold(grad, GT, Gthresh, 255, 1);

     // create laplacian and perform thresholding
     Mat lap;
     Laplacian(LGB, lap, LGB.depth(), 3);
     Mat LT; 
     threshold(lap, LT, Lthresh, 255, 0);

     // Perform masking from 1st & 2nd derivatives
     std::vector<double> ltv = CVMatToGrid(&LT).Grid;
     std::vector<double> gtv = CVMatToGrid(&GT).Grid;
     std::vector<double> Class(ltv.size());

     for(size_t i = 0; i < Class.size(); i++){
          if(gtv[i] != 0 & ltv[i] != 0){
               Class[i] = 100;
          } 
          if(gtv[i] != 0 & Grid[i] == 0){
               Class[i] = 50;     
          }
     }
     
     // testing here...
     std::vector<double> LV = CVMatToGrid(&lap).Grid;
     std::vector<double> G2 = Grid;
     NormalizeGrid(&G2);
     NormalizeGrid(&LV);
     std::vector<double> test = GridScaleByGrid(&LV, &G2);
/*     Mat testM = GridToCVMat(&test, GCS[0], false, 0.0);
     Mat eq;

     std::cout << testM.at<double>(400,400) << " " <<  testM.at<double>(300,300) << " " << testM.at<double>(500,500)  << " " << testM.at<double>(200,500) << std::endl;
     imshow("Original", testM);
     waitKey(0);
     testM.convertTo(testM, CV_8UC1);
     std::cout << testM.at<double>(400,400) << " " <<  testM.at<double>(300,300) << " " << testM.at<double>(500,500)  << " " << testM.at<double>(200,500) << std::endl;
    imshow("converted", testM);
     waitKey(0);
     equalizeHist(testM, eq);
     imshow("equal", eq);
     waitKey(0);
     eq.convertTo(eq, CV_64F);  
     imshow("conv back", eq);
     waitKey(0); */
     // to here...


     // create objects for output
     SPRUCE_HolO OP;
     OP.DL = GCS[0];
     OP.Gsig = Gsig;
     OP.Lsig = Lsig;
     OP.Gthresh = Gthresh;
     OP.Lthresh = Lthresh;
     OP.GridSize = GS;
     OP.Raw = Grid;
     OP.GGB = CVMatToGrid(&GGB).Grid;
     OP.LGB = CVMatToGrid(&LGB).Grid;
     OP.Grad = GV;
     OP.Lap = CVMatToGrid(&lap).Grid;
// this is correct... just testing     OP.GradT = CVMatToGrid(&GT).Grid;
     OP.GradT = test;
     OP.LapT = CVMatToGrid(&LT).Grid;
     OP.HumHol = Class; 

     return(OP);
}



/* Calculates the percent of a SPRUCE chamber that is classified as hollow based on an elevation threshold that is the minimum chamber elevation plus the range of chamber elevation times the threshold.. (e.g., for a threshold of 0.4 (40% max hummock height) hollow would be points with elevations less than:
min(elevation) + (max(elevation) - min(elevatioN))*thresh 
     inputs:
          VIN = a pointer to the grid representing SPRUCE chamber elevations
          thresh = the percent of the elevation range to use as the threshold (needs to be between 0 and 1)
     outputs:
          OP = vector of the classified grid.. 100 = hummock, 0 = hollow
*/
std::vector<double> HollowByPercentEleRange(std::vector<double> *VIN, double thresh){
     std::vector<double> OP(VIN->size());
     double Min = *min_element(std::begin(*VIN), std::end(*VIN));
     double Max = *max_element(std::begin(*VIN), std::end(*VIN));
     double Range = Max - Min;
     double EleThresh = Min + (Range*thresh);
     for(size_t i = 0; i < OP.size(); i++){
          if(VIN->at(i) >= EleThresh){
               OP[i] = 100;
          } else {
               OP[i] = 0;    
          }   
     }
     return(OP);
}


struct HolInd{
     std::vector<double> LV;  // Original laplacian values before weighting
     std::vector<double> GV;  // Original gradient values before weighting
     std::vector<double> HI;  // Calculated "Hollow Index".. in grid format
     std::vector<double> EC;  // Elevation after values sigmoid weighting.. in grid format
     std::vector<double> LC;  // Laplacian values.. in grid format
     std::vector<double> GC;  // Gradient values.. in grid format
     std::vector<double> ECR; // Vector of elevation values.. in vector format with all imputed cells removed.. for calculating summary statistics
     std::vector<double> LCR; // Vector of gradient values.. in vector format with all imputed cells removed.. for calculating summary statistics
     std::vector<double> GCR; // Vector of laplacian values.. in vector format with all imputed cells removed.. for calculating summary statistics
     std::vector<double> ECC; // Vector of weights from the sigmoid function for elevation on an equally spaced vector from min to max 
     std::vector<double> LCC; // Vector of weights form the sigmoid function for laplacian on an equally spaced vector from min to max
     std::vector<double> GCC; // Vector of wieghts from the sigmoid function for gradient on an equally spaced vector from min to max
     std::vector<double> EVals; // Vector of elevation values that weights in ECC are associated with
     std::vector<double> LVals; // Vector of laplacian values that weights in LCC are associated with
     std::vector<double> GVals; // Vector of gradient values that weights in GCC are associated with
};



/* Calculates the "Hollow Index" based on elevations, gradient, and laplacian.. multiplicative version 
     inputs:
          IP = pointer to the GridDerivatives object to be operated on
          ignore = a value that is ignored, typically the value that is imputed for areas outside SPRUCE plots
     outputs:
          OP = HolInd structure containing "Hollow Index" grid and relevant precursors  
*/
HolInd HollowIndex(GridDerivatives *IP, double ignore, double EleT){
     HolInd OP;
     OP.HI.resize(IP->Raw.size());
     OP.LV.resize(IP->Raw.size());
     OP.GV.resize(IP->Raw.size());
     uint16_t breaks = 2000; 
     for(size_t i = 0; i < IP->Raw.size(); i++){
          OP.LV[i] = IP->Laplacian[i];
          OP.GV[i] = IP->Gradient[i]; 
          if(IP->Raw[i] != ignore){
               OP.ECR.push_back(IP->Raw[i]);
               OP.GCR.push_back(IP->Gradient[i]);
               OP.LCR.push_back(IP->Laplacian[i]);
          }
     }
     std::sort(OP.ECR.begin(), OP.ECR.end());
     std::sort(OP.LCR.begin(), OP.LCR.end());
     std::sort(OP.GCR.begin(), OP.GCR.end());   

     double EStep = ( OP.ECR.back() - OP.ECR.front() ) / breaks;
     double LStep = ( OP.LCR.back() - OP.LCR.front() ) / breaks;
     double GStep = ( OP.GCR.back() - OP.GCR.front() ) / breaks;

     OP.EVals.resize(breaks); OP.LVals.resize(breaks); OP.GVals.resize(breaks);
     OP.EVals[0] = OP.ECR.front(); OP.LVals[0] = OP.LCR.front(); OP.GVals[0] = OP.GCR.front();

     for(size_t i = 1; i < breaks; i++){
          OP.EVals[i] = OP.EVals[i-1] + EStep;
          OP.LVals[i] = OP.LVals[i-1] + LStep;
          OP.GVals[i] = OP.GVals[i-1] + GStep;
     }

     double Zmin = *std::min_element(std::begin(OP.ECR), std::end(OP.ECR));
     double Zmax = *std::max_element(std::begin(OP.ECR), std::end(OP.ECR));
     double Zmid = Zmin + ((Zmax - Zmin) / 2);
     double Zmean = VecMean(&OP.ECR);
     double ZSD = VecSD(&OP.ECR);

     double Gmin = *std::min_element(std::begin(OP.GCR), std::end(OP.GCR)); 
     double Gmax = *std::max_element(std::begin(OP.GCR), std::end(OP.GCR)); 
     double Gmid = Gmin + ((Gmax - Gmin) / 2);
     double GSD = VecSD(&OP.GCR);

     double Lmin = *std::min_element(std::begin(OP.LCR), std::end(OP.LCR));
     double Lmax = *std::max_element(std::begin(OP.LCR), std::end(OP.LCR));
     double Lmid = Lmin + ((Lmax - Lmin) / 2);
     double LSD = VecSD(&OP.LCR);

     OP.EC = LogisticWeights(&IP->Raw, Zmean, 2, 1 / ZSD, EleT, true); 
     OP.LC = LogisticWeights(&IP->Laplacian, 0.0, 2, 1 / LSD, 0.0, false);     
     OP.GC = LogisticWeights(&IP->Gradient, 45, 1, 0.0416, 1.5, true);

     OP.ECC = LogisticWeights(&OP.EVals, Zmean, 2, 1 / ZSD , EleT, true); 
     OP.LCC = LogisticWeights(&OP.LVals, 0.0, 2, 1 / LSD, 0.0, false);     
     OP.GCC = LogisticWeights(&OP.GVals, 45, 1, 0.0416, 1.5, true);
     for(size_t i = 0; i < OP.HI.size(); i++){
          if(IP->Raw[i] != 0){
                OP.HI[i] = (OP.EC[i]) * (OP.LC[i]) * (OP.GC[i]);
          }
     }
     return(OP);
}



/* Calculates the "Hummock Index" based on elevations, gradient, and laplacian.. multiplicative version 
     inputs:
          IP = pointer to the GridDerivatives object to be operated on
          ignore = a value that is ignored, typically the value that is imputed for areas outside SPRUCE plots
     outputs:
          OP = HolInd structure containing "Hummock Index" grid and relevant precursors (spelling) 
*/
HolInd HummockIndex(GridDerivatives *IP, double ignore){
     HolInd OP;
     OP.HI.resize(IP->Raw.size());
     OP.LV.resize(IP->Raw.size());
     OP.GV.resize(IP->Raw.size());
     uint16_t breaks = 2000; 
     for(size_t i = 0; i < IP->Raw.size(); i++){
          OP.LV[i] = IP->Laplacian[i];
          OP.GV[i] = IP->Gradient[i]; 
          if(IP->Raw[i] != ignore){
               OP.ECR.push_back(IP->Raw[i]);
               OP.GCR.push_back(IP->Gradient[i]);
               OP.LCR.push_back(IP->Laplacian[i]);
          }
     }
     std::sort(OP.ECR.begin(), OP.ECR.end());
     std::sort(OP.LCR.begin(), OP.LCR.end());
     std::sort(OP.GCR.begin(), OP.GCR.end());   

     double EStep = ( OP.ECR.back() - OP.ECR.front() ) / breaks;
     double LStep = ( OP.LCR.back() - OP.LCR.front() ) / breaks;
     double GStep = ( OP.GCR.back() - OP.GCR.front() ) / breaks;

     OP.EVals.resize(breaks); OP.LVals.resize(breaks); OP.GVals.resize(breaks);
     OP.EVals[0] = OP.ECR.front(); OP.LVals[0] = OP.LCR.front(); OP.GVals[0] = OP.GCR.front();

     for(size_t i = 1; i < breaks; i++){
          OP.EVals[i] = OP.EVals[i-1] + EStep;
          OP.LVals[i] = OP.LVals[i-1] + LStep;
          OP.GVals[i] = OP.GVals[i-1] + GStep;
     }

     double Zmin = *std::min_element(std::begin(OP.ECR), std::end(OP.ECR));
     double Zmax = *std::max_element(std::begin(OP.ECR), std::end(OP.ECR));
     double Zmid = Zmin + ((Zmax - Zmin) / 2);
     double Zmean = VecMean(&OP.ECR);
     double ZSD = VecSD(&OP.ECR);

     double Gmin = *std::min_element(std::begin(OP.GCR), std::end(OP.GCR)); 
     double Gmax = *std::max_element(std::begin(OP.GCR), std::end(OP.GCR)); 
     double Gmid = Gmin + ((Gmax - Gmin) / 2);
     double GSD = VecSD(&OP.GCR);

     double Lmin = *std::min_element(std::begin(OP.LCR), std::end(OP.LCR));
     double Lmax = *std::max_element(std::begin(OP.LCR), std::end(OP.LCR));
     double Lmid = Lmin + ((Lmax - Lmin) / 2);
     double LSD = VecSD(&OP.LCR);

     OP.EC = LogisticWeights(&IP->Raw, Zmean, 2, 1 / ZSD, -1.0, false); 
     OP.LC = LogisticWeights(&IP->Laplacian, 0.0, 2, 1 / LSD, 2.0, true);     
     OP.GC = LogisticWeights(&IP->Gradient, 45, 1, 0.0416, 1.5, true);

     OP.ECC = LogisticWeights(&OP.EVals, Zmean, 2, 1 / ZSD , -1.0, false);
     OP.LCC = LogisticWeights(&OP.LVals, 0.0, 2, 1 / LSD, 2.0, true);     
     OP.GCC = LogisticWeights(&OP.GVals, 45, 1, 0.0416, 1.5, true);
     for(size_t i = 0; i < OP.HI.size(); i++){
          if(IP->Raw[i] != 0){
                OP.HI[i] = (OP.EC[i]) * (OP.LC[i]) * (OP.GC[i]); 
          }
     }
     return(OP);
}





/* Calculates the "Hollow Index" based on elevations, gradient, and laplacian.. linear combinatin version
     inputs:
          IP = pointer to the GridDerivatives object to be operated on
          ignore = a value that is ignored, typically the value that is imputed for areas outside SPRUCE plots
     outputs:
          OP = HolInd structure containing "Hollow Index" grid and relevant precursors (spelling) 


*/
HolInd HollowIndexLC(GridDerivatives *IP, double EQuant, double LScale, double GScale, bool LocalLap, bool LocalGrad){
     HolInd OP;
     OP.HI.resize(IP->Raw.size());
     uint16_t breaks = 2000; 
     for(size_t i = 0; i < IP->Raw.size(); i++){
          if(IP->Raw[i] != LScale){
               OP.ECR.push_back(IP->Raw[i]);
               OP.GCR.push_back(IP->Gradient[i]);
               OP.LCR.push_back(IP->Laplacian[i]);
          }
     }
     std::sort(OP.ECR.begin(), OP.ECR.end());
     std::sort(OP.LCR.begin(), OP.LCR.end());
     std::sort(OP.GCR.begin(), OP.GCR.end());   

     double EStep = ( OP.ECR.back() - OP.ECR.front() ) / breaks;
     double LStep = ( OP.LCR.back() - OP.LCR.front() ) / breaks;
     double GStep = ( OP.GCR.back() - OP.GCR.front() ) / breaks;

     OP.EVals.resize(breaks); OP.LVals.resize(breaks); OP.GVals.resize(breaks);
     OP.EVals[0] = OP.ECR.front(); OP.LVals[0] = OP.LCR.front(); OP.GVals[0] = OP.GCR.front();

     for(size_t i = 1; i < breaks; i++){
          OP.EVals[i] = OP.EVals[i-1] + EStep;
          OP.LVals[i] = OP.LVals[i-1] + LStep;
          OP.GVals[i] = OP.GVals[i-1] + GStep;
     }

     double Zmin = *std::min_element(std::begin(OP.ECR), std::end(OP.ECR));
     double Zmax = *std::max_element(std::begin(OP.ECR), std::end(OP.ECR));
     double Zmid = Zmin + ((Zmax - Zmin) / 2);
     double Zmean = VecMean(&OP.ECR);
     double ZSD = VecSD(&OP.ECR);

     double Gmin = *std::min_element(std::begin(OP.GCR), std::end(OP.GCR)); 
     double Gmax = *std::max_element(std::begin(OP.GCR), std::end(OP.GCR)); 
     double Gmid = Gmin + ((Gmax - Gmin) / 2);
     double GSD = VecSD(&OP.GCR);

     double Lmin = *std::min_element(std::begin(OP.LCR), std::end(OP.LCR));
     double Lmax = *std::max_element(std::begin(OP.LCR), std::end(OP.LCR));
     double Lmid = Lmin + ((Lmax - Lmin) / 2);
     double LSD = VecSD(&OP.LCR);

     OP.EC = LogisticWeights(&IP->Raw, Zmean, 2, 1 / 2*ZSD, 1.0, true); 
     OP.LC = LogisticWeights(&IP->Laplacian, 0.0, 2, 1 / 2*LSD, 0.0, false);     
     OP.GC = LogisticWeights(&IP->Gradient, 45, 1, 0.0416, 1, true);

     OP.ECC = LogisticWeights(&OP.EVals, Zmean, 2, 1 / 2*ZSD , 1.0, true); 
     OP.LCC = LogisticWeights(&OP.LVals, 0.0, 2, 1 / 2*LSD, 0.0, false);     
     OP.GCC = LogisticWeights(&OP.GVals, 45, 1, 0.0416, 1, true);
     for(size_t i = 0; i < OP.HI.size(); i++){
          if(IP->Raw[i] != 0){
                OP.HI[i] = (OP.EC[i]) + (OP.LC[i]) + (OP.GC[i]);
          }
     }
     return(OP);
}


/* makes a point cloud where xyz are elevation gradient and laplacian and GPS time stores hollow index for coloring the point cloud

///
LasFile HollowIndexParamPC(std::vector<double> *EV, std::vector<double> *LV, std::vector<double> *GV, uint32_t breaks, LasFile * LasTemplate){
     double Emin = *std::min_element(std::begin(*EV), std::end(*EV));   
     double Emax = *std::max_element(std::begin(*EV), std::end(*EV));
     double Lmin = *std::min_element(std::begin(*LV), std::end(*LV));
     double Lmax = *std::max_element(std::begin(*LV), std::end(*LV));
     double Gmin = *std::min_element(std::begin(*GV), std::end(*GV));
     double Gmax = *std::max_element(std::begin(*GV), std::end(*GV));
     double Estep = (Emax - Emin) / breaks;
     double Lstep = (Lmax - Lmin) / breaks;
     double Gstep = (Gmax - Gmin) / breaks;
     std::vector<double> ES(breaks);
     std::vector<double> LS(breaks);
     std::vector<double> GS(breaks);
     for(size_t i = 0; i < breaks; i++){
         ES[i] = Emin + i*Estep;
         LS[i] = Lmin + i*Lstep;
         GS[i] = Gmin + i*Gstep;
     }

     std::vector<double> EW = LogisticWeights(&ES, VecMean(EV), 2, 1.5 / VecSD(EV), 2.0, true);
     std::vector<double> LW = LogisticWeights(&LS, 0.0, 2, 1.732 / VecSD(LV), 0.0, false);
     std::vector<double> GW = LogisticWeights(&GS, 25, 1, 0.0516, 1.5, true);

     LasFile OP = *LasTemplate;
     OP.Data.resize(*OP.OffsetToPoints + (pow(breaks,3) * (*OP.PointDataRecordLength)));


        
}

*/


/* takes a LasFile representing the bog surface of a SPRUCE plot and a time series of water table depths and calculates the distance to water table for each timestep



*/
