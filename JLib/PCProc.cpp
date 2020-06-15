#include<vector> 
#include<string>
#include<cmath>
#include<stdint.h>
#include"Misc.h"
#include"LasReadWrite.h"
#include<iterator>
#include<PCProc.h>
#include<iostream>

// Converts coordinates in double (UTM) to the scaled int32_t from a LasFile
int64_t CoordDblToLong(double Pin, LasFile *Lasf, char Coord){
   int64_t Pout = 0;
   if(Coord == 'x' || Coord == 'X'){
     Pout = (Pin - *Lasf->XOffset) / *Lasf->XScaleFactor;
   } else if(Coord == 'y' || Coord == 'Y'){
     Pout = (Pin - *Lasf->YOffset) / *Lasf->YScaleFactor;
   } else if(Coord == 'z' || Coord == 'Z'){
     Pout = (Pin - *Lasf->ZOffset) / *Lasf->ZScaleFactor;
   } else {
     std::cout << "!!! INVALID COORDINATE !!!" << std::endl;
   }
   return(Pout);
}


// Converts all point coordinates in double to (UTM) to scaled int32_t from a LasFile
std::vector<int32_t> PointDblToLong(double X, double Y, double Z, LasFile *Lasf){
   std::vector<int32_t> Pout(3);
   Pout[0] = (X - *Lasf->XOffset) / *Lasf->XScaleFactor;
   Pout[1] = (Y - *Lasf->YOffset) / *Lasf->YScaleFactor;
   Pout[2] = (Z - *Lasf->ZOffset) / *Lasf->ZScaleFactor;
   return(Pout);   
}

// Converts coordinates in int32_ts from LasFile to double (UTM) using the files scalings & offsets
double CoordLongToDbl(int32_t Pin, LasFile *Lasf, char Coord){
   double Pout = 0;
   if(Coord == 'x' || Coord == 'X'){
      Pout = (Pin * (*Lasf->XScaleFactor) + *Lasf->XOffset);
   } else if(Coord == 'y' || Coord == 'Y'){
      Pout = (Pin * (*Lasf->YScaleFactor) + *Lasf->YOffset);
   } else if(Coord == 'z' || Coord == 'Z'){
      Pout = (Pin * (*Lasf->ZScaleFactor) + *Lasf->ZOffset);
   } else {
      std::cout << "!!! INVALID COORDINATE !!!" << std::endl;
   }
   return(Pout);
}

// Converts all point coordinates in scaled int32_t from a LasFile to double (UTM)
std::vector<double> PointLongToDbl(int32_t X, int32_t Y, int32_t Z, LasFile *Lasf){
   std::vector<double> Pout(3);
   Pout[0] = (X * (*Lasf->XScaleFactor) + *Lasf->XOffset);
   Pout[1] = (Y * (*Lasf->YScaleFactor) + *Lasf->YOffset);
   Pout[2] = (Z * (*Lasf->ZScaleFactor) + *Lasf->ZOffset);
   return(Pout);
}


/* Given a vector of voxel IDs  .. output 3x the length of input orderd x1,y1,z1,x2,y2,z2..etc
   input:  voxels = a pointer to a vector of Voxel IDs to be converted to coordinates
       XGC = # of x grid cells used in voxelization
       YGC = # of y grid cells used in voxeliation
       VoxS = the size of voxel used in voxeliation
       Lasf = a pointer to the LasFile from which voxelization occured
   output: voxr = a vector of double containing xyz coordinates for each occupied voxel in format x1,y1,z1,x2,y2,z2..etc
*/
std::vector<double> VoxCoord(std::vector<uint32_t> *voxels, int32_t XGC, int32_t YGC, double VoxS, LasFile *LasF) {
   int32_t zm, ym, xm;
   double zp, yp, xp, Xcellsize, Ycellsize, CellH;
   Xcellsize = VoxS; Ycellsize = VoxS; CellH = VoxS;
   std::vector<double> voxr(voxels->size()*3);
   for (uint32_t i = 0; i < voxels->size(); i++) {
      zm = floor(voxels->at(i) / (XGC*YGC));
      ym = floor((voxels->at(i) - (zm*XGC*YGC)) / XGC);
      xm = voxels->at(i) - (zm*XGC*YGC + ym*XGC);
      xp = *LasF->MinX + ((xm * Xcellsize) + (Xcellsize / 2));  // centered x value of voxel
      yp = *LasF->MinY + ((ym * Ycellsize) + (Ycellsize / 2)); // centered y value of voxel
      zp = *LasF->MinZ + ((zm * CellH) + (CellH / 2));      // centered z value of voxel
      voxr[i*3] = xp; 
      voxr[i*3 + 1] = yp;
      voxr[i*3 + 2] = zp;
   }
   return(voxr);
}


/*   converts min/max double values from header to a scaled vector of int32_ts, so it is comparable to point records w/out scaling
   order is: 0:MinX,1:MaxX,2:MinY,3:MaxY,4:MinZ,5:MaxZ                                     */
std::vector<int32_t> ScaleMinMax(LasFile *LasIn){
   std::vector<int32_t> OPV(6);
   OPV[0] = (*LasIn->MinX - *LasIn->XOffset) / *LasIn->XScaleFactor;
   OPV[1] = (*LasIn->MaxX - *LasIn->XOffset) / *LasIn->XScaleFactor;
   OPV[2] = (*LasIn->MinY - *LasIn->YOffset) / *LasIn->YScaleFactor;
   OPV[3] = (*LasIn->MaxY - *LasIn->YOffset) / *LasIn->YScaleFactor;
   OPV[4] = (*LasIn->MinZ - *LasIn->ZOffset) / *LasIn->ZScaleFactor;
   OPV[5] = (*LasIn->MaxZ - *LasIn->ZOffset) / *LasIn->ZScaleFactor;
   return(OPV);
}



// !!!!!!!!!!!!!!! this does not work properly for point clouds saved with low resolution (i.e., x-y-z scalefactors do not retain several decimal places) need to fix by converting to doubles before
// takes a las file and voxel size and returns a vector with numbers of cells in each coordinate order: 0 = X, 1 = Y, 2 = Z
std::vector<uint32_t> GridCounts(LasFile *Lasf, float VSize){
   int32_t Xcellsize, Ycellsize, Zcellsize;
   Xcellsize = VSize / *Lasf->XScaleFactor; Ycellsize = VSize /  *Lasf->YScaleFactor; Zcellsize = VSize /  *Lasf->ZScaleFactor;
   std::vector<int32_t> MinMax = ScaleMinMax(Lasf);
   std::vector<uint32_t> GC(3);
   GC[0] = ceil(((MinMax[1] + Xcellsize) - MinMax[0]) / Xcellsize);
   GC[1] = ceil(((MinMax[3] + Ycellsize) - MinMax[2]) / Ycellsize);
   GC[2] = ceil(((MinMax[5] + Zcellsize) - MinMax[4]) / Zcellsize);
   return(GC);
}


/* takes a pointer to a LasFile and returns a vector<bool> of voxel presence/absence
   inputs: Lasf = a pointer to the LasFile to be voxelized
       Vsize = the voxel size (meters)
   output: Bout = a 1D projected vector of bools; true = occupied & false = !occupied; 1D projection hierachry = z > y > x 

*/
VoxBO Voxels_bool(LasFile *Lasf, float VSize){
   std::vector<uint32_t> GCS = GridCounts(Lasf,VSize);
   std::vector<bool> Bout(GCS[0]*GCS[1]*GCS[2],false);
   int32_t xpos, ypos, zpos, Xcellsize, Ycellsize, Zcellsize;
   Xcellsize = VSize / *Lasf->XScaleFactor; Ycellsize = VSize /  *Lasf->YScaleFactor; Zcellsize = VSize /  *Lasf->ZScaleFactor;
   std::vector<int32_t> MinMax = ScaleMinMax(Lasf);
   int32_t* CP;
   int64_t VC = 0;
   for(uint32_t i = *Lasf->OffsetToPoints;i <  Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
       CP = (int32_t*) &Lasf->Data[i];
       xpos = gridpos(MinMax[0], Xcellsize, *CP);
       ypos = gridpos(MinMax[2], Ycellsize, *(CP+1));
       zpos = gridpos(MinMax[4], Zcellsize, *(CP+2));
      if(Bout[vxl(xpos,ypos,zpos,GCS[0],GCS[1])] != true){
         VC++;   
      } 
       Bout[vxl(xpos,ypos,zpos,GCS[0],GCS[1])] = true;
   }
//   std::cout << VC << std::endl;
   VoxBO OP;
   OP.BV = Bout;
   OP.CVV = VC;
   return(OP);
} 


/* takes a pointer to a LasFile and returns a vector<bool> of voxel presence/absence on all points above or below a threshold
   inputs: Lasf = a pointer to the LasFile to be voxelized
       Vsize = the voxel size (meters)
      Thresh = elevation threshold to be used (meters)
      UpDown = bool that tells whether to voxelize points above or below the threshold: true = voxelize points above thresh,
         false = voxelize points below thresh
   output: Bout = a 1D projected vector of bools; true = occupied & false = !occupied; 1D projection hierachry = z > y > x 

*/
VoxBO Voxels_bool_thresh(LasFile *Lasf, float VSize, double Thresh, bool UpDown){
   std::vector<uint32_t> GCS = GridCounts(Lasf,VSize);
   std::vector<bool> Bout(GCS[0]*GCS[1]*GCS[2],false);
   int32_t xpos, ypos, zpos, Xcellsize, Ycellsize, Zcellsize;
   Xcellsize = VSize / *Lasf->XScaleFactor; Ycellsize = VSize /  *Lasf->YScaleFactor; Zcellsize = VSize /  *Lasf->ZScaleFactor;
   std::vector<int32_t> MinMax = ScaleMinMax(Lasf);
   int32_t* CP;
   int64_t VC = 0;
   if(UpDown){
      for(uint32_t i = *Lasf->OffsetToPoints;i <  Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
          CP = (int32_t*) &Lasf->Data[i];
          if(*(CP+2) > CoordDblToLong(Thresh, Lasf, 'z')){
            xpos = gridpos(MinMax[0], Xcellsize, *CP);
               ypos = gridpos(MinMax[2], Ycellsize, *(CP+1));
             zpos = gridpos(MinMax[4], Zcellsize, *(CP+2));
            if(Bout[vxl(xpos,ypos,zpos,GCS[0],GCS[1])] != true){
               VC++;   
            }
             Bout[vxl(xpos,ypos,zpos,GCS[0],GCS[1])] = true;
         }
      }
   } else {
      for(uint32_t i = *Lasf->OffsetToPoints;i <  Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
          CP = (int32_t*) &Lasf->Data[i];
          if(*(CP+2) < CoordDblToLong(Thresh, Lasf, 'z')){
            xpos = gridpos(MinMax[0], Xcellsize, *CP);
               ypos = gridpos(MinMax[2], Ycellsize, *(CP+1));
               zpos = gridpos(MinMax[4], Zcellsize, *(CP+2));
            if(Bout[vxl(xpos,ypos,zpos,GCS[0],GCS[1])] != true){
               VC++;   
            }
               Bout[vxl(xpos,ypos,zpos,GCS[0],GCS[1])] = true;
         }
      }
   }
   VoxBO OP;
   OP.BV = Bout;
   OP.CVV = VC;
   return(OP);
}



/* takes a pointer to a LasFile and returns a vector voxels (by ID #) that are occupied
   inputs: Lasf = pointer to the LasFile to be voxelized
       VSize = voxel size (m)
       sort = whether the output vector that should sorted so the "lowest" numbered voxel is the first value in the vector AND reduced so there are no replicates
       VV = sorted vector of voxel IDs 
       OcVox = unsorted vector of voxel IDs
*/
std::vector<uint32_t> Voxels_ID(LasFile* Lasf, float VSize, bool sort){
   std::vector<uint32_t> GCS = GridCounts(Lasf,VSize);
std::cout << GCS[0] << " " << GCS[1] << " " << GCS[2] << std::endl;
std::cout << VSize << std::endl;
   std::vector<uint32_t> OcVox(*Lasf->NumberOfPointRecords);
   int32_t xpos, ypos, zpos, Xcellsize, Ycellsize, Zcellsize;
      Xcellsize = VSize / *Lasf->XScaleFactor; Ycellsize = VSize /  *Lasf->YScaleFactor; Zcellsize = VSize /  *Lasf->ZScaleFactor;
   std::vector<int32_t> MinMax = ScaleMinMax(Lasf);
   int32_t* CP;
   //std::cout << *Lasf->OffsetToPoints << " " << Lasf->Data.size() << " " << *Lasf->PointDataRecordLength << std::endl;
   for(uint32_t i = *Lasf->OffsetToPoints; i <  Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
       CP = (int32_t*) &Lasf->Data[i];

       xpos = gridpos(MinMax[0], Xcellsize, *CP);
       ypos = gridpos(MinMax[2], Ycellsize, *(CP+1));
       zpos = gridpos(MinMax[4], Zcellsize, *(CP+2));
      OcVox[(i - *Lasf->OffsetToPoints) / *Lasf->PointDataRecordLength] = vxl(xpos,ypos,zpos,GCS[0],GCS[1]);
   }
   if(sort){
       std::sort(OcVox.begin(),OcVox.end());   
       std::vector<uint32_t> VV;
       VV.push_back(OcVox[0]);
       for(uint32_t i = 0; i < OcVox.size(); i++){
           if(OcVox[i]!= VV.back()){
                 VV.push_back(OcVox[i]);
           }
       }
      return(VV);
   } else {
   return(OcVox);
   }
}


/* takes a pointer to a LasFile and returns a vector of sorted and reduced voxels (by ID #) that are occupied
   inputs: Lasf = pointer to the LasFile to be voxelized
       VSize = voxel size (m)
       sort = whether the output vector that should sorted so the "lowest" numbered voxel is the first value in the vector
      Thresh = elevation threshold to be used (meters)
      UpDown = bool that tells whether to voxelize points above or below the threshold: true = voxelize points above thresh,
         false = voxelize points below thresh
   outputs:     
      VV = sorted vector of voxel IDs 
       OcVox = unsorted vector of voxel IDs
*/
std::vector<uint32_t> Voxels_ID_Thresh(LasFile* Lasf, float VSize, bool sort, bool UpDown, double Thresh){
   std::vector<uint32_t> GCS = GridCounts(Lasf,VSize);
   std::vector<uint32_t> OcVox(*Lasf->NumberOfPointRecords);
   int32_t xpos, ypos, zpos, Xcellsize, Ycellsize, Zcellsize;
   Xcellsize = VSize / *Lasf->XScaleFactor; Ycellsize = VSize /  *Lasf->YScaleFactor; Zcellsize = VSize /  *Lasf->ZScaleFactor;
   std::vector<int32_t> MinMax = ScaleMinMax(Lasf);
   int32_t* CP;
   if(Thresh){
      for(uint32_t i = *Lasf->OffsetToPoints;i <  Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
          CP = (int32_t*) &Lasf->Data[i];
         if(*(CP+2) > CoordDblToLong(Thresh, Lasf, 'z')){
             xpos = gridpos(MinMax[0], Xcellsize, *CP);
               ypos = gridpos(MinMax[2], Ycellsize, *(CP+1));
               zpos = gridpos(MinMax[4], Zcellsize, *(CP+2));
                 OcVox[(i - *Lasf->OffsetToPoints) / *Lasf->PointDataRecordLength] = vxl(xpos,ypos,zpos,GCS[0],GCS[1]);
         }
      }
   } else {
      for(uint32_t i = *Lasf->OffsetToPoints;i <  Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
          CP = (int32_t*) &Lasf->Data[i];
         if(*(CP+2) < CoordDblToLong(Thresh, Lasf, 'z')){
             xpos = gridpos(MinMax[0], Xcellsize, *CP);
               ypos = gridpos(MinMax[2], Ycellsize, *(CP+1));
               zpos = gridpos(MinMax[4], Zcellsize, *(CP+2));
                 OcVox[(i - *Lasf->OffsetToPoints) / *Lasf->PointDataRecordLength] = vxl(xpos,ypos,zpos,GCS[0],GCS[1]);
         }
      }
   }
   if(sort){
       std::sort(OcVox.begin(),OcVox.end());   
       std::vector<uint32_t> VV;
       VV.push_back(OcVox[0]);
      for(uint32_t i = 0; i < OcVox.size(); i++){
           if(OcVox[i]!= VV.back()){
            VV.push_back(OcVox[i]);
           }
       }
      std::sort(VV.begin(), VV.end());
      std::cout << OcVox.size() << " " << VV.size() << std::endl;
      for(size_t i = 0; i < VV.size() - 1; i++){
         if(VV[i] == VV[i+1]){
            std::cout << i << " " << VV[i] << " "  << i + 1 << " " << VV[i + 1] << std::endl;
         }
      }
      return(VV);
   } else {
   return(OcVox);
   }
}


/* converts second las file to have same scaling factors and offsets as first, if Min = TRUE the min xyzs in the header are set to match,
this is so they can be voxelized and voxel numbers will match for comparison; if Max = TRUE the max xyzs in the header are set to match;
by setting Min = true & Max = False you save a large amount of memory when using SVV_bool   
   inputs: ConFrom = LasFile to convert ConTo's scaling to 
       ConTo = LasFile to convert to scaling from ConFrom
       Min =   true = set minimums of bounding boxes equal.. useful for when voxelizing multiple LasFiles and comparing voxels 
       Max =   true = set maximums of bounding boxes equal 
   output: N/A.. the ConTo file will now have the same scaling and mins/maxes (if true) as ConFrom  
*/
void LasScaleConv(LasFile *ConFrom, LasFile *ConTo, bool Min, bool Max){
   int32_t* IP;
   int32_t cords[3];
   for (uint32_t i = *ConTo->OffsetToPoints; i < ConTo->Data.size(); i+= *ConTo->PointDataRecordLength){
       IP = (int32_t*) &ConTo->Data[i];
       cords[0] = (((((*IP) * (*ConTo->XScaleFactor)) + *ConTo->XOffset) - *ConFrom->XOffset) / *ConFrom->XScaleFactor);
       cords[1] = ((((*(IP + 1) * (*ConTo->YScaleFactor)) + *ConTo->YOffset) - *ConFrom->YOffset) / *ConFrom->YScaleFactor);
       cords[2] = ((((*(IP + 2) * (*ConTo->ZScaleFactor)) + *ConTo->ZOffset) - *ConFrom->ZOffset) / *ConFrom->ZScaleFactor);
       memcpy(&ConTo->Data[i], &cords, 3*sizeof(int32_t));
   }
   *ConTo->XScaleFactor = *ConFrom->XScaleFactor; *ConTo->XOffset = *ConFrom->XOffset;
   *ConTo->YScaleFactor = *ConFrom->YScaleFactor; *ConTo->YOffset = *ConFrom->YOffset;
   *ConTo->ZScaleFactor = *ConFrom->ZScaleFactor; *ConTo->ZOffset = *ConFrom->ZOffset;
   if(Min){
       *ConTo->MinX = *ConFrom->MinX;
       *ConTo->MinY = *ConFrom->MinY;
       *ConTo->MinZ = *ConFrom->MinZ;
   }
   if(Max){
       *ConTo->MaxX = *ConFrom->MaxX;
       *ConTo->MaxY = *ConFrom->MaxY;
       *ConTo->MaxZ = *ConFrom->MaxZ;

   }
}


/* From a presence/absence vector<bool> voxels, measures the shrub voxel volumes within a defined space around surface points given output
vectors from Voxel_bool function (after using LasScaleConver min = true & max = false
Fast, but could be memory expensive depending on the PC point cloud
   Inputs: SP = voxel id from Voxels_ID
       BV = a pointer to a bool vector for the shrub point cloud from Voxels_bool
       XGC = number of X cells in the cube (1D)
       YGC = number of Y cells in the cube (1D)
   Output; SV = Shrub volume around voxel(SP)
   */
uint32_t SVV_bool(uint32_t SP, std::vector<bool> *BV, uint32_t XGC, uint32_t YGC){
   uint32_t SV = 0;
   for (uint16_t z = 5; z < 51 ; z++){
       for (int16_t y = -20; y < 21; y++){
           for (int16_t x = -20; x < 21; x++){
                 if(BV->at(SP + z*XGC*YGC + y*XGC + x) == true){
                     SV++;
                 }
           }
       }
   }
   return(SV);
}


/* Same as "SVV_bool" but adapted for reference points  */
uint32_t SVV_boolR(uint32_t SP, std::vector<bool> *BV, uint32_t XGC, uint32_t YGC){
   uint32_t SV = 0;
   for (int16_t z = -51; z < 5 ; z++){
       for (int16_t y = -20; y < 21; y++){
           for (int16_t x = -20; x < 21; x++){
                 if(BV->at(SP + z*XGC*YGC + y*XGC + x) == true){
                     SV++;
                 }
           }
       }
   }
   return(SV);
}


/* From a vector of occupied voxel IDs, measures the shrub voxel volumes within a defined space around sphagnum surface points given output vectors
from Voxel_ID funciton (after using LasScaleConv with min = true)
       Slow, very, but possibly saves memory depending on the PC point cloud
   inputs: SP = voxel ID of the sphagnum point being evaluated
       PCF = a pointer to the output of Voxel_ID vector from the "shrub" cloud"
       XGC = # of cells in the X grid
       YGC = # of cells in the Y grid
   output: SV = Number of occupied voxels in the given volume around the sphagnum point
   */
uint32_t SVV_ID(uint32_t SP, std::vector<uint32_t> *PCF, uint32_t XGC, uint32_t YGC){
   uint32_t SV = 0;
   uint32_t Comps[3];
   Comps[0] = floor(SP/(XGC*YGC));
   Comps[1] = floor((SP - (Comps[0]*XGC*YGC))/(YGC));
   Comps[2] = SP - ((Comps[0]*XGC*YGC) + (Comps[1]*YGC));
   for (uint32_t i = 0; i < PCF->size(); i++){
       for (uint16_t z = 5; z  < 51; z++){
           for (int16_t y = -20; y < 21; y++){
                 for (int16_t x = -20; x < 21; x++){
                     if( (Comps[0]*XGC*YGC + z*(XGC*YGC) + Comps[1]*YGC + y*(XGC) + Comps[2] + x) == PCF->at(i)){
                           SV++;
                     }
                 }
           }
       }
   }
   return(SV);
}


/* produces shrub volume for every point in a LasFile
   inputs: LF = pointer to the LasFile to be processed
       BV = vector of bool from Voxels_bool performed on LasFile containing vegetation (raw point cloud)
       XGC = number of X cells in the cube (1D)
       YGC = number of Y cells in the cube (1D)
       Csize = cell size in meters

   outpout: OPV = a vector of voxel volumes for each point in LF (in the order they are in LF)
*/
std::vector<uint32_t> LasShrubVol(LasFile *LF, std::vector<bool> *BV, uint32_t XGC, uint32_t YGC, double Csize){
   std::vector<uint32_t> OPV;
   OPV.resize(*LF->NumberOfPointRecords);
   uint32_t MaxVol = 45*40*40;
   std::vector<uint32_t> VID = Voxels_ID(LF, Csize, false);
   for (uint32_t i = 0 ; i < VID.size(); i++){
         OPV[i] = SVV_bool(VID[i], BV, XGC, YGC);
       if(MaxVol < OPV[i]) {
                 std::cout << "!! ~ WARNING ~ !! ~~~  MAX VOLUME EXCEDED  ~~~ !! ~ WARNING ~ !!"  << std::endl;
       }
   }
   return(OPV);
}


/* takes a pointer to a LasFile and returns a vector of all Z(elevation) values in UTM
   inputs: LF = a pointer to the LasFile to extract elevations from
   output: OPV = vector of elevations from LasFile
*/
std::vector<double> ElevProfileR(LasFile * LF){
   std::vector<double> OPV(*LF->NumberOfPointRecords);
   uint32_t EZI;
   for (uint32_t i = *LF->OffsetToPoints; i < LF->Data.size(); i += *LF->PointDataRecordLength){
       EZI = (i - *LF->OffsetToPoints) / *LF->PointDataRecordLength;
       OPV[EZI] = CoordLongToDbl(*(int32_t*) &LF->Data[i], LF, 'z');
   }

   return(OPV);
}


/* puts points from LasFile into a 2D grid of specified cell size, if sort == true will sort each grid cell by Z value.. size of 
   inputs: Lasf =  pointer to LasFile to be processed
       GSize = desired cell size in meters
       sort = bool if true sorts grid cells by z values... currently sorting not working..
   output: OP.Grid = 2D grid of all points in order x1,y1,z1,x2,y2,z2 etc... sorted or unsorted depending on value of 'sort'
       OP.Lasf = A pointer to the LasFile from which OP.Grid was created
*/
Las2DO Las2DGrid(LasFile * Lasf, double GSize, bool sort){
   Las2DO OP;
   OP.Lasf = Lasf;
   std::vector<uint32_t> GCS = GridCounts(Lasf, GSize);
   OP.Grid.resize(GCS[0]*GCS[1]);
   OP.Inds.resize(GCS[0]*GCS[1]);
   int32_t xpos, ypos, Xcellsize, Ycellsize;
   Xcellsize = GSize / *Lasf->XScaleFactor; Ycellsize = GSize /  *Lasf->YScaleFactor;
   std::vector<int32_t> MinMax = ScaleMinMax(Lasf);
   int32_t* CP;
   for(uint32_t i = *Lasf->OffsetToPoints;i <  Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
       CP = (int32_t*) &Lasf->Data[i];
       xpos = gridpos(MinMax[0], Xcellsize, *CP);
       ypos = gridpos(MinMax[2], Ycellsize, *(CP+1));
       OP.Grid[GridPos2D(xpos, GCS[1] - 1 - ypos, GCS[0])].push_back(*CP);
       OP.Grid[GridPos2D(xpos, GCS[1] - 1 - ypos, GCS[0])].push_back(*(CP+1));
       OP.Grid[GridPos2D(xpos, GCS[1] - 1 - ypos, GCS[0])].push_back(*(CP+2));
      OP.Inds[GridPos2D(xpos, GCS[1] - 1 - ypos, GCS[0])].push_back(i);
   }
   std::vector<std::vector<int32_t> > SG(OP.Grid.size());
   if(sort){
      for(uint32_t i = 0; i < OP.Grid.size(); i++){
           std::sort(OP.Grid[i].begin(),OP.Grid[i].end());
           for(uint32_t j = 2; j < OP.Grid[i].size()-1; j += 3){
                 if(OP.Grid[i].size() > 0){
                     SG[i].push_back(OP.Grid[i][2]);
                     if(OP.Grid[i][j]!=OP.Grid[i][j+1]){
                           SG[i].push_back(OP.Grid[i][j+1]);
                     }
                 }
           }
       }
   OP.Grid = SG;
   }
   return(OP);
}


/* puts elevations from LasFile into a 2D vector 
   inputs: Lasf = a pointer to the LasFile to process
       Gsize = desired cell size in meters
       sort = bool, if true sorts cell Zs in each cell
       Ascend = bool if true sorts ascending (OP[0] = smallest), if false sorts descending (OP[0] = largest)
   outputs: a 2D vector that stores elevations of point in that grid cell (in meters)
   
           NOTE: if sort == false, Ascend inconsequential
            if sort = true, indexes are not sorted
*/
Las2DZO Las2DZs(LasFile *Lasf, double GSize, bool sort, bool Ascend){
   Las2DZO OP;
   std::vector<uint32_t> GCS = GridCounts(Lasf, GSize);
   OP.Z.resize(GCS[0]*GCS[1]); OP.inds.resize(GCS[0]*GCS[1]);
   int32_t xpos, ypos, Xcellsize, Ycellsize;
   Xcellsize = GSize / *Lasf->XScaleFactor; Ycellsize = GSize /  *Lasf->YScaleFactor;
   std::vector<int32_t> MinMax = ScaleMinMax(Lasf);
   int32_t* CP;
   for(uint32_t i = *Lasf->OffsetToPoints; i <  Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
      CP = (int32_t*) &Lasf->Data[i];
       xpos = gridpos(MinMax[0], Xcellsize, *CP);
       ypos = gridpos(MinMax[2], Ycellsize, *(CP+1));
         OP.Z[GridPos2D(xpos, GCS[1] - ypos - 1, GCS[0])].push_back(CoordLongToDbl(*(CP+2), Lasf, 'z'));
      OP.inds[GridPos2D(xpos, GCS[1] - ypos - 1, GCS[0])].push_back(i);
   }
   if(sort){
       for (uint32_t i = 0; i < OP.Z.size(); i++){
           if(Ascend){
                 std::sort(OP.Z[i].begin(), OP.Z[i].end());

           } else {
                 std::sort(OP.Z[i].rbegin(), OP.Z[i].rend()); }
       }
   }
   return(OP);
}


/* takes LasFile and returns maximum returns (UTM) in each grid cell. reprojected into 1D. Has a cleaning option to deal with "ghost" points
   inputs: Lasf = a pointer to the LasFile to be processed
       Gsize = grid size in meters
       clean = bool of whether to perform cleaning operation
       tol = the tolerance to be used for the cleaning operation in meters
   output: OP = vector of maximum returns in each cell in meters
           NOTE = if clean == false, tol is inconsequential
           NOTE = cleaning not very good.. needs work
*/
std::vector<double> Las2DZMax(LasFile *Lasf, double Gsize, bool clean, double tol){
   Las2DZO Z2D =  Las2DZs(Lasf, Gsize, true, false);
   std::vector<double> OP(Z2D.Z.size());
   if(clean){
       for(uint32_t i = 0; i < Z2D.Z.size(); i++){
           if(Z2D.Z[i].size() > 0){
                 OP[i] = Z2D.Z[i][0];
                 for(uint32_t j = 0; j < Z2D.Z[i].size() - 1; j++){
                     if(Z2D.Z[i][j] < Z2D.Z[i][j + 1] + tol){
                           OP[i] = Z2D.Z[i][j + 1];
                     }
                 }
           }
       }
   } else {
       for(uint32_t i = 0; i < Z2D.Z.size(); i++){
           if(Z2D.Z[i].size() > 0){
                 OP[i] = Z2D.Z[i][0];
           }
       }
   }
   return(OP);
}


/* takes LasFile and returns minimum returns (UTM) in each grid cell. reprojected into 1D. Has a cleaning option to deal with "ghost" points
   inputs: Lasf = a pointer to the LasFile to be processed
       Gsize = grid size in meters
       clean = bool of whether to perform cleaning operation
       tol = the tolerance to be used for the cleaning operation in meters
   output: OP = vector of minimum returns in each cell in meters
           NOTE = if clean == false, tol is inconsequential
           NOTE = cleaning not very good.. needs work
*/
std::vector<double> Las2DZMin(LasFile *Lasf, double Gsize, bool clean, double tol){
   Las2DZO Z2D =  Las2DZs(Lasf, Gsize, true, false);
   std::vector<double> OP(Z2D.Z.size());
   if(clean){
       for(uint32_t i = 0; i < Z2D.Z.size(); i++){
           if(Z2D.Z[i].size() > 0){
                 OP[i] = Z2D.Z[i][0];
                 for(uint32_t j = 0; j < Z2D.Z[i].size() - 1; j++){
                     if(Z2D.Z[i][j] > Z2D.Z[i][j + 1] + tol){
                           OP[i] = Z2D.Z[i][j + 1];
                     }
                 }
           }
       }
   } else {
      for(uint32_t i = 0; i < Z2D.Z.size(); i++){
         if(Z2D.Z[i].size() > 0){
            OP[i] = Z2D.Z[i].back();
         }
      }
   }
   return(OP);
}

/* Performs Las2DZMax and returns a LasFile with the max points
   inputs: Lasf = pointer to LasFile to be performed on
       Gsize = grid size in meters
       clean = bool of whether to perform cleaning operation
       tol = the tolerance to be used for the cleaning operation in meters
   output: OP = vector of maximum returns in each cell in meters
           NOTE = if clean == false, tol is inconsequential      
*/
LasFile Las2DZMaxLasOut(LasFile *Lasf, double Gsize, bool clean, double tol){
   std::vector<double> TV = Las2DZMax(Lasf, Gsize, clean, tol);
   LasFile OP;
   return(OP);  //  THIS IF OBVIOUSLY NOT COMPLETE...
}


/* Performs Las2DZMin and returns a LasFile with the min points
   inputs: Lasf = pointer to LasFile to be performed on
       Gsize = grid size in meters
       clean = bool of whether to perform cleaning operation
       tol = the tolerance to be used for the cleaning operation in meters
   output: OP = vector of maximum returns in each cell in meters
           NOTE = if clean == false, tol is inconsequential      
*/
LasFile Las2DZMinLasOut(LasFile *Lasf, double Gsize, bool clean, double tol){
   
   std::vector<uint32_t> GCS = GridCounts(Lasf, Gsize); 
   std::vector<uint64_t> MinIDs(GCS[0] * GCS[1]);
   std::vector<int32_t> MinMax = ScaleMinMax(Lasf);
   int32_t * CP, MP;
   int32_t xpos, ypos, Xcellsize, Ycellsize, xypos;
   Xcellsize = Gsize / *Lasf->XScaleFactor;
   Ycellsize = Gsize / *Lasf->YScaleFactor;
   for(size_t i = *Lasf->OffsetToPoints; i < Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
   CP = (int32_t*) &Lasf->Data[i]; 
   xpos = gridpos(MinMax[0], Xcellsize, *CP); 
   ypos = gridpos(MinMax[2], Ycellsize, *(CP+1)); 
   xypos = GridPos2D(xpos, ypos, GCS[0]);
   if(!MinIDs[xypos]){
      MinIDs[xypos] = i;
   } else if (*(CP + 2) < *(int32_t*) &Lasf->Data[MinIDs[xypos] + 2*sizeof(int32_t)]){
       MinIDs[xypos] = i;
   }
   }
   LasFile OP;
   OP.Data.resize(*Lasf->OffsetToPoints + MinIDs.size()*(*Lasf->PointDataRecordLength));
   memcpy(&OP.Data[0], &Lasf->Data[0], *Lasf->OffsetToPoints);
   LasFormat(&OP);
   size_t ctr = 0;
   for(size_t i = 0; i < MinIDs.size(); i++){
   if(MinIDs[i]){
      ctr++;
      memcpy(&OP.Data[*Lasf->OffsetToPoints + ctr*(*Lasf->PointDataRecordLength)], &Lasf->Data[MinIDs[i]], *Lasf->PointDataRecordLength);
   }
   }
   OP.Data.resize(*OP.OffsetToPoints + ctr*(*Lasf->PointDataRecordLength));
   return(OP);  
}



/* Takes vector of Max returns from Las2DZMax and convolves with a exponential smoothing kernel 
   inputs: Lasf = a pointer to the LasFile from which ZM came
       Gsize = the grid cell size used to produce ZM
       clean = bool for if cleaning operation should be perfomred
       tol = the tolerance to be used to clean (meters)
   outputs: OP = convolved grid values (2D projected to 1D)
*/
std::vector<double> ExpFilt(LasFile * Lasf, double Gsize, bool clean, double tol){
   std::vector<double> ZM = Las2DZMax(Lasf, Gsize, clean, tol);
   std::vector<double> OP(ZM.size());
   std::vector<uint32_t> GCS = GridCounts(Lasf, Gsize);
   std::vector<double> tmp;
   int32_t yg, xg;
   double dist;
   uint32_t pc;
   for(uint32_t i = 0; i < ZM.size(); i++){
       yg = floor(i / GCS[0]);
       xg = i - yg*GCS[0];
       pc = 0;
       for(int32_t j = -15; j < 16; j++){
           for(int32_t k = -15; k < 16; k++){
                 if( ( (((yg*GCS[0]) + (j * GCS[0])) + (xg + k) ) < OP.size()) & ( (((yg*GCS[0]) + (j * GCS[0])) + (xg + k) ) > 0)     ){
                     dist = DP2PSI(0, 0, j, k) * Gsize;
                     if(ZM[((yg*GCS[0]) + (j * GCS[0])) + (xg + k)] != 0){
                           OP[i] = OP[i] + ((ZM[((yg*GCS[0]) + (j * GCS[0])) + (xg + k)]) - *Lasf->MinZ) * (1/pow(abs((15*Gsize)-dist)+1, .5));
                           pc = pc + 1;
                     }
                 }
           }
       }
       OP[i] = OP[i] / pc;
       //std::cout << pc << " " << OP[i] << std::endl;
   }
   return(OP);
}


/* Takes a LasFile pointer and performs a gaussian convolution on the maximum returns in the grid
   inputs: Lasf = a pointer to the LasFile from which ZM came
       Gsize = the grid cell size used to produce ZM
       clean = bool for if cleaning operation should be perfomred
       tol = the tolerance to be used to clean (meters)
   outputs: OP = convolved grid values (2D projected to 1D)
*/
std::vector<double> LasMaxGausFilt(LasFile * Lasf, double Gsize, bool clean, double tol){
   std::vector<double> ZM = Las2DZMax(Lasf, Gsize, clean, tol);
   std::vector<double> OP(ZM.size());
   std::vector<uint32_t> GCS = GridCounts(Lasf, Gsize);
   int32_t yg, xg;
   int32_t sig = 1;
   for(uint32_t i = 0; i < ZM.size(); i++){
       yg = floor(i / GCS[0]);
       xg = i - yg*GCS[0];
       for(int32_t j = -4*sig; j < 4*sig + 1; j++){
           for(int32_t k = -4*sig; k < 4*sig + 1; k++){
                 if( ( (((yg*GCS[0]) + (j * GCS[0])) + (xg + k) ) < OP.size()) && ( (((yg*GCS[0]) + (j * GCS[0])) + (xg + k) ) > 0)    ){
                     if(ZM[((yg*GCS[0]) + (j * GCS[0])) + (xg + k)] != 0){
                           OP[i] = OP[i] + ((ZM[((yg*GCS[0]) + (j * GCS[0])) + (xg + k)]) - *Lasf->MinZ) * (1/(2*M_PI*(sig*sig))) * exp(-((k*k + j*j) / (2*sig*sig)));
                     } 
                 }
           }
       }
   }
   return(OP);
}


/* Takes LasFile pointer and performs a laplacian of gaussain (LoG) convolution on maximum grid returns
   inputs: Lasf = a pointer to the LasFile from which ZM came
       Gsize = the grid cell size used to produce ZM
       clean = bool for if cleaning operation should be perfomred
       tol = the tolerance to be used to clean (meters)
   outputs: OP = convolved grid values (2D projected to 1D)
*/
std::vector<double> LapFilt(LasFile * Lasf, double Gsize, bool clean, double tol){
   std::vector<double> ZM = Las2DZMax(Lasf, Gsize, false, tol);
   std::vector<double> OP(ZM.size());
   std::vector<uint32_t> GCS = GridCounts(Lasf, Gsize);
   std::vector<double> tmp;
   int32_t yg, xg;
   uint32_t pc;
   int32_t sig = 2;
   for(uint32_t i = 0; i < ZM.size(); i++){
       yg = floor(i / GCS[0]);
       xg = i - yg*GCS[0];
       pc = 0;
       for(int32_t j = -4*sig; j < 4*sig + 1; j++){
           for(int32_t k = -4*sig; k < 4*sig + 1; k++){
                 if( ( (((yg*GCS[0]) + (j * GCS[0])) + (xg + k) ) < OP.size()) & ( (((yg*GCS[0]) + (j * GCS[0])) + (xg + k) ) > 0)     ){
                     if(ZM[((yg*GCS[0]) + (j * GCS[0])) + (xg + k)] != 0){
                           OP[i] = OP[i] + ((ZM[((yg*GCS[0]) + (j * GCS[0])) + (xg + k)]) - *Lasf->MinZ) * ((-1/(M_PI*pow(sig,4))) * (1 - ((k*k + j*j) / (2*sig*sig))) * exp(-((k*k + j*j) / (2*sig*sig))));
                           pc = pc + 1;
                     }
                 }
           }
       }
   }
   if(clean){
       for(uint32_t i = 0; i < OP.size(); i++){
           if(OP[i] < -.5){
                 OP[i] = 100;
           } else {
                 OP[i] = 0;
           }
       }

   }
   return(OP);
}

/* Adjusts the coordinate minimum(s) LasFile generated by voxelization or griding that are going to be placed in a new grid or voxel grid. Need to do this to this to make center of grid cells correspond to where points have been placed during rasterization/voxelization.
   inputs: 
      Lasf = pointer to the LasFile being evaluated 
      GS = size of the grid cells being used
      x = bool.. if true -> adjust MinX
      y = bool.. if true -> adjust MinY
      z = bool.. if true -> adjust MinZ
     ~~~~~~~~~~~ NOTE!!! this modifies the min values in memory NOT in the actual file, if the file is written LasClean should be run to correctly assign the minimum values in the Las header.  ~~~~~~~~~~~~~~

*/
void VoxelGridMinAdjust(LasFile *Lasf, double GS, bool x, bool y, bool z){
   if(x){
   *Lasf->MinX = *Lasf->MinX - (GS/2);
   }
   if(y){
   *Lasf->MinY = *Lasf->MinY - (GS/2);
   }
   if(z){
   *Lasf->MinZ = *Lasf->MinZ - (GS/2);
   }
}


/* Adjusts the Min/Max values in a LasHeader. Useful for adding padding to a rasterized point cloud. x/y/zpad add grid cells before and after the original data. Option for performing VoxelGridMinAdjust.
   inputs:
      Lasf = pointer to to the LasFile being evaluated
      GS = size of grid cells being used
      adj = if true, performs VoxelGridMinAdjust
      xpad = number of grid cells to pad orignial data with in the x direction. This number of grid cells will be added to both sides of the data.
      ypad = number of grid cells to pad orignial data with in the y direction. This number of grid cells will be added to both sides of the data.
      zpad = number of grid cells to pad orignial data with in the z direction. This number of grid cells will be added to both sides of the data.
   ~~~~~~~~~~~ NOTE!!! this modifies the min & max values in memory NOT in the actual file, if the file is written LasClean should be run to correctly assign the minimum & maximum values in the Las header.  ~~~~~~~~~~~~~~

*/
void LasGridPad(LasFile *Lasf, double GS, bool adj, int32_t xpad, int32_t ypad, int32_t zpad){
   if(adj){
      VoxelGridMinAdjust(Lasf, GS, true, true, true);
   }
   *Lasf->MinX = *Lasf->MinX - xpad*GS;
   *Lasf->MinY = *Lasf->MinY - ypad*GS;
   *Lasf->MinZ = *Lasf->MinZ - zpad*GS;
   *Lasf->MaxX = *Lasf->MaxX + xpad*GS;
   *Lasf->MaxY = *Lasf->MaxY + ypad*GS; 
   *Lasf->MaxZ = *Lasf->MaxZ + zpad*GS;
}


/* Sets the classification of points from a LasFile based on a classification grid (vector) (i.e., places points in a grid and classify based on "some" metric)
   intputs:
   Lasf = a pointer to the LasFile to be classified 
   CV = a pointer to the vector (grid) representing classsification.. needs to be a vector of doubles here
   GS = the grid cell size the classification grid was made using
   outputs: NA
*/
void LasClassByClassGrid(LasFile *Lasf, std::vector<double> *CV, double GS){
   int32_t ClassByte = 15;
   std::vector<int32_t> PointGridID = Las2DGridPositions(Lasf, GS);
   int32_t Ofs = *Lasf->OffsetToPoints + ClassByte; // THIS is not correct needs to be adjusted to where the "class" is specified in las 1.2
   unsigned char Class1 = 0;
   unsigned char Class2 = 1;
   for(size_t i = 0; i < PointGridID.size(); i++){
      if(CV->at(PointGridID[i])){
         memcpy(&Lasf->Data[Ofs + i*(*Lasf->PointDataRecordLength)], &Class1, sizeof(unsigned char));
      } else {
         memcpy(&Lasf->Data[Ofs + i*(*Lasf->PointDataRecordLength)], &Class2, sizeof(unsigned char));
      }
   }
}


/* Returns a vector of 2D grid positions (Projected in 1D) for each point in a LasFile
   inputs:
   Lasf = a pointer to the las file to be operated on
   GS = grid cell size
   
   outputs:
   OP = vector containing 2D grid positions for each point record in order they occur in the LasFile
*/
std::vector<int32_t> Las2DGridPositions(LasFile *Lasf, double GS){
   std::vector<uint32_t> GCS = GridCounts(Lasf, GS);
   std::vector<int32_t> OP(*Lasf->NumberOfPointRecords); 
   int32_t xpos, ypos, Xcellsize, Ycellsize;
   Xcellsize = GS / *Lasf->XScaleFactor;
   Ycellsize = GS / *Lasf->YScaleFactor;
   std::vector<int32_t> MinMax = ScaleMinMax(Lasf);
   int32_t* CP;
   for(uint32_t i = *Lasf->OffsetToPoints; i <  Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
      CP = (int32_t*) &Lasf->Data[i];
      xpos = gridpos(MinMax[0], Xcellsize, *CP);
      ypos = gridpos(MinMax[2], Ycellsize, *(CP+1));
      OP[ (i - *Lasf->OffsetToPoints) / *Lasf->PointDataRecordLength] = GridPos2D(xpos, ypos, GCS[0]);
   } 
   return(OP);   
}

/* takes a vector of type T (needs to be numeric) and inserts the values into the user data field of a LasFile 
   inputs:
   Lasf = pointer to the LasFile the user data is to be filled in
   VIN = pointer to the vector of values to put into user data
   GS = grid cell size
   outputs: NA
*/
template <typename T> 
void LasScalarForUsrData(LasFile *Lasf, std::vector<T> *VIN){
   if(VIN->size() != (Lasf->Data.size() - *Lasf->OffsetToPoints) / *Lasf->PointDataRecordLength){
      std::cout << "!!!!! WARNING! vector length does not match the number of point records in the LasFile !!!!!!! " << std::endl;
   } else {    
      uint32_t UDataByte = 17; // just for las 1.2 and whatever point record version cloud compare outputs 
      int32_t Ofs = *Lasf->OffsetToPoints + UDataByte;
      T vmax = *max_element(std::begin(*VIN), std::end(*VIN));
      T vmin = *min_element(std::begin(*VIN), std::end(*VIN)); 
      T R = vmax - vmin;
      double SF = R / 255;
      std::vector<unsigned char> Udat(VIN->size()); 
      for(size_t i = 0; i < Udat.size(); i++){
         Udat[i] = (VIN->at(i) - vmin) / SF;
      }
      for(size_t i = 0; i < Udat.size(); i++){
         memcpy(&Lasf->Data[Ofs + i*( *Lasf->PointDataRecordLength)], &Udat[i], sizeof(unsigned char));
      }
   }
}



/* puts a the values from a vector<double> into the GPS time field of a LasFile
   inputs:
   Lasf = pointer to the LasFile the GPS time is to be overwritten
   VIN = pointer to the vector<double> to place in LasFile GPS time
   GS = grid cell size of VIN   

   outputs: NA

*/
//template <typename T> 
void LasScalarForGPSTime(LasFile *Lasf, std::vector<double> *VIN){
   if(VIN->size() != (Lasf->Data.size() - *Lasf->OffsetToPoints) / *Lasf->PointDataRecordLength){
      std::cout << "!!!!! WARNING! vector length does not match the number of point records in the LasFile !!!!!!! " << std::endl;
   } else {  
      uint32_t GPSbyte = 20; // just for las 1.2 and whatever point record version cloud compare outputs 
      int32_t Ofs = *Lasf->OffsetToPoints + GPSbyte;
      for(size_t i = 0; i < VIN->size(); i++){
         memcpy(&Lasf->Data[Ofs + i*(*Lasf->PointDataRecordLength)], &VIN->at(i), sizeof(double));
      }
   }
}


/* Extracts points from a LasFile that are contained in a polygon given its vertices 2d vertices (e.i., x & y coordinates)
   inputs:
   Lin = pointer to the LasFile to be processed
   nvert = number of vertices in the polygon being assessed
   XV = pointer to x components of polygon vertices
   YV = pointer to y components of polygon vertices
   outputs:
   OP = new LasFile containing only points within the polygon


*/
LasFile ExtractLasInPolygon(LasFile * Lin, int32_t nvert, double * XV, double * YV){
   LasFile OP;
   OP.Data.resize(Lin->Data.size());
   int32_t * ptr;
   ptr = (int32_t*) &Lin->Data[*Lin->OffsetToPoints]; 
   int32_t In;
   size_t np = 0;
   std::vector<double> Tmp; 
   memcpy(&OP.Data[0], &Lin->Data[0], *Lin->OffsetToPoints);
   LasFormat(&OP);
   for(size_t i = *Lin->OffsetToPoints; i < Lin->Data.size(); i+= *Lin->PointDataRecordLength){
      ptr = (int32_t*) &Lin->Data[i];
      Tmp = PointLongToDbl(*ptr, *(ptr + 1), *(ptr + 2), Lin);
      In = PointInPolygon(8, XV, YV, Tmp[0], Tmp[1]);
      if(In){
         memcpy(&OP.Data[*OP.OffsetToPoints + np*(*OP.PointDataRecordLength)], &Lin->Data[i], *OP.PointDataRecordLength);
         np++;
      }
   }
   OP.Data.resize(*OP.OffsetToPoints + np*(*OP.PointDataRecordLength));
   *OP.NumberOfPointRecords = np; 
   return(OP);
}



/* Calculates the voxel volume in a given radius around supplied xy points
   I



*/
std::vector<uint32_t> PointVoxVolRadius(LasFile *Lasf, std::vector<double> * PV, double rad){
   std::vector<uint32_t> OP(PV->size()/2);
   LasFormat(Lasf);
   int32_t * CP;
   size_t ctr;
   for(size_t i = 0; i < PV->size(); i += 2){
      ctr = 0;
      for(size_t j = *Lasf->OffsetToPoints; j < Lasf->Data.size(); j += *Lasf->PointDataRecordLength){
         CP = (int32_t*) &Lasf->Data[j];
         if( DP2P(CoordLongToDbl(*CP, Lasf, 'x'), CoordLongToDbl(*(CP + 1), Lasf, 'y'), PV->at(i), PV->at(i+1)) < rad ){
            ctr++;
         }
      }
      OP[i/2] = ctr;
   }
   return(OP);
}


/* Calculates the voxel volume in a given radius around supplied xy points
   I



*/
LasFile PointVoxVolRadiusLas(LasFile *Lasf, std::vector<double> * PV, double rad){
   LasFile OP;
   OP.Data.resize(Lasf->Data.size()*3);
//   std::vector<double> PIDs((Lasf->Data.size() - *Lasf->OffsetToPoints) / *Lasf->PointDataRecordLength);
   std::vector<double> PIDs;
   size_t NoP = (Lasf->Data.size() - *Lasf->OffsetToPoints) / *Lasf->PointDataRecordLength;
   LasFormat(Lasf);
   memcpy(&OP.Data[0], &Lasf->Data[0], *Lasf->OffsetToPoints);
   LasFormat(&OP);
   LasFormat(Lasf);
   uint32_t otp = *OP.OffsetToPoints;
   uint32_t pdrl = *OP.PointDataRecordLength;
   int32_t * CP;
   size_t ctr = 0;
   size_t RC = 2;
   std::vector<char> tmp;
   for(size_t i = 0; i < PV->size(); i += 2){
      std::cout << (i/2) + 1 << std::endl;
      for(size_t j = *Lasf->OffsetToPoints; j < Lasf->Data.size(); j += *Lasf->PointDataRecordLength){
         CP = (int32_t*) &Lasf->Data[j];
         if( DP2P(CoordLongToDbl(*CP, Lasf, 'x'), CoordLongToDbl(*(CP + 1), Lasf, 'y'), PV->at(i), PV->at(i+1)) < rad ){
           if(*OP.OffsetToPoints + ctr*(*OP.PointDataRecordLength) > OP.Data.size()){
              OP.Data.resize(RC*NoP*pdrl + otp);
              LasFormat(&OP);
              
              std::cout << *OP.OffsetToPoints + ctr*(*OP.PointDataRecordLength) << " " << OP.Data.size() << std::endl;
           }
           memcpy(&OP.Data[*OP.OffsetToPoints + ctr*(*OP.PointDataRecordLength)], &Lasf->Data[j], *OP.PointDataRecordLength);
           PIDs.push_back((i/2) + 1);
           ctr++;
         }
      }
   } 
   OP.Data.resize(*OP.OffsetToPoints + ctr*(*OP.PointDataRecordLength));
   LasScalarForGPSTime(&OP, &PIDs);
   return(OP);
}



/* Voxelizes a las file and outputs the voxelized version
   inputs:
      LasIn = point to the LasFile to be voxelized
      VS = size of voxel to be voxelized with

   outputs: 
      LasOut = the voxelized point cloud      

*/

LasFile VoxelizeLas2Las(LasFile * LasIn, double VS){
   LasFile LasOut;
   LasOut.Data.resize(*LasIn->OffsetToPoints);
   memcpy(&LasOut.Data[0], &LasIn->Data[0], *LasIn->OffsetToPoints);
   LasFormat(&LasOut);
   std::vector<uint32_t> GCS = GridCounts(LasIn, VS);
   std::vector<uint32_t> VID = Voxels_ID(LasIn, VS, true);
   std::vector<double> VC = VoxCoord(&VID, GCS[0], GCS[1], VS, LasIn);
   std::vector<int32_t> tmpv;
   LasOut.Data.resize(*LasIn->OffsetToPoints + (*LasIn->PointDataRecordLength * VID.size()));
   LasFormat(&LasOut);
   for(size_t i = 0; i < VC.size() - 3; i+=3){
      tmpv = PointDblToLong(VC[i], VC[i + 1], VC[i + 2], LasIn);
      memcpy(&LasOut.Data[*LasOut.OffsetToPoints + ((i/3)*(*LasOut.PointDataRecordLength))], &tmpv[0], tmpv.size()*sizeof(int32_t));
   }
   return(LasOut);
}













/*

colors point cloud of input domain for Hollow index... needs work 

*/
int IDK(LasFile *Lasf, std::vector<double> *PV, std::string BWF, double rad, double gs){
   int OP;
   std::vector<uint32_t> GCS = GridCounts(Lasf, gs);
   std::vector<bool> FullVec(GCS[0] * GCS[1], true);
   std::vector<double> BWV = ReadTxt(BWF);
   std::vector<double> XV(BWV.size()/2);
   std::vector<double> YV(XV.size());
   for(size_t i = 0; i < BWV.size()/2; i+=2){
      XV[i/2] = BWV[i];
      YV[(i/2) + 1] = BWV[i+1]; 
   }
   LasFile BG;
   BG.Data.resize(GCS[0] * GCS[1] * (*Lasf->PointDataRecordLength));
   LasFormat(&BG);
   memcpy(&BG.Data[0], &Lasf->Data[0], *Lasf->OffsetToPoints);
   int32_t X, Y;
   for(size_t i = 0; i < GCS[0]; i++){
      for(size_t j = 0; j < GCS[1]; j++){
         X = CoordDblToLong(*BG.MinX, &BG, 'x') * (i*(0.01 * *BG.XScaleFactor));
         Y = CoordDblToLong(*BG.MinY, &BG, 'y') * (j*(0.01 * *BG.YScaleFactor));
         memcpy(&BG.Data[*BG.OffsetToPoints + (i*(*BG.PointDataRecordLength)) + (j*(*BG.PointDataRecordLength))], &X, sizeof(int32_t));
         memcpy(&BG.Data[*BG.OffsetToPoints + (i*(*BG.PointDataRecordLength)) + (j*(*BG.PointDataRecordLength)) + sizeof(int32_t)], &X, sizeof(int32_t));
      }
   }


   LasFile IBW = ExtractLasInPolygon(Lasf, XV.size(), &XV[0], &YV[0]);
   std::vector<uint32_t> VoxVols = PointVoxVolRadius(&IBW, PV, rad);
std::cout << (FullVec.size() * 1000) / pow(10,9) << std::endl;     

   return(OP);
}



/*  Given a las file of reference points and an associated point cloud, finds the nearest point to all references points and calculates the z difference, also outputs a lasfile of the extracted closet points from the point cloud


*/
REO Las2LasZErrorClosestXY(LasFile * RP, LasFile * PC){
   REO GO;
   LasFile LasO;
   LasO.Data.resize(RP->Data.size());
   memcpy(&LasO.Data[0], &PC->Data[0], *PC->OffsetToPoints);
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
      for(uint32_t j = *PC->OffsetToPoints; j < PC->Data.size(); j += *PC->PointDataRecordLength){
         sp = (int32_t*) &PC->Data[j];
         SCoords[0] = CoordLongToDbl(*sp, PC, 'x');
         SCoords[1] = CoordLongToDbl(*(sp + 1), PC, 'y');
         SCoords[2] = CoordLongToDbl(*(sp + 2), PC, 'z');
         dist = DP2P(RCoords[0], RCoords[1], SCoords[0], SCoords[1]);
         if(mindist){
            if(dist < mindist){
               mindist = dist;
               MinI = j;
               spm = (int32_t*) &PC->Data[j];
               closestz = SCoords[2];
            }
         } else {
            mindist = dist;
            MinI = j;
            spm = (int32_t*) &PC->Data[j];
            closestz = SCoords[2];
         }
      }
      memcpy(&LasO.Data[i], &PC->Data[MinI],*PC->PointDataRecordLength);
      error[ctr] = closestz - (RCoords[2]);
      RefEle[ctr] = RCoords[2];
      ReConEle[ctr] = closestz;
      ctr++;
   }
   GO.RefEle = RefEle;
   GO.PCEle = ReConEle;
   GO.LasF = LasO;
   GO.Error = error;
   return(GO);
}







