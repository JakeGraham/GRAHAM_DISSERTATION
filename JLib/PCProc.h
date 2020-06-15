#ifndef PCProc_h
#define PCProc_h
#include<vector> 
#include<iostream>
#include<string>
#include<cmath>
#include<stdint.h>
#include"Misc.h"
#include"LasReadWrite.h"
#include<iterator>


// Converts coordinates in double (UTM) to the scaled int32_t from a LasFile
int64_t CoordDblToLong(double Pin, LasFile *Lasf, char Coord);


// Converts all point coordinates in double to (UTM) to scaled int32_t from a LasFile
std::vector<int32_t> PointDblToLong(double X, double Y, double Z, LasFile *Lasf);


// Converts coordinates in int32_ts from LasFile to double (UTM) using the files scalings & offsets
double CoordLongToDbl(int32_t Pin, LasFile *Lasf, char Coord);


// Converts all point coordinates in scaled int32_t from a LasFile to double (UTM)
std::vector<double> PointLongToDbl(int32_t X, int32_t Y, int32_t Z, LasFile *Lasf);


/* Given a vector of xyz coordinates in as doubles .. output 3x the length of input orderd x1,y1,z1,x2,y2,z2..etc
        input:  voxels = a pointer to a vector of Voxel IDs to be converted to coordinates
                XGC = # of x grid cells used in voxelization
                YGC = # of y grid cells used in voxeliation
                VoxS = the size of voxel used in voxeliation
                Lasf = a pointer to the LasFile from which voxelization occured
        output: voxr = a vector of double containing xyz coordinates for each occupied voxel in format x1,y1,z1,x2,y2,z2..etc
*/
std::vector<double> VoxCoord(std::vector<uint32_t> *voxels, int32_t XGC, int32_t YGC, double VoxS, LasFile *LasF);


/*       converts min/max double values from header to a scaled vector of int32_ts, so it is comparable to point records w/out scaling
        order is: 0:MinX,1:MaxX,2:MinY,3:MaxY,4:MinZ,5:MaxZ                                                             */
std::vector<int32_t> ScaleMinMax(LasFile *LasIn);


// takes a las file and voxel size and returns a vector with numbers of cells in each coordinate order: 0 = X, 1 = Y, 2 = Z
std::vector<uint32_t> GridCounts(LasFile *Lasf, float VSize);


/* Structure to be output but Voxels_bool and Voxels_bool_thresh

*/
struct VoxBO{
	std::vector<bool> BV; // Binary represenation of point cloud projected in 1D
	uint64_t CVV; // voxel volume of the entire point cloud
};


/* takes a pointer to a LasFile and returns a vector<bool> of voxel presence/absence
        inputs: Lasf = a pointer to the LasFile to be voxelized
                Vsize = the voxel size (meters)
        output: Bout = a 1D projected vector of bools; true = occupied & false = !occupied; 1D projection hierachry = z > y > x 

*/
VoxBO Voxels_bool(LasFile *Lasf, float VSize);
 

/* takes a pointer to a LasFile and returns a vector<bool> of voxel presence/absence on all points above or below a threshold
        inputs: Lasf = a pointer to the LasFile to be voxelized
                Vsize = the voxel size (meters)
		Thresh = elevation threshold to be used (meters)
		UpDown = bool that tells whether to voxelize points above or below the threshold: true = voxelize points above thresh,
			false = voxelize points below thresh
        output: Bout = a 1D projected vector of bools; true = occupied & false = !occupied; 1D projection hierachry = z > y > x 

*/
VoxBO Voxels_bool_thresh(LasFile *Lasf, float VSize, double Thresh, bool UpDown);


/* takes a pointer to a LasFile and returns a vector of sorted and reduced voxels (by ID #) that are occupied
        inputs: Lasf = pointer to the LasFile to be voxelized
                VSize = voxel size (m)
                sort = whether the output vector that should sorted so the "lowest" numbered voxel is the first value in the vector
                VV = sorted vector of voxel IDs 
                OcVox = unsorted vector of voxel IDs
*/
std::vector<uint32_t> Voxels_ID(LasFile* Lasf, float VSize, bool sort);


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
std::vector<uint32_t> Voxels_ID_Thresh(LasFile* Lasf, float VSize, bool sort, bool UpDown, double Thresh);


/* converts second las file to have same scaling factors and offsets as first, if Min = TRUE the min xyzs in the header are set to match,
this is so they can be voxelized and voxel numbers will match for comparison; if Max = TRUE the max xyzs in the header are set to match;
by setting Min = true & Max = False you save a large amount of memory when using SVV_bool    
        inputs: ConFrom = LasFile to convert ConTo's scaling to 
                ConTo = LasFile to convert to scaling from ConFrom
                Min =   true = set minimums of bounding boxes equal.. useful for when voxelizing multiple LasFiles and comparing voxels 
                Max =   true = set maximums of bounding boxes equal 
        output: N/A.. the ConTo file will now have the same scaling and mins/maxes (if true) as ConFrom  
*/
void LasScaleConv(LasFile *ConFrom, LasFile *ConTo, bool Min, bool Max);


/* From a presence/absence vector<bool> voxels, measures the shrub voxel volumes within a defined space around surface points given output
vectors from Voxel_bool function (after using LasScaleConver min = true & max = false
Fast, but could be memory expensive depending on the PC point cloud
        Inputs: SP = voxel id from Voxels_ID
                BV = a pointer to a bool vector for the shrub point cloud from Voxels_bool
                XGC = number of X cells in the cube (1D)
                YGC = number of Y cells in the cube (1D)
        Output; SV = Shrub volume around voxel(SP)
        */
uint32_t SVV_bool(uint32_t SP, std::vector<bool> *BV, uint32_t XGC, uint32_t YGC);


/* Same as "SVV_bool" but adapted for reference points  */
uint32_t SVV_boolR(uint32_t SP, std::vector<bool> *BV, uint32_t XGC, uint32_t YGC);


/* From a vector of occupied voxel IDs, measures the shrub voxel volumes within a defined space around sphagnum surface points given output vectors
from Voxel_ID funciton (after using LasScaleConv with min = true)
                Slow, very, but possibly saves memory depending on the PC point cloud
        inputs: SP = voxel ID of the sphagnum point being evaluated
                PCF = a pointer to the output of Voxel_ID vector from the "shrub" cloud"
                XGC = # of cells in the X grid
                YGC = # of cells in the Y grid
        output: SV = Number of occupied voxels in the given volume around the sphagnum point
        */
uint32_t SVV_ID(uint32_t SP, std::vector<uint32_t> *PCF, uint32_t XGC, uint32_t YGC);


/* produces shrub volume for every point in a LasFile
        inputs: LF = pointer to the LasFile to be processed
                BV = vector of bool from Voxels_bool performed on LasFile containing vegetation (raw point cloud)
                XGC = number of X cells in the cube (1D)
                YGC = number of Y cells in the cube (1D)
                Csize = cell size in meters

        outpout: OPV = a vector of voxel volumes for each point in LF (in the order they are in LF)
*/
std::vector<uint32_t> LasShrubVol(LasFile *LF, std::vector<bool> *BV, uint32_t XGC, uint32_t YGC, double Csize);


/* takes a pointer to a LasFile and returns a vector of all Z(elevation) values in UTM
        inputs: LF = a pointer to the LasFile to extract elevations from
        output: OPV = vector of elevations from LasFile
*/
std::vector<double> ElevProfileR(LasFile * LF);


// stucture for Las2DGrid output
struct Las2DO {
        std::vector<std::vector<int32_t> > Grid;
	std::vector<std::vector<size_t> > Inds;
        LasFile *Lasf;
};


/* puts points from LasFile into a 2D grid of specified cell size, if sort == true will sort each grid cell by Z value.. size of 
        inputs: Lasf =  pointer to LasFile to be processed
                GSize = desired cell size in meters
                sort = bool if true sorts grid cells by z values... currently sorting not working..
        output: OP.Grid = 2D grid of all points in order x1,y1,z1,x2,y2,z2 etc... sorted or unsorted depending on value of 'sort'
                OP.Lasf = A pointer to the LasFile from which OP.Grid was created
*/
Las2DO Las2DGrid(LasFile * Lasf, double GSize, bool sort);


/* structure to be outout by Las2DZs
*/
struct Las2DZO{
	std::vector<std::vector<double> > Z;  // Z values
	std::vector<std::vector<size_t> > inds; // Indices
};

/* puts elevations from LasFile into a 2D vector 
        inputs: Lasf = a pointer to the LasFile to process
                Gsize = desired cell size in meters
                sort = bool, if true sorts cell Zs in each cell
                Ascend = bool if true sorts ascending (OP[0] = smallest), if false sorts descending (OP[0] = largest)
        outputs: a 2D vector that stores elevations of point in that grid cell (in meters)
        
                        NOTE: if sort == false, Ascend inconsequential
				if sort = true, indexes are not sorted
*/
Las2DZO Las2DZs(LasFile *Lasf, double GSize, bool sort, bool Ascend);


/* takes LasFile and returns maximum returns (UTM) in each grid cell. reprojected into 1D. Has a cleaning option to deal with "ghost" points
        inputs: Lasf = a pointer to the LasFile to be processed
                Gsize = grid size in meters
                clean = bool of whether to perform cleaning operation
                tol = the tolerance to be used for the cleaning operation in meters
        output: OP = vector of maximum returns in each cell in meters
                        NOTE = if clean == false, tol is inconsequential
                        NOTE = cleaning not very good.. needs work
*/
std::vector<double> Las2DZMax(LasFile *Lasf, double Gsize, bool clean, double tol);


/* takes LasFile and returns minimum returns (UTM) in each grid cell. reprojected into 1D. Has a cleaning option to deal with "ghost" points
        inputs: Lasf = a pointer to the LasFile to be processed
                Gsize = grid size in meters
                clean = bool of whether to perform cleaning operation
                tol = the tolerance to be used for the cleaning operation in meters
        output: OP = vector of minimum returns in each cell in meters
                        NOTE = if clean == false, tol is inconsequential
                        NOTE = cleaning not very good.. needs work
*/
std::vector<double> Las2DZMin(LasFile *Lasf, double Gsize, bool clean, double tol);


/* Performs Las2DZMax and returns a LasFile with the max points
        inputs: Lasf = pointer to LasFile to be performed on
                Gsize = grid size in meters
                clean = bool of whether to perform cleaning operation
                tol = the tolerance to be used for the cleaning operation in meters
        output: OP = vector of maximum returns in each cell in meters
                        NOTE = if clean == false, tol is inconsequential        
*/
LasFile Las2DZMaxLasOut(LasFile *Lasf, double Gsize, bool clean, double tol);


/* Takes vector of Max returns from Las2DZMax and convolves with a exponential smoothing kernel 
        inputs: Lasf = a pointer to the LasFile from which ZM came
                Gsize = the grid cell size used to produce ZM
                clean = bool for if cleaning operation should be perfomred
                tol = the tolerance to be used to clean (meters)
        outputs: OP = convolved grid values (2D projected to 1D)
*/
std::vector<double> ExpFilt(LasFile * Lasf, double Gsize, bool clean, double tol);


/* Takes a LasFile pointer and performs a gaussian convolution on the maximum returns in the grid
        inputs: Lasf = a pointer to the LasFile from which ZM came
                Gsize = the grid cell size used to produce ZM
                clean = bool for if cleaning operation should be perfomred
                tol = the tolerance to be used to clean (meters)
        outputs: OP = convolved grid values (2D projected to 1D)
*/
std::vector<double> LasMaxGausFilt(LasFile * Lasf, double Gsize, bool clean, double tol);


/* Takes LasFile pointer and performs a laplacian of gaussain (LoG) convolution on maximum grid returns
        inputs: Lasf = a pointer to the LasFile from which ZM came
                Gsize = the grid cell size used to produce ZM
                clean = bool for if cleaning operation should be perfomred tol = the tolerance to be used to clean (meters) outputs: OP = convolved grid values (2D projected to 1D)
*/
std::vector<double> LapFilt(LasFile * Lasf, double Gsize, bool clean, double tol);


/* Adjusts the coordinate minimum(s) LasFile generated by voxelization or griding that are going to be placed in a new grid or voxel grid. Need to do this to this to make center of grid cells correspond to where points have been placed during rasterization/voxelization.
     inputs: 
             Lasf = pointer to the LasFile being evaluated 
             GS = size of the grid cells being used
             x = bool.. if true -> adjust MinX
             y = bool.. if true -> adjust MinY
             z = bool.. if true -> adjust MinZ
            ~~~~~~~~~~~ NOTE!!! this modifies the min values in memory NOT in the actual file, if the file is written LasClean should be run to correctly assign the minimum values in the Las header.  ~~~~~~~~~~~~~~

*/
void VoxelGridMinAdjust(LasFile *Lasf, double GS, bool x, bool y, bool z);


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
void LasGridPad(LasFile *Lasf, double GS, bool adj, int32_t xpad, int32_t ypad, int32_t zpad);


/* Colors the points of a LasFile based on a classification grid (vector) (i.e., places points in a grid and colors them based on "some" metric)
     intputs:
          Lasf = a pointer to the LasFile to be colored
          CV = a pointer to the vector (grid) representing classsification.. needs to be a vector of doubles here
          GS = the grid size the classification grid was made using
     outputs: NA
*/
void LasColorByGridClass(LasFile *Lasf, std::vector<double> *CV, double GS);


/* Returns a vector of 2D grid positions (Projected in 1D) for each point in a LasFile
     inputs:
          Lasf = a pointer to the las file to be operated on
          GS = grid cell size
     
     outputs:
          OP = vector containing 2D grid positions for each point record in order they occur in the LasFile
*/
std::vector<int32_t> Las2DGridPositions(LasFile *Lasf, double GS);


/*  puts a scalar vector into the GPS time field in a LasFile

*/
void LasScalarForGPSTime(LasFile *Lasf, std::vector<double> *VIN, double GS);


/* Calculates the voxel volume in a given radius around supplied xy points 
   I 
 
 
 
*/ 
std::vector<uint32_t> PointVoxVolRadius(LasFile *Lasf, std::vector<double> * PV, double rad);

/* Voxelizes a las file and outputs the voxelized version
   inputs:
      LasIn = point to the LasFile to be voxelized
      VS = size of voxel to be voxelized with

   outputs: 
      LasOut = the voxelized point cloud      

*/

LasFile VoxelizeLas2Las(LasFile * LasIn, double VS);



/* structure for Las2LasZErrorClosestXY to output
*/
struct REO{
   LasFile LasF;
   std::vector<double> Error;
   std::vector<double> RefEle;
   std::vector<double> PCEle;
};

/*  Given a las file of reference points and an associated point cloud, finds the nearest point to all references points and calculates the z difference, also outputs a lasfile of the extracted closet points from the point cloud


*/
REO Las2LasZErrorClosestXY(LasFile * RP, LasFile * PC);

#endif
