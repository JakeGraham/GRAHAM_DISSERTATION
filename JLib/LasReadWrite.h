#ifndef LasReadWrite_h
#define LasReadWrite_h

#include<vector>
#include<iostream>
#include<fstream>
#include<stdint.h>
#include<cstring>
#include<istream>

// Structure with members frequently used, and raw data stored as a character vector
struct LasFile {
        std::vector<char> Data;
        uint32_t *OffsetToPoints, *NumberOfPointRecords;
        uint16_t *PointDataRecordLength;
        double *XScaleFactor, *YScaleFactor, *ZScaleFactor, *XOffset, *YOffset, *ZOffset, *MaxX, *MinX, *MaxY, *MinY, *MaxZ, *MinZ;
};

// "Cleans" LasFile.. sets min max values and rescales/offsets
void CleanLas(LasFile *Lasf);


// Takes a pointer to a LasFile and outputs as a Las 1.2
void WriteLas(LasFile *Lasf, bool clean, std::string Name);


// **************!!!!!!!!!!!!!!!        NEEDS WORK.. NUMBER OF POINT RECORDS ETC
// Formats pointers to member variables when LasFile is initialized internally and not read via I/O
void LasFormat(LasFile *LasOut);


// Reads a .las file and puts it in to a LasFile structure .. see above
LasFile ReadLas(char File[150]);


#endif
