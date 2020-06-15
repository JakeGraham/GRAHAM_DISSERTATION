#include<vector>
#include<iostream>
#include<fstream>
#include<stdint.h>
#include<cstring>
#include<istream>
#include"LasReadWrite.h"

// "Cleans" LasFile.. sets min max values and rescales/offsets
void CleanLas(LasFile *Lasf){
        std::vector<int32_t> MinMax(6);
        for(size_t i = *Lasf->OffsetToPoints; i < Lasf->Data.size(); i += *Lasf->PointDataRecordLength){
                if(i == *Lasf->OffsetToPoints){
                        MinMax[0] = *(int32_t*) &Lasf->Data[i];
                        MinMax[1] = *(int32_t*) &Lasf->Data[i];
                        MinMax[2] = *(int32_t*) &Lasf->Data[i + 1];
                        MinMax[3] = *(int32_t*) &Lasf->Data[i + 1];
                        MinMax[4] = *(int32_t*) &Lasf->Data[i + 2];
                        MinMax[5] = *(int32_t*) &Lasf->Data[i + 2];
                } else {
                        if(MinMax[0] < *(int32_t*) &Lasf->Data[i]) {MinMax[0] = *(int32_t*) &Lasf->Data[i];}
                        if(MinMax[1] > *(int32_t*) &Lasf->Data[i]) {MinMax[1] = *(int32_t*) &Lasf->Data[i];}
                        if(MinMax[2] < *(int32_t*) &Lasf->Data[i + 1]) {MinMax[2] = *(int32_t*) &Lasf->Data[i + 1];}
                        if(MinMax[3] > *(int32_t*) &Lasf->Data[i + 1]) {MinMax[3] = *(int32_t*) &Lasf->Data[i + 1];}
                        if(MinMax[4] < *(int32_t*) &Lasf->Data[i + 2]) {MinMax[4] = *(int32_t*) &Lasf->Data[i + 2];}
                        if(MinMax[5] > *(int32_t*) &Lasf->Data[i + 2]) {MinMax[5] = *(int32_t*) &Lasf->Data[i + 2];}
                }
        }
}


// Takes a pointer to a LasFile and outputs as a Las 1.2
void WriteLas(LasFile *Lasf, bool clean, std::string Name){
        std::ofstream outfile;
        outfile.open(Name, std::ios::out | std::ios::binary);
        if(clean){
                CleanLas(Lasf);
        }
        outfile.write((char*)&Lasf->Data[0], Lasf->Data.size());
       outfile.close();
}




// **************!!!!!!!!!!!!!!!        NEEDS WORK.. NUMBER OF POINT RECORDS ETC
// Formats pointers to member variables when LasFile is initialized internally and not read via I/O
void LasFormat(LasFile *LasOut){
        // Get offsets and scaling factors & other useful header data..
        LasOut->OffsetToPoints = (uint32_t*) &LasOut->Data[96]; LasOut->PointDataRecordLength = (uint16_t*) &LasOut->Data[105];
        LasOut->NumberOfPointRecords = (uint32_t*) &LasOut->Data[107];
        LasOut->XScaleFactor = (double*) &LasOut->Data[131]; LasOut->YScaleFactor = (double*) &LasOut->Data[139]; LasOut->ZScaleFactor = (double*) &LasOut->Data[147];
        LasOut->XOffset = (double*) &LasOut->Data[155]; LasOut->YOffset = (double*) &LasOut->Data[163]; LasOut->ZOffset = (double*) &LasOut->Data[171];
        LasOut->MaxX = (double*) &LasOut->Data[179]; LasOut->MinX = (double*) &LasOut->Data[187]; LasOut->MaxY = (double*) &LasOut->Data[195];
        LasOut->MinY = (double*) &LasOut->Data[203]; LasOut->MaxZ = (double*) &LasOut->Data[211]; LasOut->MinZ = (double*) &LasOut->Data[219];
}


// Reads a .las file and puts it in to a LasFile structure .. see above
LasFile ReadLas(char File[150]) {
        LasFile LasOut;

        // read LAS file
        std::fstream file(File, std::ios::in | std::ios::out | std::ios::binary);
        uint32_t OStP, NoPR; uint16_t PRL;
        file.seekg(96); file.read((char*)&OStP, sizeof(uint32_t));
        file.seekg(105); file.read((char *)&PRL, sizeof(uint16_t));
        file.seekg(107); file.read((char*)&NoPR, sizeof(uint32_t));

        // Read file into memory
        uint64_t FileSize = OStP + (NoPR*PRL);
        LasOut.Data.resize(FileSize);
        file.seekg(0); file.read((char*)&LasOut.Data[0], FileSize);
        file.close();

        // Get offsets and scaling factors & other useful header data..
        LasOut.OffsetToPoints = (uint32_t*) &LasOut.Data[96]; LasOut.PointDataRecordLength = (uint16_t*) &LasOut.Data[105];
        LasOut.NumberOfPointRecords = (uint32_t*) &LasOut.Data[107];
        LasOut.XScaleFactor = (double*) &LasOut.Data[131]; LasOut.YScaleFactor = (double*) &LasOut.Data[139]; LasOut.ZScaleFactor = (double*) &LasOut.Data[147];
        LasOut.XOffset = (double*) &LasOut.Data[155]; LasOut.YOffset = (double*) &LasOut.Data[163]; LasOut.ZOffset = (double*) &LasOut.Data[171];
        LasOut.MaxX = (double*) &LasOut.Data[179]; LasOut.MinX = (double*) &LasOut.Data[187]; LasOut.MaxY = (double*) &LasOut.Data[195];
        LasOut.MinY = (double*) &LasOut.Data[203]; LasOut.MaxZ = (double*) &LasOut.Data[211]; LasOut.MinZ = (double*) &LasOut.Data[219];
        return LasOut;
}



