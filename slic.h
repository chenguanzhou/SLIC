#ifndef SLIC_H
#define SLIC_H

#include <iostream>
#include <vector>
#include <gdal_priv.h>

typedef std::vector<double> FeatureVector;

class Slic
{
public:
    Slic(GDALDataset* &poSrcDS,GDALDataset* &poDstDS,int regionSize,int regularizer);
    bool StartSegmentation();

private:
    //Params
    GDALDataset*    _poSrcDS;  //Any Data Type
    GDALDataset*    _poDstDS;  //UInt32
    int             _regionSize;
    int             _regularizer;

    //Temp var
    int             _width;
    int             _height;
    int             _bandCount;
    int             _M;         //SLIC starts by dividing the image domain into a regular grid with M×N tiles
    int             _N;
    std::vector< std::vector< FeatureVector > > _centerVector;
};

#endif //SLIC_H
