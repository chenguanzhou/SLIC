#ifndef SLIC_H
#define SLIC_H

#include <iostream>
#include <vector>
#include <gdal_priv.h>
#include <gdal_alg.h>

typedef std::vector<double> FeatureVector;

class Slic
{
public:
    Slic(GDALDataset* &poSrcDS,GDALDataset* &poDstDS,int regionSize,double regularizer);
    bool StartSegmentation();

private:
    //Params
    GDALDataset*    _poSrcDS;  //Any Data Type
    GDALDataset*    _poDstDS;  //GDT_Int32
    int             _regionSize;
    double          _regularizer;

    //Temp var
    int             _width;
    int             _height;
    int             _bandCount;

    int             _M;         //SLIC starts by dividing the image domain into a regular grid with M×N tiles
    int             _N;    
    GDALDataType    _dataType;  //GDALDatatype
    int             _dataSize;  //Data size per pixel

    std::vector< FeatureVector  >       _centerVector;
    std::vector< std::vector<double> >  _normalizationParam;

    //Temp Funcs
    bool    _InitData();
    void    _InitCenterFeature (int i, int j, FeatureVector &featureVec);
    void    _GenerateSuperpixels();
    int     _GetNearestCenter(const std::vector<int> &candidateCenterID, const FeatureVector &featureVec);
    double  _ComputeDistance(const FeatureVector &vec1,const FeatureVector &vec2);
    void    _ComputeNewCenterVector();
    void    _ElininateSmallSuperpixel();
};

#endif //SLIC_H
