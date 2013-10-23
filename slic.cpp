#include "slic.h"



Slic::Slic(GDALDataset *&poSrcDS, GDALDataset *&poDstDS, int regionSize, int regularizer)
    :_regionSize(regionSize),_regularizer(regularizer)
{
    _poSrcDS = poSrcDS;
    _poDstDS = poDstDS;

    _width      = _poSrcDS->GetRasterXSize();
    _height     = _poSrcDS->GetRasterYSize();
    _bandCount  = _poSrcDS->GetRasterCount();
}

bool Slic::StartSegmentation()
{

}
