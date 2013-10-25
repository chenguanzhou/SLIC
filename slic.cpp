#include "slic.h"
#include <cmath>
#include <cstdio>
#include <assert.h>
#include <limits>
#include <map>
#include <time.h>
#include <omp.h>

typedef unsigned char uchar;


Slic::Slic(GDALDataset *&poSrcDS, GDALDataset *&poDstDS, int regionSize, double regularizer)
    :_regionSize(regionSize),_regularizer(regularizer)
{
    _poSrcDS = poSrcDS;
    _poDstDS = poDstDS;

    _width      = _poSrcDS->GetRasterXSize();
    _height     = _poSrcDS->GetRasterYSize();
    _bandCount  = _poSrcDS->GetRasterCount();

    _dataType   = _poSrcDS->GetRasterBand(1)->GetRasterDataType();
    _dataSize   = GDALGetDataTypeSize(_dataType)/8;
}

bool Slic::StartSegmentation()
{
    if (_InitData()==false)
        return false;

    _GenerateSuperpixels();
    _ElininateSmallSuperpixel();

    return true;
}

bool Slic::_InitData()
{
    // Check the params
    if (_poSrcDS == NULL|| _poDstDS == NULL)
    {
        std::cerr<<"Input image or output image invalid !"<<std::endl;
        return false;
    }

    if (_regionSize<0 || _regularizer<0)
    {
        std::cerr<<"Parameter regionSize and regularizer must bigger than 0!"<<std::endl;
        return false;
    }

    // Init some vars
    _M = static_cast<int> (static_cast<double>(_width)  / _regionSize + 0.5);
    _N = static_cast<int> (static_cast<double>(_height) / _regionSize + 0.5);

    // Init normalization params
    for (int k=0;k<_bandCount;++k)
    {
        GDALRasterBand* poBand =  _poSrcDS->GetRasterBand(k+1);
        int     bGotMin, bGotMax;
        double  adfMinMax[2];
        adfMinMax[0] = poBand->GetMinimum( &bGotMin );
        adfMinMax[1] = poBand->GetMaximum( &bGotMax );
        if( ! (bGotMin && bGotMax) )
            GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
        std::vector<double>  adfNormalizationParam(2);
        adfNormalizationParam[0] = 1./(adfMinMax[1]-adfMinMax[0]);
        adfNormalizationParam[1] = adfMinMax[0]/(adfMinMax[1]-adfMinMax[0]);
        _normalizationParam.push_back(adfNormalizationParam);
    }

    // Init centers
    for (int i=_regionSize/2;i<_height ;i+=_regionSize)
    {
        for (int j=_regionSize/2;j<_width ;j += _regionSize)
        {
            FeatureVector featureVec;
            _InitCenterFeature(i,j,featureVec);
            _centerVector.push_back(featureVec);
        }
    }

    return true;
}

void Slic::_InitCenterFeature(int i, int j, FeatureVector &featureVec)
{
    featureVec.clear();
    uchar *buffer = new uchar[5*5*_dataSize*_bandCount];
    _poSrcDS->RasterIO(GF_Read,j-2,i-2,5,5,buffer,5,5,_dataType,_bandCount,0,0,0,0);

    double minEdge = std::numeric_limits<double>::max();
    double minN;
    double minM;
    for (int n=1;n<4;++n)
    {
        for(int m=1;m<4;++m)
        {
            double edge = 0;
            edge += (buffer[(n-1)*5+m]-buffer[(n+1)*5+m])*(buffer[(n-1)*5+m]-buffer[(n+1)*5+m]);
            edge += (buffer[n*5+m-1]-buffer[n*5+m+1])*(buffer[n*5+m-1]-buffer[n*5+m+1]);

            if (edge < minEdge)
            {
                minEdge = edge;
                minN = n;
                minM = m;
            }
        }
    }

    int centerIndex = minN*5+minM;

    uchar* p = buffer;
    for(int k=0;k<_bandCount;++k,p += (5*5*_dataSize))
    {
        featureVec.push_back( SRCVAL(p,_dataType,centerIndex)/_regularizer);
    }
    featureVec.push_back(static_cast<double>(j+minM-2)/_regionSize);       //x
    featureVec.push_back(static_cast<double>(i+minN-2)/_regionSize);       //y
    delete []buffer;
}

static omp_lock_t lock;

void Slic::_GenerateSuperpixels()
{
    const int MAX_ITER = 10;
    for(int I=0;I<MAX_ITER;++I)
    {
        clock_t t1 = clock();
        std::cout<<"This is the "<<I<<"th circulation:"<<std::endl;

        omp_init_lock(&lock);

        #pragma omp parallel for
        for(int n=0;n<_N;++n)
        {
            for(int m=0;m<_M;++m)
            {
                // Init
                int nXOff  = m*_regionSize;
                int nYOff  = n*_regionSize;
                int nXSize = m==(_M-1)? _width  - m*_regionSize : _regionSize;
                int nYSize = n==(_N-1)? _height - n*_regionSize : _regionSize;

                uchar *bufferSrc = new uchar[nXSize*nYSize*_dataSize*_bandCount];
                int *bufferDst = new int[nXSize*nYSize];

                omp_set_lock(&lock);
                _poSrcDS->RasterIO(GF_Read,nXOff,nYOff,nXSize,nYSize,bufferSrc,nXSize,nYSize,_dataType,_bandCount,0,0,0,0);
                _poDstDS->GetRasterBand(1)->RasterIO(GF_Read,nXOff,nYOff,nXSize,nYSize,bufferDst,nXSize,nYSize,GDT_Int32,0,0);
                omp_unset_lock(&lock);

                std::vector< int > candidateCenterID;
                for(int i=-1;i<2;++i)
                {
                    for(int j=-1;j<2;++j)
                    {
                        if ((n+i)>=0 && (n+i)<_N && (m+j)>=0 && (m+j)<_M)
                            candidateCenterID.push_back( (n+i) * _M + (m+j) ) ;
                    }
                }

                // GetFeatureInfo
                FeatureVector featureVec(_bandCount + 2);
                for(int i=0,index=0;i<nYSize;++i)
                {
                    for(int j=0;j<nXSize;++j,++index)
                    {
                        uchar* p = bufferSrc;
                        for(int k=0;k<_bandCount;++k,p += nXSize*nYSize*_dataSize)
                        {
                            featureVec[k]= SRCVAL(p,_dataType,index)/_regularizer;
                        }
                        featureVec[_bandCount]   = static_cast<double>(nXOff+j)/_regionSize;  //x
                        featureVec[_bandCount+1] = static_cast<double>(nYOff+i)/_regionSize; //y
                        bufferDst[i*nXSize+j] = _GetNearestCenter(candidateCenterID,featureVec);
                    }
                }

                omp_set_lock(&lock);
                _poDstDS->GetRasterBand(1)->RasterIO(GF_Write,nXOff,nYOff,nXSize,nYSize,bufferDst,nXSize,nYSize,GDT_Int32,0,0);
                omp_unset_lock(&lock);

                delete []bufferSrc;
                delete []bufferDst;
            }

        }

        omp_destroy_lock(&lock);
        _ComputeNewCenterVector();
        std::cout<<"This circle cost: "<<(clock()-t1)/1000.<<"s"<<std::endl;
    }


}

int Slic::_GetNearestCenter(const std::vector< int > &candidateCenterID,const FeatureVector &featureVec)
{
    double minDis = std::numeric_limits<double>::max();
    int nMinID = 0;
    for (std::vector< int >::const_iterator iter = candidateCenterID.begin();iter!=candidateCenterID.end();++iter)
    {
        double dis = _ComputeDistance(_centerVector[*iter],featureVec);
        if ( dis < minDis)
        {
            minDis = dis;
            nMinID = *iter;
        }
    }
    return nMinID;
}

double Slic::_ComputeDistance(const FeatureVector &vec1, const FeatureVector &vec2)
{
    assert(vec1.size() == vec2.size());
    double dis = 0;
    for (unsigned int i=0;i<vec1.size();++i)
    {
        dis += (vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
    }
    return dis;
}

void Slic::_ComputeNewCenterVector()
{
    std::vector<int> centerCounts(_centerVector.size(),0);

    uchar *bufferSrc = new uchar[_width*_dataSize*_bandCount];
    int *bufferDst = new int[_width];

    for (unsigned int i=0;i<_centerVector.size();++i)
    {
        _centerVector[i].resize(_bandCount+2,0.);
    }

    for (int i=0;i<_height;++i)
    {
        _poSrcDS->RasterIO(GF_Read,0,i,_width,1,bufferSrc,_width,1,_dataType,_bandCount,0,0,0,0);
        _poDstDS->GetRasterBand(1)->RasterIO(GF_Read,0,i,_width,1,bufferDst,_width,1,GDT_Int32,0,0);

        for (int j=0;j<_width;++j)
        {
            uchar* p = bufferSrc;
            int nCenterID = bufferDst[j];
            centerCounts[nCenterID]++;
            for(int k=0;k<_bandCount;++k,p+=_width*_dataSize)
            {
                _centerVector[nCenterID][k] += SRCVAL(p,_dataType,j);
            }

            _centerVector[nCenterID][_bandCount]    += static_cast<double>(j);
            _centerVector[nCenterID][_bandCount+1]  += static_cast<double>(i);
        }
    }

    for (unsigned int i=0;i<_centerVector.size();++i)
    {
        if (centerCounts[i]==0)
            std::cerr<<"The "<<i<<"th superpixel has no pixels!"<<std::endl;
        for(int k=0;k<_bandCount;++k)
        {
            _centerVector[i][k] = _centerVector[i][k]/centerCounts[i]/_regularizer;
        }
        _centerVector[i][_bandCount]    /= centerCounts[i]*_regionSize;
        _centerVector[i][_bandCount+1]  /= centerCounts[i]*_regionSize;
    }

    delete []bufferSrc;
    delete []bufferDst;


}

void Slic::_ElininateSmallSuperpixel()
{
    int* bufferOldResult= new int[_width];
    int* bufferResult   = new int[_width];
    int* bufferUp       = new int[_width];
    int* bufferCurrent  = new int[_width];
    int* bufferDown     = new int[_width];
    GDALRasterBand* poBand = _poDstDS->GetRasterBand(1);

    poBand->RasterIO(GF_Read,0,0,_width,1,bufferUp,_width,1,GDT_Int32,0,0);
    poBand->RasterIO(GF_Read,0,1,_width,1,bufferCurrent,_width,1,GDT_Int32,0,0);

    for (int i=1;i<_height-1;++i)
    {
        poBand->RasterIO(GF_Read,0,i+1,_width,1,bufferDown,_width,1,GDT_Int32,0,0);
        bufferResult[0] = bufferCurrent[0];
        bufferResult[_width-1] = bufferCurrent[_width-1];
        for (int j=1;j<_width-1;++j)
        {
            std::map<int,int> counter;
            counter[bufferUp[j-1]]++;
            counter[bufferUp[j]]++;
            counter[bufferUp[j+1]]++;
            counter[bufferCurrent[j-1]]++;
            counter[bufferCurrent[j]]++;
            counter[bufferCurrent[j+1]]++;
            counter[bufferDown[j-1]]++;
            counter[bufferDown[j]]++;
            counter[bufferDown[j+1]]++;
            int maxCount = -1;
            int maxID = 0;
            for(std::map< int,int >::iterator iter = counter.begin();iter != counter.end();++iter)
            {
                if (iter->second>maxCount)
                {
                    maxCount    = iter->second;
                    maxID       = iter->first;
                }
            }
            //            if (counter[bufferCurrent[j]] == maxCount)
            //                bufferResult[j] = bufferCurrent[j];
            //            else
            bufferResult[j] = maxID;
        }

        if (i>=2)
            poBand->RasterIO(GF_Write,0,i-1,_width,1,bufferOldResult,_width,1,GDT_Int32,0,0);

        int* p          = bufferUp;
        bufferUp        = bufferCurrent;
        bufferCurrent   = bufferDown;
        bufferDown      = p;

        p               = bufferOldResult;
        bufferOldResult = bufferResult;
        bufferResult    = p;
    }

    poBand->RasterIO(GF_Write,0,_height-2,_width,1,bufferResult,_width,1,GDT_Int32,0,0);

    delete []bufferResult ;
    delete []bufferUp ;
    delete []bufferCurrent ;
    delete []bufferDown ;
    delete []bufferOldResult ;
}



