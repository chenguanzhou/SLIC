#include "slic.h"
#include <gdal_priv.h>
#include <iostream>
#include <time.h>

//G:\Document\data\clip-武大\江南集中区-后时相.tif G:\Document\data\clip-武大\后时相superpixel.tif 40 100
//G:\0804_1w_mosaic_218.img G:\0804.tif 40 100

int main(int argc , char** argv)
{
    if (argc<4)
    {
        std::cerr<<"Parameters is less than 4!"<<std::endl;
        return -1;
    }

    clock_t t1 = clock();

    GDALAllRegister();
    CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
    GDALDataset* poSrcDS = (GDALDataset*) GDALOpen(argv[1],GA_ReadOnly);
    if(poSrcDS == NULL)
    {
        std::cerr<<"Open File "<<argv[1]<<" Failed!"<<std::endl;
        return -1;
    }
    GDALDriver* poDriver = (GDALDriver*)GDALGetDriverByName("GTiff");
    GDALDataset* poDstDS =  poDriver->Create(argv[2],poSrcDS->GetRasterXSize(),poSrcDS->GetRasterYSize(),1,GDT_Int32,NULL);
    if (poDstDS == NULL)
    {
        std::cerr<<"Create File "<<argv[2]<<" Failed!"<<std::endl;
        return -1;
    }
    poDstDS->SetProjection(poSrcDS->GetProjectionRef());
    double adfGeoTransform[6];
    poSrcDS->GetGeoTransform(adfGeoTransform);
    poDstDS->SetGeoTransform(adfGeoTransform);

    Slic slic(poSrcDS,poDstDS,atoi(argv[3]),atof(argv[4]));
    slic.StartSegmentation();

    GDALClose(poSrcDS);
    GDALClose(poDstDS);


    std::cout<<"cost "<<(clock()-t1)/1000.<<"s"<<std::endl;
    return 0;
}
