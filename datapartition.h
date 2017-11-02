#ifndef DATAPARTITION_H
#define DATAPARTITION_H

#include "Eigen\Dense"

using namespace Eigen;

typedef struct{
    VectorXd RefVal;
    MatrixXd Response;
    VectorXd Variable;
}spectraDataStruct;

typedef struct{
    VectorXi partitionRatio;
    spectraDataStruct CalData;
    spectraDataStruct ValidData;
    spectraDataStruct TestData;
}partitionedDataStruct;

partitionedDataStruct datapartition(spectraDataStruct OriginData,VectorXi partitionRatio);


#endif // DATAPARTITION_H
