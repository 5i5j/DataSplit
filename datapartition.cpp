#include "datapartition.h"
#include "pseudoinverse.h"
#include <vector>

using namespace std;

int whichMax(VectorXd ssy)
{
    int noMax;
    for(int i=0;i<ssy.size();i++)
    {
        if(ssy(i) == ssy.maxCoeff())
        {
            noMax = i;
        }
    }

    return noMax;
}

partitionedDataStruct datapartition(spectraDataStruct OriginData,VectorXi partitionRatio)
{
    vector<VectorXd> X;
    vector<double> Y;
    for(int i=0;i<OriginData.Response.rows();i++)
    {
        X.push_back(OriginData.Response.row(i));
        Y.push_back(OriginData.Response(i));
    }

    int nRows=OriginData.Response.rows();
    int nCols=OriginData.Response.cols();

    int CalSampNumb;
    int ValSampNumb;
    int TestSampNumb;

    CalSampNumb=floor(partitionRatio(0)/100*nRows);

    if(floor(partitionRatio(2)/100*nRows)>0)
    {
        ValSampNumb=floor(partitionRatio(1)/100*nRows);
        TestSampNumb=nRows-CalSampNumb-ValSampNumb;
    }
    else
    {
        ValSampNumb=nRows-CalSampNumb;
        TestSampNumb=0;
    }

    MatrixXd Xcal(CalSampNumb,nCols);
    VectorXd Ycal(CalSampNumb);
    MatrixXd Xval(ValSampNumb,nCols);
    VectorXd Yval(ValSampNumb);
    MatrixXd Xtest(ValSampNumb,nCols);
    VectorXd Ytest(ValSampNumb);

    MatrixXd Xmean=OriginData.Response.colwise().sum()/nRows;
    vector<VectorXd> SelectedX;
    MatrixXd Xsubspace=Xmean;

    vector<VectorXd>::iterator itX;
    vector<double>::iterator itY;

    for(int i=0;i<CalSampNumb;i++)
    {
        MatrixXd ResidueSpec;
        MatrixXd eyeI;
        eyeI.setIdentity(nCols,nCols);
        MatrixXd matX(X.size(),nCols);
        for(int j=0;j<X.size();j++)
        {
            matX.row(j)=X[j];
        }
        ResidueSpec=matX*(eyeI-pinv(Xsubspace)*Xsubspace);
        VectorXd SampResidue;
        SampResidue=(ResidueSpec*ResidueSpec.adjoint()).diagonal();
        int MinIndex=whichMax(SampResidue);

        SelectedX.push_back(matX.row(MinIndex));
        Xsubspace.resize(SelectedX.size(),nCols);
        for(int j=0;j<SelectedX.size();j++)
        {
            Xsubspace.row(j)=SelectedX[j];
        }
        Xcal.row(i)=X[MinIndex];
        Ycal(i)=Y[MinIndex];
        itX=X.begin();
        X.erase(itX+MinIndex);
        itY=Y.begin();
        Y.erase(itY+MinIndex);
    }

    MatrixXd tmpX(X.size(),nCols);
    for(int i=0;i<X.size();i++)
    {
        tmpX.row(i)=X[i];
    }

    /* ******************************** */
    Xmean=tmpX.colwise().sum()/X.size();
    SelectedX.clear();
    Xsubspace=Xmean;
    if(floor(partitionRatio(2)/100*nRows)>0)
    {
        for(int i=0;i<ValSampNumb;i++)
        {
            MatrixXd ResidueSpec;
            MatrixXd eyeI;
            eyeI.setIdentity(nCols,nCols);
            MatrixXd matX(X.size(),nCols);
            for(int j=0;j<X.size();j++)
            {
                matX.row(j)=X[j];
            }
            ResidueSpec=matX*(eyeI-pinv(Xsubspace)*Xsubspace);
            VectorXd SampResidue;
            SampResidue=(ResidueSpec*ResidueSpec.adjoint()).diagonal();
            int MinIndex=whichMax(SampResidue);
            SelectedX.push_back(matX.row(MinIndex));
            Xsubspace.resize(SelectedX.size(),nCols);
            for(int j=0;j<SelectedX.size();j++)
            {
                Xsubspace.row(j)=SelectedX[j];
            }
            Xval.row(i)=X[MinIndex];
            Yval(i)=Y[MinIndex];
            itX=X.begin();
            X.erase(itX+MinIndex);
            itY=Y.begin();
            Y.erase(itY+MinIndex);
        }

        for(int i=0;i<X.size();i++)
        {
            Xtest.row(i)=X[i];
            Ytest[i]=Y[i];
        }
    }
    else
    {
        for(int i=0;i<X.size();i++)
        {
            Xval.row(i)=X[i];
            Yval(i)=Y[i];
        }
    }

    partitionedDataStruct PartitionedData;
    PartitionedData.CalData.Variable=OriginData.Variable;
    PartitionedData.CalData.Response=Xcal;
    PartitionedData.CalData.RefVal=Ycal;

    PartitionedData.ValidData.Variable=OriginData.Variable;
    PartitionedData.ValidData.Response=Xval;
    PartitionedData.ValidData.RefVal=Yval;

    PartitionedData.TestData.Variable=OriginData.Variable;
    PartitionedData.TestData.Response=Xtest;
    PartitionedData.TestData.RefVal=Ytest;

    return PartitionedData;
}

