#include <iostream>
#include <fstream>
#include <iostream>
#include <vector>
//#include "datapartition.h"

#include "Eigen\Dense"

using namespace Eigen;

using namespace std;

int main(int argc, char *argv[])
{
//    MatrixXd X,Y;

//    X.resize(150,4);

//    ifstream inputTD("TrainingData.txt");
//    for (int i = 0; i < 150; i++)
//    {
//        for (int j = 0; j < 4; j++)
//        {
//            inputTD >> X(i,j) ;
//        }
//    }
//    inputTD.close();

//    cout << X(0,0) << endl;

    VectorXi v(2);
    v(0)=1;
    v(1)=2;

    VectorXd V1(3),V2(3),V3(3);
    V1(0)=1;
    V1(1)=2;
    V1(2)=3;
    V2(0)=4;
    V2(1)=5;
    V2(2)=6;
    V3(0)=7;
    V3(1)=8;
    V3(2)=9;

    vector<VectorXd> X;
    X.push_back(V1);
    X.push_back(V2);
    X.push_back(V3);

//    cout<<X[0]<<endl;
//    cout<<X[1]<<endl;
//    cout<<X[2]<<endl;

    vector<VectorXd>::iterator it=X.begin()+1;

//    X.erase(it);

    cout<<X.size()<<endl;

    return 0;
}
