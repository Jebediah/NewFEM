#include <cstdlib>
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include "matrixops.h"
#include "node.h"

using namespace std;

int main()
{
    ofstream Data;
    Data.open ("Data.txt");
    
/*------------------Initialize Matrices, Nodes, and Variables-----------------*/
    
    double**M = CreateMatrix(6);                                                //Mass matrix
    double**C = CreateMatrix(6);                                                //Dampener matrix
    double**K = CreateMatrix(6);                                                //Spring matrix
    double**summat = CreateMatrix(6);                                           //Placeholder matrix for summing
    double**summat2 = CreateMatrix(3);                                          //Placeholder matrix for summing
    double F[6];                                                                //Force vector
    double G[6];                                                                //Gravitational force on all nodes
    double a[6];                                                                //Acceleration vector
    double v[6];                                                                //Velocity vector
    double uc[6];                                                               //Current displacement vector
    double up[6];                                                               //Previous displacement vector
    double un[6];                                                               //Next displacement vector
    double sumvec[6];                                                           //Placeholder vector for summing
    double sumvec2[3];                                                          //Placeholder vector for summing
    double sumvec3[3];
    double deltat = 0.0001;                                                     //Delta t
    double deltat2 = deltat*deltat;                                             //Delta t squared 
    Node n0;                                                                    //Node 0
    Node n1;                                                                    //Node 1
    Node n2;                                                                    //Node 2
    Node n3;                                                                    //Node 3
    Node n4;                                                                    //Node 4
    Node n5;                                                                    //Node 5
    
/*------------------Set Matrices, Call initProb3 on Nodes---------------------*/
    
    n0.initProb3(0);
    n1.initProb3(1);
    n2.initProb3(2);
    n3.initProb3(3);
    n4.initProb3(4);
    n5.initProb3(5);
    M[0][0] = n0.getselffactor(2);                                              //Set mass matrix        
    M[1][1] = n1.getselffactor(2);
    M[2][2] = n2.getselffactor(2);
    M[3][3] = n3.getselffactor(2);
    M[4][4] = n4.getselffactor(2);
    M[5][5] = n5.getselffactor(2);
    C[1][1] += n0.getconnfactor(0,1);
    C[1][0] -= n0.getconnfactor(0,1);
    C[0][1] -= n0.getconnfactor(0,1);
    C[1][1] += n1.getconnfactor(0,1);
    C[2][2] += n1.getconnfactor(0,1);
    C[2][1] -= n1.getconnfactor(0,1);
    C[1][2] -= n1.getconnfactor(0,1);
    C[1][1] += n1.getconnfactor(1,1);
    C[3][3] += n1.getconnfactor(1,1);
    C[3][1] -= n1.getconnfactor(1,1);
    C[1][3] -= n1.getconnfactor(1,1);
    C[2][2] += n2.getconnfactor(0,1);
    C[4][4] += n2.getconnfactor(0,1);
    C[4][2] -= n2.getconnfactor(0,1);
    C[2][4] -= n2.getconnfactor(0,1);
    C[3][3] += n3.getconnfactor(0,1);
    C[4][4] += n3.getconnfactor(0,1);
    C[4][3] -= n3.getconnfactor(0,1);
    C[3][4] -= n3.getconnfactor(0,1);
    C[4][4] += n4.getconnfactor(0,1);
    C[5][5] += n4.getconnfactor(0,1);
    C[5][4] -= n4.getconnfactor(0,1);
    C[4][5] -= n4.getconnfactor(0,1);
    K[0][0] += n0.getconnfactor(0,0);                                           //Set spring matrix
    K[1][1] += n0.getconnfactor(0,0);
    K[1][0] -= n0.getconnfactor(0,0);
    K[0][1] -= n0.getconnfactor(0,0);
    K[1][1] += n1.getconnfactor(0,0);
    K[2][2] += n1.getconnfactor(0,0);
    K[2][1] -= n1.getconnfactor(0,0);
    K[1][2] -= n1.getconnfactor(0,0);
    K[1][1] += n1.getconnfactor(1,0);
    K[3][3] += n1.getconnfactor(1,0);
    K[3][1] -= n1.getconnfactor(1,0);
    K[1][3] -= n1.getconnfactor(1,0);
    K[2][2] += n2.getconnfactor(0,0);
    K[4][4] += n2.getconnfactor(0,0);
    K[4][2] -= n2.getconnfactor(0,0);
    K[2][4] -= n2.getconnfactor(0,0);
    K[3][3] += n3.getconnfactor(0,0);
    K[4][4] += n3.getconnfactor(0,0);
    K[4][3] -= n3.getconnfactor(0,0);
    K[3][4] -= n3.getconnfactor(0,0);
    K[4][4] += n4.getconnfactor(0,0);
    K[5][5] += n4.getconnfactor(0,0);
    K[5][4] -= n4.getconnfactor(0,0);
    K[4][5] -= n4.getconnfactor(0,0);
    G[0] = n0.getselffactor(2)*(-9.80665);
    G[1] = n1.getselffactor(2)*(-9.80665);
    G[2] = n2.getselffactor(2)*(-9.80665);
    G[3] = n3.getselffactor(2)*(-9.80665);
    G[4] = n4.getselffactor(2)*(-9.80665);
    G[5] = n5.getselffactor(2)*(-9.80665);
    cout << "Mass matrix:" << endl;
    MPrint(M,6,6);
    cout << endl << "Dampener matrix:" << endl;
    MPrint(C,6,6);
    cout << endl << "Spring matrix:" << endl;
    MPrint(K,6,6);
    cout << endl;
    for (int i=0; i<6; i++)
    {
        uc[i] = 0;
    }
    for (int i=0; i<6; i++)
    {
        up[i] = 0;
    }
    for (int i=0; i<6; i++)
    {
        un[i] = 0;
    }
    
/*------------------Solution Start--------------------------------------------*/
    
    Data << setw(8) << "0";                                                     //Initial print to file
    for (int printcount=0; printcount<6; printcount++)                      
    {
        Data << setw(18) << uc[printcount];
    }
    Data << endl;
    
    for (double perm=0.0001; perm<=7.5; perm+=0.0001)                           //Permutations begin (0.0001s to 7.5s)
    {
        for (int i=0; i<6; i++)
        {
            F[i] = 0;
        }
        for (int i=0; i<6; i++)
        {
            F[i] += G[i];
        }
        
        subm2(K,M,summat,6,6,1,2/(deltat2));
        mulsquvec(summat,uc,sumvec,6);
        for (int i=0; i<6; i++)
        {
            F[i] -= sumvec[i];
        }
        subm2(M,C,summat,6,6,1/deltat2,(1/(2*deltat)));
        mulsquvec(summat,up,sumvec,6);
        for (int i=0; i<6; i++)
        {
            F[i] -= sumvec[i];
        }
        addm2(M,C,summat,6,6,1/deltat2,(1/(2*deltat)));
        
        sumvec2[0] = F[1];
        sumvec2[1] = F[2];
        sumvec2[2] = F[4];
        summat2[0][0] = summat[1][1];
        summat2[0][1] = summat[1][2];
        summat2[0][2] = summat[1][4];
        summat2[1][0] = summat[2][1];
        summat2[1][1] = summat[2][2];
        summat2[1][2] = summat[2][4];
        summat2[2][0] = summat[4][1];
        summat2[2][1] = summat[4][2];
        summat2[2][2] = summat[4][4];
        
        lud(summat2,sumvec2,3,sumvec3);
        
        un[1] = sumvec3[0];
        un[2] = sumvec3[1];
        un[4] = sumvec3[2];
        for (int i=0; i<6; i++)
        {
            up[i] = uc[i];
        }
        for (int i=0; i<6; i++)
        {
            uc[i] = un[i];
        }
        
        cout << endl << perm;
        Data << setw(8) << perm;                                                //Print to file displacements
        for (int printcount=0; printcount<6; printcount++)                      
        {
            Data << setw(18) << uc[printcount];
        }
        Data << endl;
    }
    system("pause");
}
