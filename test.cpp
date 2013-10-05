#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "matrixops.h"
#include "node.h"
/*
#define DTau (2.5e-5)
#define DTSQUARED (6.25e-10)
#define DAMPING 1
#define STIFFNESS 100
*/
using namespace std;

//Global Variables
double DTau;
double TFinal;
double DTauSqd;
double TwiceDTau;
double *CurDisplacement;
double *LastDisplacement;
double *NextDisplacement;

struct ForceSinusoid
{
    public:
    double frequency;
    double amplitude;
};

struct Force
{
    public:
    int type; //0 for constant force, 1 for sinusoidal, 2 for ramp
    ForceSinusoid FS;
    double curval;
    double multiplier;
};

void ForcePrint(Force*array,int num)
{
    cout << endl;
    for(int i = 0;i<num;i++)
    {
        cout << array[i].curval << endl;
    }
}
//Select is used to choose double derivative, derivative, or proportional (2,1, or 0)
void selfrefAssemble(Node*list,double**m,int select, int NCNT, int size, int DOF)
{
    zero(m,size);
    for(int i=0;i<NCNT;i++)
    {
        if(!(list[i].getselffactor(select) <= 1e-5))
        {
            for(int j=0; j<DOF;j++)
            {
                m[i*DOF+j][i*DOF+j]=list[i].getselffactor(select);
            }
        }

    }
}

//Select is used to choose double derivative, derivative, or proportional (2,1, or 0)
//It is assumed that the dependant variable (ex. force) is only transferred axially
//transformation matrix is passed as a parameter so that it can be used multiple
//times without recalculating
void coupledAssemble(Node*list,double**c,int select, int NCNT, int DoF)
{
    zero(c,DoF*NCNT);
    double**transformation = CreateMatrix(DoF*2);
    for(int cnt=0;cnt<NCNT;cnt++)
    {
        for(int cnum=0;cnum<list[cnt].getconn();cnum++)
        {
            if(list[cnt].getconnfactor(cnum,select) >= 1e-5)
            {
                list[cnt].gettransformation(list[list[cnt].connto(cnum)], transformation);
                addLocToGlo(c,transformation,cnt,list[cnt].connto(cnum),DoF,list[cnt].getconnfactor(cnum,select));
            }
        }
    }
    DeleteMatrix(DoF*2,transformation);
}

// matrices named for M, C, and K for visualization purposes
void GAssemble(int size,double**M,double**C,double**K,double*G,Force*forces)
{
    double **temp = CreateMatrix(size), **temp2 = CreateMatrix(size);//may make this global so that memory is not constantly being allocated and deallocated
    //MPrint(temp,size,size);
    mulsca(M,temp,size,size,(2*DTauSqd));
    //MPrint(temp,size,size);
    subm(temp,K,temp,size,size);
    //MPrint(temp,size,size);
    mulsquvec(temp,CurDisplacement,G,size);
    //MPrint(temp,size,size);
    mulsca(M,temp,size,size,-(DTauSqd));
    //MPrint(temp,size,size);
    mulsca(C,temp2,size,size,(TwiceDTau));
    //MPrint(temp2,size,size);
    addm(temp,temp2,temp,size,size);
    //MPrint(temp,size,size);
    mulsquvec(temp,LastDisplacement,temp2[0],size);
    //MPrint(temp2,size,size);
    for(int i = 0;i<size;i++)
    {
        G[i] += forces[i].curval + temp2[0][i];
    }
    DeleteMatrix(size,temp);
    DeleteMatrix(size,temp2);
    temp = temp2 = NULL;
}
void AAssemble(int size,double**M,double**C,double**A)
{
    double **temp = CreateMatrix(size);
    mulsca(M,temp,size,size,DTauSqd);
    mulsca(C,A,size,size,TwiceDTau);
    addm(A,temp,A,size,size);
    DeleteMatrix(size,temp);
}


int main()
{
    //Initialization.
    int NCNT, DoF, numFixed = 0;
    int *nodesfixed;
    Node*nodes = NULL;
    Force*exforce = NULL;
    int size;
    ofstream node3 ("node3.txt");//used solely for problem 4. Could easily be generalized

    //set up timestep values
    cout << "This Program solves a system of 2nd order (no bending allowed) constant term ODEs using FEM" << endl;
    cout << "Enter the timestep in seconds." << endl;
    cin >> DTau;
    cout << "Enter last time to solve for." << endl;
    cin >> TFinal;
    int cycles = ceil(TFinal/DTau)+1;
    DTauSqd = 1/(DTau*DTau); //Precalculate repeatedly used values
    TwiceDTau = 1/(DTau*2);

    cout<<"Enter 1 to automatically set up the system in problem 4"<<endl;
    cin >> NCNT; //bad place to store this but only temporary

    if(NCNT == 1)
    {
        NCNT = 7;
        DoF = 2;
        size = 14;
        nodes = new Node[NCNT];
        for (int i=0;i<NCNT;i++)
        {
            nodes[i].initProb4(i);
        }
        exforce = new Force[14];
        for(int i=0;i<14;i++)
        {
            if(i == 12)
            {
                exforce[i].type = 1;
                exforce[i].FS.amplitude = 10;
                exforce[i].FS.frequency = 30;
                exforce[i].curval = 0;
            }
            else
            {
                exforce[i].type = 0;
                exforce[i].curval = 0;
            }
        }
        numFixed = 2;
        nodesfixed = new int[4];
        nodesfixed[0] = 0;
        nodesfixed[1] = 1;
        nodesfixed[2] = 8;
        nodesfixed[3] = 9;
    }

    else
    {
        cout << "enter # of nodes:          ";
        cin >> NCNT;
        cout << "Choose 1 or 2 degrees of freedom: ";
        cin >> DoF;
        nodes = new Node[NCNT];
        for (int i=0;i<NCNT;i++)
        {
            nodes[i].init(i,DoF);
        }
        size = NCNT*DoF;
        exforce = new Force[size];
        for(int i = 0;i<size;i++)
        {
            cout << "Input Type of Forced Dependant for Node" << i << endl << "0 for constant" << endl << "1 for sinusoidally varying" << endl << "2 for ramp" << endl;
            cin >> exforce[i].type;
            if(exforce[i].type == 0)
            {
                cout<< "Input Constant term" << endl;
            }
            else if(exforce[i].type == 1)
            {
                cout << "Input Amplitude" << endl;
                cin >> exforce[i].FS.amplitude;
                cout << "Input Frequency" << endl;
                cin >> exforce[i].FS.frequency;
            }
            else
            {
                cout << "Input Factor" <<endl;
                cin >> exforce[i].multiplier;
            }
        }

        //if time will reorganize this so it's not so inefficient
        for (int cnt = 0;cnt<NCNT;cnt++)
        {
            if(nodes[cnt].isFixed())
            {
                numFixed++;
            }
        }
        nodesfixed = new int[numFixed*2];
        for (int cnt = 0;cnt<NCNT;cnt++)
        {
            int cntr = 0;
            if(nodes[cnt].isFixed())
            {
                nodesfixed[cntr] = cnt*DoF;
                nodesfixed[cntr+1] = cnt*DoF + 1;
            }
        }
    }

    //allocate and zero out
    CurDisplacement = new double[size];
    LastDisplacement = new double[size];
    NextDisplacement = new double[size];
    double *G = new double[size];
    double **A = CreateMatrix(size);
    double **Ks = CreateMatrix(size);
    double **Cs = CreateMatrix(size);
    double **Ms = CreateMatrix(size);
    double **Kc = CreateMatrix(size);
    double **Cc = CreateMatrix(size);
    double **Mc = CreateMatrix(size);
    double updatevector[DoF];
    double **submatrix = CreateMatrix(size-(DoF*numFixed));
    double *subG = new double[size-(DoF*numFixed)];
    double *subnextdis = new double[size-(DoF*numFixed)];
    double *velocity = new double[size];
    double *acceleration = new double[size];

    for(int i=0;i<size;i++)
    {
        CurDisplacement[i] = 0;
        LastDisplacement[i] = 0;
        NextDisplacement[i] = 0;
        velocity[i] = 0;
        acceleration[i] = 0;
        exforce[i].curval = 0;
    }

    //first the selfreferential elements
    //they don't change
    selfrefAssemble(nodes,Ks,0,NCNT,size,DoF);
    selfrefAssemble(nodes,Cs,1,NCNT,size,DoF);
    selfrefAssemble(nodes,Ms,2,NCNT,size,DoF);

    int fixindx;
    int lj;
    for (int i = 0;i<cycles;i++)
    {
        //now the coupled elements
        coupledAssemble(nodes,Kc,0,NCNT,DoF);
        coupledAssemble(nodes,Cc,1,NCNT,DoF);
        coupledAssemble(nodes,Mc,2,NCNT,DoF);

        addm(Kc,Ks,Kc,size,size);
        addm(Cc,Cs,Cc,size,size);
        addm(Mc,Ms,Mc,size,size);

        //MPrint(Kc,size,size);
        //MPrint(Cc,size,size);
        //MPrint(Mc,size,size);
        //createsubmatrix(Cc, submatrix,size,nodesfixed);
        //MPrint(submatrix,(size-(DoF*numFixed)),(size-(DoF*numFixed)));
        //Update Force
        fixindx = 0;
        for(int j=0;j<size;j++)
        {
            if(exforce[j].type != 0) //if the force is external varying
            {
                if(exforce[j].type == 1) //sinusoid
                {
                    exforce[j].curval = (exforce[j].FS.amplitude)*sin((exforce[j].FS.frequency)*((i)*DTau));//+(mulonerow(Mc,acceleration,j,size)+(mulonerow(Cc,velocity,j,size))+(mulonerow(Kc,CurDisplacement,j,size)));
                }
                else //ramp
                {
                    exforce[j].curval = (exforce[j].multiplier)*((i)*DTau);//-(mulonerow(Mc,acceleration,j,size)-(mulonerow(Cc,velocity,j,size))-(mulonerow(Kc,CurDisplacement,j,size)));
                }
            }
            //else if (j == nodesfixed[fixindx])
            //
                //PrintV(CurDisplacement,size);
                //PrintV(velocity, size);
                //PrintV(acceleration, size);
                //MPrint(Kc,size,size);
                //MPrint(Cc,size,size);
                //MPrint(Mc,size,size);
                //exforce[j].curval = (mulonerow(Mc,acceleration,j,size)+(mulonerow(Cc,velocity,j,size))+(mulonerow(Kc,CurDisplacement,j,size)));
                //fixindx++;
            //}
        }
        ForcePrint(exforce,size);

        GAssemble(size,Mc,Cc,Kc,G,exforce);
        AAssemble(size,Mc,Cc,A);

        //PrintV(G,size);
        //Before solving, remove fixed nodes
        createsubmatrix(A, submatrix,size,nodesfixed);
        createsubG(G,subG,size,nodesfixed);

        PrintV(subG,(size-(numFixed*DoF)));

        //MPrint(A,size,size);
        //MPrint(submatrix,(size-(DoF*numFixed)),(size-(DoF*numFixed)));

        lud(submatrix,subG,(size-(DoF*numFixed)),subnextdis);

        //put next displacements into main vector
        fixindx = 0;
        lj = 0;

        //PrintV(subnextdis,(size-(DoF*numFixed)));
        for(int j=0;j<(size-(DoF*numFixed));j++)
        {
            while(nodesfixed[fixindx] == lj)
            {
                lj++;
                fixindx++;
                if(lj >= size)
                {
                    continue;
                }
            }
            NextDisplacement[lj] = subnextdis[j];
            lj++;
        }
        PrintV(NextDisplacement,size);

        //update velocity and acceleration
        for(int j=0;j<size;j++)
        {
            velocity[j] = (NextDisplacement[j]-LastDisplacement[j])*TwiceDTau;
            acceleration[j] = ((NextDisplacement[j]-2*CurDisplacement[j]+LastDisplacement[j])*DTauSqd);
        }
        if(i%100 == 0)
        {
            node3 << nodes[2].GetPosition(0) << " , " << nodes[2].GetPosition(1) << " , "
                << velocity[4] << " , " << velocity[5] << " , "
                << acceleration[4] << " , " << acceleration[5]
                << " , " << (i*DTau) << endl;
        }

        //Shift Indices
        for(int j=0;j<NCNT;j++)
        {
            for(int d=0;d<DoF;d++)
            {
                updatevector[d] = NextDisplacement[(DoF*j+d)];
                LastDisplacement[(j*DoF + d)] = CurDisplacement[(j*DoF + d)];
                CurDisplacement[(j*DoF + d)] = NextDisplacement[(j*DoF + d)];
            }
            nodes[j].UpdateDelta(updatevector);
        }


    }
    cout << "no seg faults" << endl;
    cin.get();
    return 0;
}
