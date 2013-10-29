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
    for(int i=0;i<size;i++)
    {
        G[i] = forces[i].curval;
        for(int j=0;j<size;j++)
        {
            G[i] += ((2*M[i][j]/DTauSqd)-K[i][j])*CurDisplacement[j] + ((C[i][j]/TwiceDTau) - (M[i][j]/DTauSqd))*LastDisplacement[j];
        }
    }
}
void AAssemble(int size,double**M,double**C,double**A)
{
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            A[i][j] = (M[i][j]/DTauSqd+C[i][j]/TwiceDTau);
        }
    }
}


int main()
{
    //Initialization.
    int NCNT, DoF, numFixed = 0;
    int scratch;
    int *nodesfixed;
    Node*nodes = NULL;
    Force*exforce = NULL;
    bool prob3 = false, prob4 = false, freaksweep = false, grav = false;
    int size;
    ofstream output ("Output.txt");

    //set up timestep values
    cout << "This Program solves a system of 2nd order (no bending allowed) constant term ODEs using FEM" << endl;
    cout << "Enter the timestep in seconds." << endl;
    cin >> DTau;
    cout << "Enter last time to solve for." << endl;
    cin >> TFinal;
    cout << "Enter 1 for gravity, 0 for no gravity"
        << endl << "(Gravity is assumed to act in the negative y direction)" << endl;
    cin >> grav;

    cout<<"Enter 4 to automatically set up the system in problem 4"<<endl;
    cin >> scratch;

    if(scratch == 4)
    {
        prob4 = true;
        NCNT = 7;
        DoF = 2;
        size = 14;
        nodes = new Node[NCNT];
        for (int i=0;i<NCNT;i++)
        {
            nodes[i].initProb4(i);
        }
        exforce = new Force[14];
        cout << "Enter 1 for sin force, 2 for ramp," << endl << "and 3 to perform a frequency sweep" << endl;
        cin >> scratch;
        for(int i=0;i<14;i++)
        {
            exforce[i].type = 0;
            exforce[i].curval = 0;
        }
        if(scratch == 1)
        {
            exforce[12].type = 1;
            exforce[12].FS.amplitude = 10;
            cout << "Input Frequency" << endl;
            cin >> exforce[12].FS.frequency;
        }
        else if(scratch == 2)
        {
            exforce[12].type = 2;
            exforce[12].multiplier = 20;
        }
        else if(scratch == 3)
        {
            freaksweep = true;
            exforce[12].type = 1;
            exforce[12].FS.amplitude = 10;
        }
        numFixed = 2;
        nodesfixed = new int[4];
        nodesfixed[0] = 0;
        nodesfixed[1] = 1;
        nodesfixed[2] = 8;
        nodesfixed[3] = 9;
        scratch = 4;
    }
    else
    {
        cout << "enter # of nodes:          ";
        cin >> NCNT;
        DoF = 2;
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
    double *amplitude = new double[size];
    double *mean = new double[size];
    double *peak = new double[size];

    for(int i=0;i<size;i++)
    {
        CurDisplacement[i] = 0;
        LastDisplacement[i] = 0;
        NextDisplacement[i] = 0;
        velocity[i] = 0;
        acceleration[i] = 0;
        amplitude[i] = 0;
        mean[i] = 0;
        peak[i] = -1e20;
        if(grav && ((i&1) == 1))
        {
            exforce[i].curval -= 9.80665*nodes[(i/DoF)].getselffactor(2); //exploiting integer division here
        }
    }

    //first the selfreferential elements
    //they don't change
    selfrefAssemble(nodes,Ks,0,NCNT,size,DoF);
    selfrefAssemble(nodes,Cs,1,NCNT,size,DoF);
    selfrefAssemble(nodes,Ms,2,NCNT,size,DoF);

    int fixindx;
    int lj;

    //set up time
    int cycles = ceil(TFinal/DTau)+1;
    int plotinterval = cycles/5000; //plotting millions of points is insane
    if(plotinterval<1)
    {
        plotinterval = 1;
    }
    DTauSqd = DTau*DTau; //Precalculate repeatedly used values
    TwiceDTau = DTau*2;
    double T1;

    //calculate the maximum value that may be used to get a frequency amplitude
    if(prob4)
    {
        if(exforce[12].type == 1)
        {
            T1 = TFinal - (2*pi)/exforce[12].FS.frequency;
        }
    }

    if(!freaksweep)
    {
        for (int i = 0;i<cycles;i++)
        {
            //now the coupled elements
            coupledAssemble(nodes,Kc,0,NCNT,DoF);
            coupledAssemble(nodes,Cc,1,NCNT,DoF);
            coupledAssemble(nodes,Mc,2,NCNT,DoF);

            addm(Kc,Ks,Kc,size,size);
            addm(Cc,Cs,Cc,size,size);
            addm(Mc,Ms,Mc,size,size);

            fixindx = 0;
            for(int j=0;j<size;j++)
            {
                if(exforce[j].type != 0) //if the force is external varying
                {
                    if(exforce[j].type == 1) //sinusoid
                    {
                        if(j&1)
                        {
                            exforce[j].curval = (exforce[j].FS.amplitude)*sin((exforce[j].FS.frequency)*((i)*DTau)) - grav*(9.80665*nodes[(j/DoF)].getselffactor(2));
                        }
                        else
                        {
                            exforce[j].curval = (exforce[j].FS.amplitude)*sin((exforce[j].FS.frequency)*((i)*DTau));
                        }
                    }
                    else //ramp
                    {
                        if(j&1)
                        {
                            exforce[j].curval = (exforce[j].multiplier)*((i)*DTau) - grav*(9.80665*nodes[(j/DoF)].getselffactor(2));
                        }
                        else
                        {
                            exforce[j].curval = (exforce[j].multiplier)*((i)*DTau);
                        }
                    }
                }
            }
            GAssemble(size,Mc,Cc,Kc,G,exforce);
            AAssemble(size,Mc,Cc,A);
            createsubmatrix(A, submatrix,size,nodesfixed);
            createsubG(G,subG,size,nodesfixed);
            lud(submatrix,subG,(size-(DoF*numFixed)),subnextdis);

            //put next displacements into main vector
            fixindx = 0;
            lj = 0;

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
            //update velocity and acceleration
            for(int j=0;j<size;j++)
            {
                velocity[j] = (NextDisplacement[j]-LastDisplacement[j])/TwiceDTau;
                acceleration[j] = ((NextDisplacement[j]-2*CurDisplacement[j]+LastDisplacement[j])/DTauSqd);
            }
            //check for peak and sum for mean if needed
            if(prob4)
            {
                if((exforce[12].type == 1) && (i*DTau >= T1))
                {
                    for(int k = 0;k<size;k++)
                    {
                        if(peak[k] < CurDisplacement[k])
                        {
                            peak[k] = CurDisplacement[k];
                        }
                    }
                    freqsweep(CurDisplacement,LastDisplacement,mean,DTau,exforce[12].FS.frequency,size);
                }
            }
            if((i%plotinterval) == 0)
            {
                if(prob4)
                {
                    if(exforce[12].type == 1)
                    {
                        if(i)
                        output << CurDisplacement[4] << " , " << CurDisplacement[5] << " , "
                            << velocity[4] << " , " << velocity[5] << " , "
                            << acceleration[4] << " , " << acceleration[5]
                            << " , " << (i*DTau) << endl;
                    }
                    else
                    {
                        for(int p = 0;p<NCNT;p++)
                        {
                            for(int d = 0; d<DoF;d++)
                            {
                                output << CurDisplacement[p*DoF+d] << " , " ;
                            }
                            for(int d = 0; d<DoF;d++)
                            {
                                output << velocity[p*DoF+d] << " , ";
                            }
                            for(int d = 0; d<DoF;d++)
                            {
                                output << acceleration[p*DoF+d] << " , ";
                            }
                        }
                        output <<  (i*DTau) << endl;
                    }
                }
                else if (prob3)
                {
                    output << CurDisplacement[1] << " , "<< CurDisplacement[2]
                    << " , " << CurDisplacement[4] << " , " << (i*DTau) << endl;
                }
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
        for(int i=0; i<NCNT;i++)
        {
            for(int d=0;d<DoF;d++)
            {
                amplitude[i*DoF+d] = abs(peak[i*DoF+d]-mean[i*DoF+d]);
                cout << "Amplitude for Node " << i << ", dim " << d << " is " << amplitude[i*DoF+d] << endl
                    << "Mean for Node " << i << ", dim " << d << " is " << mean[i*DoF+d] << endl
                    << "Peak for Node " << i << ", dim " << d << " is " << peak[i*DoF+d] << endl << endl;
            }
        }
    }
    else
    {
        for(exforce[12].FS.frequency = 1;exforce[12].FS.frequency < 150.1;exforce[12].FS.frequency += 0.50)
        {
            TFinal = 25 + (2*pi)/exforce[12].FS.frequency;
            cycles = ceil(TFinal/DTau)+1;
            DTauSqd = DTau*DTau; //Precalculate repeatedly used values
            TwiceDTau = DTau*2;
            for (int i = 0;i<cycles;i++)
            {
                //now the coupled elements
                coupledAssemble(nodes,Kc,0,NCNT,DoF);
                coupledAssemble(nodes,Cc,1,NCNT,DoF);
                coupledAssemble(nodes,Mc,2,NCNT,DoF);

                addm(Kc,Ks,Kc,size,size);
                addm(Cc,Cs,Cc,size,size);
                addm(Mc,Ms,Mc,size,size);

                fixindx = 0;
                for(int j=0;j<size;j++)
                {
                    if(exforce[j].type != 0) //if the force is external varying
                    {
                        if(exforce[j].type == 1) //sinusoid
                        {
                            if(j&1)
                            {
                                exforce[j].curval = (exforce[j].FS.amplitude)*sin((exforce[j].FS.frequency)*((i)*DTau)) - grav*(9.80665*nodes[(j/DoF)].getselffactor(2));
                            }
                            else
                            {
                                exforce[j].curval = (exforce[j].FS.amplitude)*sin((exforce[j].FS.frequency)*((i)*DTau));
                            }
                        }
                        else //ramp
                        {
                            if(j&1)
                            {
                                exforce[j].curval = (exforce[j].multiplier)*((i)*DTau) - grav*(9.80665*nodes[(j/DoF)].getselffactor(2));
                            }
                            else
                            {
                                exforce[j].curval = (exforce[j].multiplier)*((i)*DTau);
                            }
                        }
                    }
                }
                GAssemble(size,Mc,Cc,Kc,G,exforce);
                AAssemble(size,Mc,Cc,A);
                createsubmatrix(A, submatrix,size,nodesfixed);
                createsubG(G,subG,size,nodesfixed);
                lud(submatrix,subG,(size-(DoF*numFixed)),subnextdis);
                //put next displacements into main vector
                fixindx = 0;
                lj = 0;

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
                //update velocity and acceleration
                for(int j=0;j<size;j++)
                {
                    velocity[j] = (NextDisplacement[j]-LastDisplacement[j])/TwiceDTau;
                    acceleration[j] = ((NextDisplacement[j]-2*CurDisplacement[j]+LastDisplacement[j])/DTauSqd);
                }
                if(i*DTau >= 25)
                {
                    for(int k = 0;k<size;k++)
                    {
                        if(peak[k] < CurDisplacement[k])
                        {
                            peak[k] = CurDisplacement[k];
                        }
                    }
                    freqsweep(CurDisplacement,LastDisplacement,mean,DTau,exforce[12].FS.frequency,size);
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
            for(int i=0; i<NCNT;i++)
            {
                for(int d=0;d<DoF;d++)
                {
                    amplitude[i*DoF+d] = abs(peak[i*DoF+d]-mean[i*DoF+d]);
                    output<<  amplitude[i*DoF+d] << " , ";
                    mean[i*DoF+d] = 0;
                    peak[i*DoF+d] = -1e20;
                    CurDisplacement[i*DoF+d] = 0;
                    LastDisplacement[i*DoF+d] = 0;
                    NextDisplacement[i*DoF+d] = 0;
                }
            }
            output << exforce[12].FS.frequency << endl;
        }
    }

    cout << "no seg faults" << endl;
    cin.get();
    return 0;
}
