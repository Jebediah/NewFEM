#define TRUE 1
#define FALSE 0

struct Connection
{
      public:
        int tonode;
        double factors[3]; //proportional, derivative and second derivative
};

class Node {
      private:
              bool *Fixed; //only boundary condition used in part 4
              double *initialstate;
              double *curstate, *paststate; //displacement or delta-charge, current and previous
              double *selffactors;
              int DegreesFreedom;
              int conncount;
              Connection *conn; //will hold which nodes are connected to each other
      public:
              void init(int,int); //this function will have the nodes either hard-coded or read from a CSV
              void initProb4(int);
              void UpdateDelta(double*);
              double GetCurState(int);
              double GetInitState(int);
              double GetPosition(int);
              bool isFixed();
              bool isFixed(int);
              int getconn();
              int connto(int);
              double get_2D_Distance(Node,int,int);
              double getFrac(Node,int,int);
              double getselffactor(int);
              double getconnfactor(int,int);
              void initProb3(int);
              void gettransformation(Node,double**);
};

/***********************CLASS FUNCTIONS***********************************/

void Node::init(int index, int DoF)
{
    DegreesFreedom = DoF;
    curstate = new double[DoF];
    paststate = new double[DoF];
    initialstate = new double[DoF];
    Fixed = new bool[DoF];
    selffactors = new double[3]; //proportional, derivative and second derivative term
    std::cout<< "Enter initial conditions for Node " << index << std::endl;
    for(int i = 0;i<DegreesFreedom;i++)
    {
        std::cout << "Initial Value (ex. charge in circuit, position in mds system) for Degree " << (i+1) << ": ";
        std::cin >> initialstate[i] ;
        paststate[i] = 0;
        curstate[i] = 0;
        std::cout<<std::endl<<"Is this degree fixed (1 if true, 0 if false)?";
        std::cin >> Fixed[i];
    }
    std::cout<<"Factor for self-referential proportional   ";
    std::cin>>selffactors[0];
    std::cout<<"Factor for self-referential derivative     ";
    std::cin>>selffactors[1];
    std::cout<<"Factor for self-referential 2nd derivative ";
    std::cin>>selffactors[2];
    std::cout<<"How many higher indexed nodes is this node connected to?"<<std::endl;
    std::cin>>conncount;
    conn = new Connection[conncount];
    for(int i = 0;i<conncount;i++)
    {
        std::cout<<"Input node index for connection"<<(i+1)<<std::endl;
        std::cin>>conn[i].tonode;
        std::cout<<"Factor for coupled proportional   ";
        std::cin>>conn[i].factors[0];
        std::cout<<"Factor for coupled derivative     ";
        std::cin>>conn[i].factors[1];
        std::cout<<"Factor for coupled 2nd derivative ";
        std::cin>>conn[i].factors[2];
    }
    std::cout<<std::endl;
}
void Node::initProb3(int index)
{
    DegreesFreedom = 1;
    curstate = new double[DegreesFreedom];
    paststate = new double[DegreesFreedom];
    initialstate = new double[DegreesFreedom];
    Fixed = new bool[DegreesFreedom];
    selffactors = new double[3];
    selffactors[0] = 0;
    selffactors[1] = 0;
    curstate[0] = 0;
    paststate[0] = 0;
    switch (index)
    {
         case 0:
            initialstate[0] = 80;
            conn = new Connection[1];
            conn[0].tonode = 1;
            conn[0].factors[0] = 1000; //stiffness k = 100N/m
            conn[0].factors[1] = 1200;   //damping   c = 1  Ns/m
            conn[0].factors[2] = 0;   //no coupled mass (whatever that is)
            selffactors[2] = 0;
            conncount = 1;
            Fixed[0] = TRUE;
            break;
         case 1:
            initialstate[0] = 60;
            conn = new Connection[2];
            conn[0].tonode = 2;
            conn[1].tonode = 3;
            conn[0].factors[0] = 1500; //stiffness k = 100N/m
            conn[0].factors[1] = 1600;   //damping   c = 1  Ns/m
            conn[0].factors[2] = 0;
            conn[1].factors[0] = 2000; //stiffness k = 100N/m
            conn[1].factors[1] = 1600;   //damping   c = 1  Ns/m
            conn[1].factors[2] = 0;
            selffactors[2] = 15;
            conncount = 2;
            Fixed[0] = FALSE;
            break;
         case 2:
            initialstate[0] = 40;
            conn = new Connection[1];
            conn[0].tonode = 4;
            conn[0].factors[0] = 1000; //stiffness k = 100N/m
            conn[0].factors[1] = 0;   //damping   c = 1  Ns/m
            conn[0].factors[2] = 0;
            selffactors[2] = 10;
            conncount = 1;
            Fixed[0] = FALSE;
            break;
         case 3:
            initialstate[0] = 40;
            conn = new Connection[1];
            conn[0].tonode = 4;
            conn[0].factors[0] = 2000; //stiffness k = 100N/m
            conn[0].factors[1] = 1600;   //damping   c = 1  Ns/m
            conn[0].factors[2] = 0;
            selffactors[2] = 0;
            conncount = 1;
            Fixed[0] = TRUE;
            break;
         case 4:
            initialstate[0] = 20;
            conn = new Connection[1];
            conn[0].tonode = 2;
            conn[0].factors[0] = 1200; //stiffness k = 100N/m
            conn[0].factors[1] = 2500;   //damping   c = 1  Ns/m
            conn[0].factors[2] = 0;
            selffactors[2] = 10;
            conncount = 1;
            Fixed[0] = FALSE;
            break;
         case 5:
            initialstate[0] = 0;
            conn = new Connection[0];
            selffactors[2] = 0;
            conncount = 0;
            Fixed[0] = TRUE;
            break;
         default:
            break;
    }
}

void Node::initProb4(int index)
{
    DegreesFreedom = 2;
    curstate = new double[DegreesFreedom];
    paststate = new double[DegreesFreedom];
    initialstate = new double[DegreesFreedom];
    Fixed = new bool[DegreesFreedom];
    selffactors = new double[3];
    selffactors[0] = 0;
    selffactors[1] = 0;
    selffactors[2] = 1; //They are have a mass of 1kg
    curstate[0] = 0;
    curstate[1] = 0;
    paststate[0] = 0;
    paststate[0] = 0;
    switch (index)
    {
         case 0:
            initialstate[0] = -200;
            initialstate[1] = 200;
            conn = new Connection[1];
            conn[0].tonode = 1;
            conn[0].factors[0] = 100; //stiffness k = 100N/m
            conn[0].factors[1] = 1;   //damping   c = 1  Ns/m
            conn[0].factors[2] = 0;   //no coupled mass (whatever that is)
            conncount = 1;
            Fixed[0] = TRUE;
            Fixed[1] = TRUE;
            break;
         case 1:
            initialstate[0] = -100;
            initialstate[1] = 100;
            conn = new Connection[1];
            conn[0].tonode = 2;
            conn[0].factors[0] = 100; //stiffness k = 100N/m
            conn[0].factors[1] = 1;   //damping   c = 1  Ns/m
            conn[0].factors[2] = 0;
            conncount = 1;
            Fixed[0] = FALSE;
            Fixed[1] = FALSE;
            break;
         case 2:
            initialstate[0] = 0;
            initialstate[1] = 0;
            conn = new Connection[2];
            conn[0].tonode = 3;
            conn[1].tonode = 5;
            conn[0].factors[0] = 100; //stiffness k = 100N/m
            conn[0].factors[1] = 1;   //damping   c = 1  Ns/m
            conn[0].factors[2] = 0;
            conn[1].factors[0] = 100; //stiffness k = 100N/m
            conn[1].factors[1] = 1;   //damping   c = 1  Ns/m
            conn[1].factors[2] = 0;
            conncount = 2;
            Fixed[0] = FALSE;
            Fixed[1] = FALSE;
            break;
         case 3:
            initialstate[0] = -100;
            initialstate[1] = -100;
            conn = new Connection[1];
            conn[0].tonode = 4;
            conn[0].factors[0] = 100; //stiffness k = 100N/m
            conn[0].factors[1] = 1;   //damping   c = 1  Ns/m
            conn[0].factors[2] = 0;
            conncount = 1;
            Fixed[0] = FALSE;
            Fixed[1] = FALSE;
            break;
         case 4:
            initialstate[0] = -200;
            initialstate[1] = -200;
            conn = NULL;
            conncount = 0;
            Fixed[0] = TRUE;
            Fixed[1] = TRUE;
            break;
         case 5:
            initialstate[0] = 100;
            initialstate[1] = 0;
            conn = new Connection[1];
            conn[0].tonode = 6;
            conn[0].factors[0] = 100; //stiffness k = 100N/m
            conn[0].factors[1] = 1;   //damping   c = 1  Ns/m
            conn[0].factors[2] = 0;   //no coupled mass (whatever that is)
            conncount = 1;
            Fixed[0] = FALSE;
            Fixed[1] = FALSE;
            break;
         case 6:
            initialstate[0] = 200;
            initialstate[1] = 0;
            conn = NULL;
            conncount = 0;
            Fixed[0] = FALSE;
            Fixed[1] = FALSE;
            break;
         default:
            break;
    }
}

void Node:: UpdateDelta(double *newPos)
{
     for(int i = 0;i<DegreesFreedom;i++)
     {
         paststate[i] =  curstate[i]; //pos[][0] is the old displacement
         curstate[i] = newPos[i];
     }
     return;
}

double Node::GetInitState(int index)
{
    return initialstate[index];
}

double Node::GetCurState(int index)
{
    return curstate[index];
}

double Node::GetPosition(int index)
{
    return curstate[index]+initialstate[index];
}
bool Node::isFixed(int index)
{
     return Fixed[index];
}
bool Node::isFixed()
{
    bool temp = TRUE;
    for(int i = 0;i<DegreesFreedom;i++)
    {
        temp = temp&Fixed[i];
    }
    return temp;
}
//to get distance if degrees describe physical orientation.
double Node::get_2D_Distance(Node other, int index1, int index2)
{
    double sum = 0;
    sum = (((*this).GetPosition(index1))-(other.GetPosition(index1)))*(((*this).GetPosition(index1))-(other.GetPosition(index1)));
    sum += (((*this).GetPosition(index2))-(other.GetPosition(index2)))*(((*this).GetPosition(index2))-(other.GetPosition(index2)));
    return sqrt(sum);
}

//this gets the fraction that one degree affects another
//If the system being solved was 2-D. X and Y would be index
//0 and 1 respectively. Cos is delta-X over length.
//This would be specified by calling Node1.getFrac(Node2,0,1)
//In a similar manner, sin would be  Node1.getFrac(Node2,1,0)
//main code needs to check degree of freedom to decide whether or not this
//function is necessary
double Node::getFrac(Node other, int index1, int index2)
{
    double delta = other.GetPosition(index1)-(*this).GetPosition(index1);
    double dis = (*this).get_2D_Distance(other,index1,index2);
    return delta/dis;
}

int Node::getconn()
{
    return conncount;
}

int Node::connto(int index)
{
    return conn[index].tonode;
}

double Node:: getconnfactor(int index,int select)
{
    return conn[index].factors[select];
}

double Node:: getselffactor(int index)
{
    return  selffactors[index];
}

void Node:: gettransformation(Node other, double**transformation)
{
    //1 degree of freedom. ex circuit, linear MDS system
    if(DegreesFreedom == 1)
    {
        if( initialstate[0] < other.GetInitState(0) )
        {
            for(int i=0;i<2;i++)
            {
                for(int j=0;j<2;j++)
                {
                    transformation[i][j] = -1 + 2*(i == j);
                }
            }
        }
        else
        {
            for(int i=0;i<2;i++)
            {
                for(int j=0;j<2;j++)
                {
                    transformation[i][j] = 1 - 2*(i == j);
                }
            }
        }
        return;
    }
    else //2 degrees of freedom
    {
        double a = (*this).getFrac(other,0,1);//cos
        double b = (*this).getFrac(other,1,0);//sin
        int temp;
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
                temp = i+j; //this variable cuts down on comparisons. Hopefully is faster though who knows.
                if ((i==0&&j==0)||(i==2&&j==2))
                   transformation[i][j] = a*a;
                else if ((i==1&&j==1)||(i==3&&j==3))
                   transformation[i][j] = b*b;
                else if (temp == 1 || temp == 5)
                   transformation[i][j] = a*b;
                else if (temp == 3)
                   transformation[i][j] = -a*b;
                else if (temp == 2)
                   transformation[i][j] = -a*a;
                else
                   transformation[i][j] = -b*b;
            }
        }
    }
    return;
}
