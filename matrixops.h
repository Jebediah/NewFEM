/************************MATRIX FUNCTIONS************************************/

//creates a square matrix of size n
double **CreateMatrix(int size)
{
    double **matrix = new double*[size];
    for(int i=0;i<size;i++)
    {
        matrix[i] = new double [size];
        for(int j=0;j<size;j++)
        {
            matrix[i][j] = 0;
        }
    }
    return matrix;
}

//Delete a matrix given number of rows
void DeleteMatrix(int rows, double**Matrix)
{
    for(int i=0;i<rows;i++)
    {
        delete [] Matrix[i];
    }
    delete [] Matrix;
    Matrix = NULL;
}

//prints out a matrix
void MPrint(double**matrix,int rows, int columns)
{
    std::cout << std::endl;
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<columns;j++)
        {
            std::cout.width(5);
            std::cout << matrix[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//given costheta and sintheta, this calculates the transformation matrix
//assumes the pointer given to it points to a matrix size 4x4
void trans(double a, double b, double** trans)
{
    double ab = a*b, asquared = a*a, bsquared = b*b; //trading memory for less cycles because I have gigabytes of space
    int temp;
    for (int i = 0; i<4; i++)
    {
        for (int j = 0; j<4; j++)
        { //time for some if statements
            temp = i+j; //this variable cuts down on comparisons. Hopefully is faster though who knows.
            if ((i==0&&j==0)||(i==2&&j==2))
               trans[i][j] = asquared;
            else if ((i==1&&j==1)||(i==3&&j==3))
               trans[i][j] = bsquared;
            else if (temp == 1 || temp == 5)
               trans[i][j] = ab;
            else if (temp == 3)
               trans[i][j] = -ab;
            else if (temp == 2)
               trans[i][j] = -asquared;
            else
               trans[i][j] = -bsquared;
        }
    }
    return;
}

//given two square matrices, this function returns the multiplication
//it assumes that the matrices are square and that it has been handed a
//correct pointer
void mulsquare(double**a,double**b,double**result,int size)
{
     for (int i = 0; i<size;i++)
     {
         for (int j = 0; j<size;j++)
         {
             result[i][j] = 0; //zero the result matrix out
             for (int k = 0; k<size;k++)
             {
                 result[i][j] += ((a[i][k])*(b[k][j]));
             }
         }
     }
     return;
}

//multiplies a square matrix of size nxn and a vector of length n
void mulsquvec(double**square, double*vector, double*result, int size)
{
     for (int i = 0; i<size; i++)
     {
         result[i] = 0;
         for (int j = 0; j<size;j++)
         {
             result[i] += square[i][j]*vector[j] ;
         }
     }
}

void zero(double**matrix,int size)
{
    for (int i = 0; i<size;i++)
     {
         for (int j = 0; j<size;j++)
         {
             matrix[i][j] = 0;
         }
     }
     return;
}

//multiplies a matrix by a scalar
void mulsca(double**matrix, double**out, int rows, int col, double scalar)
{
     for (int i=0;i<rows;i++)
     {
         for (int j=0;j<col;j++)
         {
             out[i][j] = scalar*matrix[i][j];
         }
     }
}

//adds two matrices of nxm size

void addm(double**a,double**b,double**result,int rows,int col)
{
     for (int i=0;i<rows;i++)
     {
         for (int j=0;j<col;j++)
         {
             result[i][j] = a[i][j] +b[i][j];
         }
     }
}
//General Form
void addLocToGlo(double**global,double**local,int GIndex1, int GIndex2, int DoF)
{
    GIndex1 *= DoF;
    GIndex2 *= DoF;
    for (int i=0;i<2*DoF;i++)
    {
        for (int j=0;j<2*DoF;j++)
        {
            if(i<DoF)
            {
                if(j<DoF)
                    global[GIndex1+i][GIndex1+j] += local[i][j];
                else
                    global[GIndex1+i][GIndex2+(j-DoF)] += local[i][j];
            }
            else if(j<DoF)
                global[GIndex2+(i-DoF)][GIndex1+j] += local[i][j];
            else
                global[GIndex2+(i-DoF)][GIndex2+(j-DoF)] += local[i][j];
        }
    }
    return;
}
//Multiplies the local matrix by a factor before adding it to the global matrix
void addLocToGlo(double**global,double**local,int GIndex1, int GIndex2, int DoF, double factor)
{
    GIndex1 *= DoF;
    GIndex2 *= DoF;
    for (int i=0;i<2*DoF;i++)
    {
        for (int j=0;j<2*DoF;j++)
        {
            if(i<DoF)
            {
                if(j<DoF)
                    global[GIndex1+i][GIndex1+j] += factor*local[i][j];
                else
                    global[GIndex1+i][GIndex2+(j-DoF)] += factor*local[i][j];
            }
            else if(j<DoF)
                global[GIndex2+(i-DoF)][GIndex1+j] += factor*local[i][j];
            else
                global[GIndex2+(i-DoF)][GIndex2+(j-DoF)] += factor*local[i][j];
        }
    }
    return;
}

void subm(double**a,double**b,double**result,int rows,int col)
{
     for (int i=0;i<rows;i++)
     {
         for (int j=0;j<col;j++)
         {
                result[i][j] = a[i][j] - b[i][j];
         }
     }
}

//changes a square matix into its transpose
void transpose(double**matrix,int size)
{
     double temp;
     for (int i=0;i<size;i++)
     {
         for (int j=(size-1);j>i;j--)
         {
             temp = matrix[i][j];
             matrix[i][j] = matrix[j][i];
             matrix[j][i] = temp;
         }
     }
}

double det(double**matrix, int size)
{
    double sum = 0;
    if (size > 2)
    {
        for (int root=0; root<size; root++)
        {
            double **submatrix = CreateMatrix(size-1);
            int cnti=0;
            int cntj=0;
            for (int j=0; j<size; j++)
            {
                if (j != root)
                {
                    for (int i=1; i<size; i++)
                    {
                        submatrix[cnti][cntj] = matrix[i][j];
                        cnti++;
                    }
                    cnti = 0;
                    cntj++;
                }
            }
            if ((root & 1) == 0)
            {
                sum += (matrix[0][root])*(det(submatrix, size-1));
                DeleteMatrix((size-1),submatrix);
                submatrix = NULL;
            }
            else
            {
                sum -= (matrix[0][root])*(det(submatrix, size-1));
                DeleteMatrix((size-1),submatrix);
                submatrix = NULL;
            }
        }
    }
    else if (size == 2)
    {
        sum += ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
    }
    else
    {
        sum += matrix[0][0];
    }
    return sum;
}


void EditMatrix(double**matrix, int rows, int col)
{
    std::cout << "Enter values for matrix" << std::endl;
    for (int i=0;i<rows;i++)
    {
         for (int j=0;j<col;j++)
         {
             std::cout << "(" << i << "," << j << ") = ";
             std::cin >> matrix[i][j];
             std::cout << std::endl;
         }
    }
}

//calculates the cofactor of a matrix
double**cofact(double**matrix,int size)
{
    double**temp = CreateMatrix(size);
    if (size > 2)
    {
        int cnti = 0, cntj = 0;
        for (int i=0;i<size;i++)
        {
            for (int j=0;j<size;j++)
            {
                double **submatrix = CreateMatrix(size-1);
                for(int ii = 0;ii<size;ii++)
                {
                    if(ii!=i)
                    {
                        for(int jj=0;jj<size;jj++)
                        {
                            if(jj!=j)
                            {
                            submatrix[cnti][cntj] = matrix[ii][jj];
                            cntj++;
                            }
                        }
                        cntj = 0;
                        cnti++;
                    }
                }
                cnti = 0;
                if (((i+j) & 1) == 0)
                {
                    temp[i][j] = (det(submatrix, size-1));
                }
                else
                {
                    temp[i][j] = -(det(submatrix, size-1));
                }
                DeleteMatrix((size-1),submatrix);
                submatrix = NULL;
            }
        }
    }
    else if(size==2)
    {
        temp[0][0] = matrix[1][1];
        temp[1][0] = -matrix[0][1];
        temp[0][1] = -matrix[1][0];
        temp[1][1] = matrix[0][0];
    }
    else
    {
        temp[0][0] = 1;
    }
    return temp;
}
void Copy(double**a,double**b,int size)
{
    for (int i=0;i<size;i++)
    {
        for (int j=0;j<size;j++)
        {
            a[i][j] = b[i][j];
        }
    }
}
//returns the inverse of a matrix. Looking back, writing this and the cofactor matrix was a waste of time.
void Invert(double**matrix, int size)
{
    double determinant = det(matrix,size);
    if(determinant == 0)
    {
        std::cout<< "bad assumptions"<<std::endl;
        return;
    }
    double**cofactor = cofact(matrix,size);
    transpose(cofactor,size);
    mulsca(cofactor,cofactor,size,size,(1/determinant));
    Copy(matrix,cofactor,size);
    DeleteMatrix(size,cofactor);
    cofactor = NULL;
    return;
}

void lud(double**matrix, double*b, int size, double*x)
{
     double**upper = CreateMatrix(size);
     double**lower = CreateMatrix(size);
     double y[size];
     double sum;
     for (int i=0; i<size; i++)
     {
         upper[i][i] = 1;
         lower[i][0] = matrix[i][0];
     }

     for (int perm=0; perm<((size-1)*2); perm++)
     {
         if ((perm & 1) == 0)
         {
                   for (int i=((perm/2)+1); i<size; i++)
                   {
                       sum = matrix[perm/2][i];
                       for (int j=0; j<(perm/2); j++)
                       {
                           sum -= lower[perm/2][j]*upper[j][i];
                       }
                       sum = sum/lower[perm/2][perm/2];
                       upper[perm/2][i] = sum;
                   }
         }
         else
         {
                   for (int i=((perm+1)/2); i<size; i++)
                   {
                       sum = matrix[i][(perm+1)/2];
                       for (int j=0; j<((perm+1)/2); j++)
                       {
                           sum -= lower[i][j]*upper[j][(perm+1)/2];
                       }
                       lower[i][(perm+1)/2] = sum;
                   }
         }
     }
     for(int i=0; i<size; i++)
     {
         sum = b[i];
         for(int j=0; j<i; j++)
         {
                 sum -= lower[i][j]*y[j];
         }
         sum = sum/lower[i][i];
         y[i] = sum;
     }
     for(int i=(size-1); i>-1; i--)
     {
         sum = y[i];
         for(int j=(size-1); j>i; j--)
         {
                 sum -= upper[i][j]*x[j];
         }
         x[i] = sum;
     }
     DeleteMatrix(size,upper);
     DeleteMatrix(size,lower);
}
//will merge with next function if time
void createsubmatrix(double**Large,double**sub,int size, int*fixednodes)
{
    int si = 0, sj = 0, fixrow = 0, fixcol = 0;
    for(int li=0;li<size;li++)
    {
        while(li == fixednodes[fixrow])
        {
            li++;
            if(li == size)
            {
                continue;
            }
            fixrow++;
        }
        for(int lj=0;lj<size;lj++)
        {
            while(lj == fixednodes[fixcol])
            {
                lj++;
                if(lj == size)
                {
                    continue;
                }
                fixcol++;
            }
            sub[si][sj] = Large[li][lj];
            sj++;
        }
        fixcol = 0;
        sj = 0;
        si++;
    }
}

void PrintV(double*v,int length)
{
    std::cout<<std::endl;
    for(int i=0;i<length;i++)
    {
        std::cout<< v[i] << std::endl;
    }
}

void createsubG(double*Large,double*sub,int size, int*fixednodes)
{
    int si = 0, fix = 0;
    for(int li = 0;li<size;li++)
    {
        while(li == fixednodes[fix])
        {
            li++;
            if(li == size)
            {
                continue;
            }
            fix++;
        }
        sub[si] = Large[li] ;
        si++;
    }
}

double mulonerow(double**matrix,double*vector,int row, int length)
{
    double sum = 0;
    for(int i=0; i<length; i++)
    {
        sum += matrix[row][i]*vector[i];
    }
    return sum;
}

void freqsweep(double*uc,double*up,double*result,double deltat,double freq,int size)
{
    for(int i=0; i<size; i++)
    {
        result[i] = ((deltat*freq)*(uc[i] + up[i]))/(4*pi);
    }
}
