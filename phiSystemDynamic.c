#include "phiSystemDynamic.h"

float **creatingEmptyStateMatrices(systemDynamicParameterPtr ptrSys)
{
    float **pd = (float **)malloc(ptrSys->rows * sizeof(float *));

    if (pd == NULL)
    {
        phiErrorHandler(ALLOCATION_ERROR);
        exit(EXIT_FAILURE);
    }

    if (ptrSys->rows != ptrSys->cols)
    {
        phiErrorHandler(INCONSISTENT_ROW_COLUMN);
        exit(EXIT_FAILURE);
    }

    if (ptrSys->dt <= 0)
    {
        phiErrorHandler(SAMPLING_RATE_ERROR);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < ptrSys->rows; i++)
        pd[i] = (float *)malloc(ptrSys->cols * sizeof(float));

    return pd;
}

float **creatingEmptyInputMatrices(systemDynamicParameterPtr ptrSys)
{
    float **pd = (float **)malloc(ptrSys->rows * sizeof(float *));

    if (pd == NULL)
    {
        phiErrorHandler(ALLOCATION_ERROR);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < ptrSys->rows; i++)
        pd[i] = (float *)malloc(ptrSys->cols * sizeof(float));

    return pd;
}

void writeTheMatrices(systemDynamicParameterPtr ptrSys)
{
    printf("Printing state matrices...\n");
    for (int i = 0; i < ptrSys->rows; i++)
    {
        for (int j = 0; j < ptrSys->cols; j++)
        {
            printf("%f ", ptrSys->stateMatrices[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("Printing input matrices...\n");
    for (int i = 0; i < ptrSys->rows; i++)
    {
        for (int j = 0; j < ptrSys->inputNumber; j++)
        {
            printf("%f ", ptrSys->inputMatrices[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

float **eyeMatricesCreation(systemDynamicParameterPtr ptrSys)
{
    float **pd = (float **)malloc(ptrSys->rows * sizeof(float *));

    if (pd == NULL)
    {
        phiErrorHandler(ALLOCATION_ERROR);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < ptrSys->rows; i++)
        pd[i] = (float *)malloc(ptrSys->cols * sizeof(float));

    for (int i = 0; i < ptrSys->rows; i++)
    {
        for (int j = 0; j < ptrSys->cols; j++)
        {
            if (i == j)
            {
                pd[i][j] = 1;
            }
            else
            {
                pd[i][j] = 0;
            }
        }
    }

    return pd;
}

float **creatingEmptyMatrices(int rows, int cols)
{
    float **pd = (float **)malloc(rows * sizeof(float *));

    if (pd == NULL)
    {
        phiErrorHandler(ALLOCATION_ERROR);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < rows; i++)
        pd[i] = (float *)malloc(cols * sizeof(float));

    return pd;
}

float **phiVectorMatrixMultiplication(float **firstTerm, float **SecondTerm, int row1, int col1, int row2, int col2)
{
    float **pd = creatingEmptyMatrices(row1, col2);
    float sum = 0;

    if (pd == NULL)
    {
        phiErrorHandler(ALLOCATION_ERROR);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < row1; i++)
    {
        for (int j = 0; j < col2; j++)
        {
            pd[i][j] = 0;
        }
    }

    for (int i = 0; i < row1; i++) //row of first matrix
    {
        for (int j = 0; j < col2; j++) //column of second matrix
        {
            sum = 0;
            for (int k = 0; k < col1; k++)
            {
                sum = sum + firstTerm[i][k] * SecondTerm[k][j];
            }
            pd[i][j] = sum;
        }
    }
    return pd;
}

float **phiSkalarMatrixMultiplication(float skalarTerm, float **SecondTerm, int row, int col)
{
    float **pd = creatingEmptyMatrices(row, col);

    if (pd == NULL)
    {
        phiErrorHandler(ALLOCATION_ERROR);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            pd[i][j] = skalarTerm * SecondTerm[i][j];
        }
    }

    return pd;
}

float **phiMatrixSummation(float **firstTerm, float **SecondTerm, int row, int col)
{
    float **pd = creatingEmptyMatrices(row, col);

    if (pd == NULL)
    {
        phiErrorHandler(ALLOCATION_ERROR);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            pd[i][j] = firstTerm[i][j] + SecondTerm[i][j];
        }
    }

    return pd;
}

void phiMatrixAssignment(float **assignedTerm, float **SecondTerm, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            assignedTerm[i][j] = SecondTerm[i][j];
        }
    }
}

void niteStaticSolver(systemDynamicParameterPtr ptrSys, float finalTime, const char *fileName)
{
    ////////////////////////////////////////////////
    // internal terms
    float **DtA;
    float **eyeDtA;
    float **DtB;
    float **stateMultiplicationMatrices;
    float **inputMultiplicationMatrices;
    float **eye = eyeMatricesCreation(ptrSys);
    int numberOfLength = (int)(finalTime / ptrSys->dt);

    // internal terms
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // file creation to write txt file
    FILE *fp;

    fp = fopen(fileName, "w");

    if (fp == NULL)
    {
        phiErrorHandler(FILE_OPEN_ERROR);
        exit(EXIT_FAILURE);
    }

    // file creation to write txt file
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // creating dynamic matrices
    ptrSys->stateNow = creatingEmptyMatrices(ptrSys->rows, 1);
    ptrSys->statePre = creatingEmptyMatrices(ptrSys->rows, 1);
    ptrSys->inputNow = creatingEmptyMatrices(1, 1);

    // creating dynamic matrices
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // Solving discrete form of state space

    // the whole operations
    ptrSys->inputNow[0][0] = ptrSys->inputValue;

    DtA = phiSkalarMatrixMultiplication(ptrSys->dt, ptrSys->stateMatrices, ptrSys->rows, ptrSys->cols);
    DtB = phiSkalarMatrixMultiplication(ptrSys->dt, ptrSys->inputMatrices, ptrSys->rows, ptrSys->inputNumber);
    eyeDtA = phiMatrixSummation(eye, DtA, ptrSys->rows, ptrSys->cols);

    for (int iter = 0; iter < numberOfLength; iter++)
    {
        stateMultiplicationMatrices = phiVectorMatrixMultiplication(eyeDtA, ptrSys->statePre,
                                                                    ptrSys->rows, ptrSys->cols, ptrSys->rows, 1);
        inputMultiplicationMatrices = phiVectorMatrixMultiplication(DtB, ptrSys->inputNow,
                                                                    ptrSys->rows, 1, 1, 1);

        float currentTime = iter * ptrSys->dt;
        for (int i = 0; i < ptrSys->rows; i++)
        {
            ptrSys->stateNow[i][0] = stateMultiplicationMatrices[i][0] + inputMultiplicationMatrices[i][0];
        }

        ///////////////////////////////////////////////////////////
        // printing results

        printf("State values : ");

        for (int i = 0; i < ptrSys->rows; i++)
        {
            printf("x[%d] : %f ", i, ptrSys->statePre[i][0]);
            fprintf(fp, "%f ", ptrSys->statePre[i][0]);
        }
        printf("elapsed Time : %f seconds\n", currentTime);
        fprintf(fp, "%f \n", currentTime);

        // printing results
        ///////////////////////////////////////////////////////////

        phiMatrixAssignment(ptrSys->statePre, ptrSys->stateNow, ptrSys->rows, 1);
    }

    // Solving discrete form of state space
    ////////////////////////////////////////////////

    ///////////////////////////////////////////////////////
    // free the whole memory
    phiFree(DtA, ptrSys->rows, ptrSys->cols);
    phiFree(eyeDtA, ptrSys->rows, ptrSys->cols);
    phiFree(DtB, ptrSys->rows, ptrSys->inputNumber);
    phiFree(stateMultiplicationMatrices, ptrSys->rows, 1);
    phiFree(inputMultiplicationMatrices, ptrSys->rows, 1);

    phiFree(ptrSys->stateNow, ptrSys->rows, 1);
    phiFree(ptrSys->statePre, ptrSys->rows, 1);
    phiFree(ptrSys->inputNow, ptrSys->inputNumber, 1);

    fclose(fp);

    // free the whole memory
    ///////////////////////////////////////////////////////
}

void phiFree(float **pd, int row, int col)
{
    for (int i = 0; i < row; i++)
        free(pd[i]);

    free(pd);
}

void phiExit(systemDynamicsParameter pSys)
{

    free(pSys.stateMatrices);
    free(pSys.inputMatrices);
}

///////////////////////////////////////////////////////////////////////
// demos

void ex1Demo()
{
    systemDynamicsParameter pSys;

    pSys.rows = 1;        // number of row in state space
    pSys.cols = 1;        // number of col in state space
    pSys.inputNumber = 1; // number of input in state space
    pSys.dt = 0.001;      // sampling period

    pSys.stateMatrices = creatingEmptyStateMatrices(&pSys);
    // state matrices creation
    pSys.inputMatrices = creatingEmptyInputMatrices(&pSys);
    // input matrices creation

    // matrices value assignment
    pSys.stateMatrices[0][0] = -1;
    pSys.inputMatrices[0][0] = 1;

    printf("Writing state and input matrices \n\n");

    // write the state space matrices
    writeTheMatrices(&pSys);

    printf("\n\n");

    /// solver solution
    // select the inputValue to be given to the system
    pSys.inputValue = 10.0;

    // static time parameter solver
    niteStaticSolver(&pSys, 10, "ver1MCK.txt");

    // should be called to free memory
    phiExit(pSys);
}

void ex2Demo()
{
    systemDynamicsParameter pSys;

    pSys.rows = 2;        // number of row in state space
    pSys.cols = 2;        // number of col in state space
    pSys.inputNumber = 1; // number of input in state space
    pSys.dt = 0.001;      // sampling period

    pSys.stateMatrices = creatingEmptyStateMatrices(&pSys);
    // state matrices creation
    pSys.inputMatrices = creatingEmptyInputMatrices(&pSys);
    // input matrices creation

    // matrices value assignment
    pSys.stateMatrices[0][0] = 0;
    pSys.stateMatrices[0][1] = 1;
    pSys.stateMatrices[1][0] = -0.1;
    pSys.stateMatrices[1][1] = -1;

    pSys.inputMatrices[0][0] = 0;
    pSys.inputMatrices[1][0] = 1;

    printf("Writing state and input matrices \n\n");

    // write the state space matrices
    writeTheMatrices(&pSys);

    printf("\n\n");

    /// solver solution
    // select the inputValue to be given to the system
    pSys.inputValue = 10.0;

    // static time parameter solver
    niteStaticSolver(&pSys, 10, "ver1MCK.txt");

    // should be called to free memory
    phiExit(pSys);
}

void ex3Demo()
{
    systemDynamicsParameter pSys;

    pSys.rows = 3;        // number of row in state space
    pSys.cols = 3;        // number of col in state space
    pSys.inputNumber = 1; // number of input in state space
    pSys.dt = 0.001;      // sampling period

    pSys.stateMatrices = creatingEmptyStateMatrices(&pSys);
    // state matrices creation
    pSys.inputMatrices = creatingEmptyInputMatrices(&pSys);
    // input matrices creation

    // matrices value assignment
    pSys.stateMatrices[0][0] = 0;
    pSys.stateMatrices[0][1] = 1;
    pSys.stateMatrices[0][2] = 0;

    pSys.stateMatrices[1][0] = 0;
    pSys.stateMatrices[1][1] = 0;
    pSys.stateMatrices[1][2] = 1;

    pSys.stateMatrices[2][0] = -2;
    pSys.stateMatrices[2][1] = -3;
    pSys.stateMatrices[2][2] = -4;

    pSys.inputMatrices[0][0] = 0;
    pSys.inputMatrices[1][0] = 0;
    pSys.inputMatrices[2][0] = 1;

    printf("Writing state and input matrices \n\n");

    // write the state space matrices
    writeTheMatrices(&pSys);

    printf("\n\n");

    /// solver solution
    // select the inputValue to be given to the system
    pSys.inputValue = 2.0;

    // static time parameter solver
    niteStaticSolver(&pSys, 10, "ver1MCK.txt");

    // should be called to free memory
    phiExit(pSys);
}
// demos
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// error Handler

void phiErrorHandler(int errorType)
{
    switch (errorType)
    {
    case FILE_OPEN_ERROR:
        printf("System Dynamic Parameter files cannot be created!\n");
        break;
    case INCONSISTENT_ROW_COLUMN:
        printf("The rows and columns are not consistent!\n");
        break;
    case ALLOCATION_ERROR:
        printf("Memory allocation cannot be done!\n");
        break;
    case SAMPLING_RATE_ERROR:
        printf("Sampling period cannot be assigned to either negative or zero value!\n");
        break;

    default:
        break;
    }
}

// error Handler
//////////////////////////////////////////////////////////////////////
