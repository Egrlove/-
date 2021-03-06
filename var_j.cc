#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#define R1 100.0
// #define R2 1e6
// #define R3 20.0
#define R2 20.0
#define R3 1e6

#define C3 1000e-12
#define C4 1000e-12
#define C5 2e-12
#define L 1e-3
#define E_n 10.0
#define mFt 0.026
#define It 1e-12

using namespace std;

// ���������� ���������
const int n = 4;
const int gnuplot = 1;

double norm(double x[n], double **a);
double norm_f(double **a);

void gauss_forward(double **A, int n, int m)
{
    for (int j = 0; j < n - 1; j++)
        for (int i = j + 1; i < n; i++)
        {
            double coeff = A[i][j] / A[j][j];

            for (int k = j; k < m; k++)
                A[i][k] -= A[j][k] * coeff;
        }
}

//=============================================

void gauss_reverse(double **A, int n, int m)
{
    for (int j = n - 1; j > 0; j--)
        for (int i = j - 1; i > -1; i--)
        {
            double coeff = A[i][j] / A[j][j];
            A[i][m - 1] -= A[j][m - 1] * coeff;

            for (int k = j; k > j - i - 1; k--)
                A[i][k] -= A[j][k] * coeff;
        }
}
//=============================================
void print_matr(double **a)
{
    for (int i = 0; i < n; i++)
    {

        for (int j = 0; j < n + 1; j++)
        {
            printf("%.5lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("==============================\n");
}

//=============================================

double sup_norm(double const *dV)
{
    double max = fabs(dV[0]);

    for (int i = 0; i < n; i++)
        if (fabs(dV[i]) > max)
            max = fabs(dV[i]);

    return max;
}

void Phi_plus_dPhi(double *Phi, double *dPhi)
{

    for (int i = 0; i < n; i++)
    {
        Phi[i] += dPhi[i];
        // xk[i] = 15;
    }
}

//=============================================

int main()
{
    int iter = 0;
    FILE *file = fopen("result.txt", "w");

    // ��������
    double eps = 1e-4;

    double eps1 = 1e-9;
    double eps2 = 1e-5;
    bool can_step = false;
    double dt = 0.00025;
    double now_T = 0.000001;
    double T = 0.001;

    // �������
    // double a[n][n + 1];
    double **a = (double **)calloc(n, sizeof(double *));

    for (int i = 0; i < n; i++)
        a[i] = (double *)calloc(n + 1, sizeof(double));

    // ������ ������ ������
    double b[n];

    // ������ �������
    double x[n];
    double xk[n], p[n];
    double dPhi[n];
    double Phi[n];

    for (int i = 0; i < n; i++)
    {
        Phi[i] = 0.0;
        // xk[i] = 15;
    }

    for (int i = 0; i < n; i++)
    {
        dPhi[i] = 0.0;
        // xk[i] = 15;
    }

    // �������� ��������� �����������
    double I_n1 = 0.0;
    double E = 0.0;

    double I2 = 0.0;
    double U_n1_C3 = 0.0;
    double U_n1_C4 = 0.0;
    double U_n1_C5 = 0.0;

    double fi4 = 0.0;

    xk[0] = 0.0;
    xk[1] = 0.0;
    xk[2] = 0.0;
    xk[3] = 0.0;
    printf("=======================================================================\n");
    printf("=======================================================================\n");
    printf("=======================================================================\n");

    while (now_T < T)
    {
        double P = 1e-4;

        I_n1 += (dt / L) * (-dPhi[0]);
        E = E_n * sin((2 * M_PI * now_T) / P);
        I2 = It * (exp(U_n1_C5 / mFt) - 1.0);
        double I_E = E / R1;

        // {---------------------------------

        // }---------------------------------

        do
        {
            Phi_plus_dPhi(Phi, dPhi);
            iter++;
            // {---------------------------------

            a[0][0] = (dt / L);
            a[0][1] = 0.0;
            a[0][2] = 0.0;
            a[0][3] = 0.0;

            a[1][0] = 0.0;
            a[1][1] = (C4 / dt) + (1 / R2);
            a[1][2] = -(1 / R2);
            a[1][3] = 0.0;

            a[2][0] = 0.0;
            a[2][1] = -(1 / R2);
            a[2][2] = (1 / R2) + (C5 / dt) + (1 / R3) + (It * exp(U_n1_C5 / mFt)) / mFt; // надо подумать, надо ли брать от тока I2 производную
            a[2][3] = -(C5 / dt) - (1 / R3) - (It * exp(U_n1_C5 / mFt)) / mFt;           // если нет, то это вообще можно вынести за пределы цикла do_while, что я и сделала

            a[3][0] = 0.0;
            a[3][1] = 0.0;
            a[3][2] = -(C5 / dt) - (1 / R3) - (It * exp(U_n1_C5 / mFt)) / mFt;            // надо подумать, надо ли брать от тока I2 производную
            a[3][3] = (C5 / dt) + (1 / R3) + (C3 / dt) + (It * exp(U_n1_C5 / mFt)) / mFt; // надо подумать, надо ли брать от тока I2 производную

            double I_n1_iter = I_n1 + (dt / L) * (-Phi[0]);

            double I_C3 = C3 / dt * (Phi[3] - U_n1_C3);
            double I_C4 = C4 / dt * (Phi[1] - U_n1_C4);
            double I_C5 = C5 / dt * (Phi[2] - Phi[3] - U_n1_C5);

            double I_R2 = (Phi[1] - Phi[2]) / R2;
            double I_R3 = Phi[2] - Phi[3] / R3;

            a[0][4] = -(-(I_n1_iter + (dt / L) * (-Phi[0])) + I_E);
            a[1][4] = -(-I_E + I_C4 + I_R2);
            a[2][4] = -(-I_R2 + I_C5 + I_R3 + I2);
            a[3][4] = -(-I_C5 - I_R3 - I2 + I_C3);
            // }---------------------------------
            // ��������� ������ ������ ������

            print_matr(a);
            gauss_forward(a, n, n + 1);
            gauss_reverse(a, n, n + 1);
            // print_matr(a);

            for (int i = 0; i < n; i++)
                dPhi[i] = a[i][n];

            double k = sup_norm(dPhi);

            if (k < eps2)
            {
                can_step = true;

                if (k < eps1)
                    dt *= 2;
            }
            else
            {
                can_step = false;
                dt = dt / 2.0;
            }

            if (iter > 35)
            {
                cout << "ne soshlos' =(" << endl;
                return (0);
            }

        } while (not can_step);
        // } while (sup_norm(dPhi) >= eps);
        iter = 0;
        Phi_plus_dPhi(Phi, dPhi);

        // fi4 += a[3][n];
        // U_n1_C3 = a[3][n];
        // U_n1_C4 = a[1][n];
        // U_n1_C5 = a[2][n] - a[3][n];

        U_n1_C3 = Phi[3];
        U_n1_C4 = Phi[1];
        U_n1_C5 = Phi[2] - Phi[3];

        // fprintf(file, "%.20lf    %.10lf     \n", now_T, ([3][n]));
        fprintf(file, "%.20lf    %.10lf     \n", now_T, Phi[n]);

        fflush(file);

        // printf("%s", "������ �������: [ ");
        // for (int i = 0; i < n; i++)
        // {
        //     printf("%f ", xk[i]);
        // }
        now_T += dt;
    }
    fclose(file);

    if (gnuplot)
    {
        FILE *gnu = popen("gnuplot -persist", "w");
        fprintf(gnu, "set xrange [0:0.001]\n");
        fprintf(gnu, "set grid xtics ytics\n");
        fprintf(gnu, "plot \"result.txt\" using 1:2 w li lw 2 lt rgb 'blue'\n");
        fflush(gnu);
    }
    return 0;
}

//=============================================

double norm(double x[n], double **a)
{
    double sum1 = 0;
    double sum2 = 0;

    for (int i = 0; i < n; i++)
    {
        sum1 += a[i][n] * a[i][n];
    }

    sum1 = sqrtf(sum1);

    for (int i = 0; i < n; i++)
    {
        sum2 += x[i] * x[i];
    }

    sum2 = sqrtf(sum2);

    return fabs(sum1 - sum2);
}

//=============================================

double norm_f(double **a)
{
    double sum1 = 0;
    double sum2 = 0;

    for (int i = 0; i < n; i++)
    {
        sum1 += a[i][n] * a[i][n];
    }

    sum1 = sqrtf(sum1);

    return fabs(sum1);
}

// double norm_f(double **a)
// {

//     double max = a[0][n];
//     for (int i = 0; i < n; i++)
//     {
//         if (fabs(a[i][n]) > max)
//             max = fabs(a[i][n]);
//     }

//     return max;
// }

// double max = dV[0];
// for (int i = 0; i < n; i++)
//     if (fabs(dV[i]) > max)
//         max = fabs(dV[i]);
// return max;
