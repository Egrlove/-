#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#define R1 1e-6
#define R2 20.0
#define R3 1e6
#define R4 1000.0

#define C3 1000e-12
#define C4 1000e-12
#define C5 2e-12
#define L 1e-3
#define E_n 10.0
#define mFt 0.026
#define It 1e-12
double dt[2] = {1.0e-8, 1.0e-8};
double eps[2] = {2.5e-3, 1.0e-2};

using namespace std;

// ���������� ���������
const int n = 4;
const int gnuplot = 1;

double norm(double x[n], double **a);
double norm_f(double **a);

//=============================================
void print_matr(double **a)
{
    for (int i = 0; i < n; i++)
    {

        for (int j = 0; j < n; j++)
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
int MethodGauss(int n, double **a, double *y, double *x)
{
    double max;
    int k, index;
    const double eps = 1e-21; // точность
    k = 0;

    while (k < n)
    {
        // Поиск строки с максимальным a[i][k]
        double r1 = a[k][k];
        max = fabs(r1);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            double r2 = a[i][k];
            double q1 = fabs(r2);
            if (q1 > max)
            {
                max = q1;

                index = i;
            }
            // cout<<max<<endl;
        }
        // Перестановка строк
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            return 10000;
        }
        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            double r = fabs(temp);
            if (r < eps)
                continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)
                continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
            double m = a[i][k] / a[k][k];
            //  for (size_t j = 0; j < n; j++) {
            //    a[i][j] = a[i][j] - m*a[k][j];
            //
            //  }
            //  y[i]-=m*y[k];
        }
        k++;
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; k--)
    {

        for (int i = n - 1; i > k; i--)
            y[k] -= a[k][i] * x[i];
        x[k] = y[k] / a[k][k];
    }

    //  for (k = n - 1; k >= 0; k--)
    //  {  x[k]=y[k]/a[k][k];
    //    for (size_t i = n-1; i > k; i--) {
    //      x[k]-=x[i]*a[k][i]/a[k][k];
    //    }
    //  }
}

//=============================================

double Norma2(double *f, int n)
{
    double q = 0;
    for (size_t i = 0; i < n; i++)
    {
        q += f[i] * f[i];
    }
    q = sqrt(q);
    return q;
}

double pr2(double *dfn2, double *dfn1, double *dfn, int n)
{
    double pp;
    double *p1 = new double[n];
    for (size_t i = 0; i < n; i++)
    {
        p1[i] = 0.5 * (((dfn2[i] - dfn1[i]) / dt[0] - (dfn1[i] - dfn[i]) / dt[1]) / dt[1]) * dt[1] * dt[1];
    }
    pp = Norma2(p1, n);
    delete p1;
    return pp;
}

int dtmod(double p)
{
    if (p < eps[0] && dt[1] < 1.0e-8)
    {
        dt[1] = 2 * dt[0];
        // cout << "uvelichenie shaga" << endl;
    }
    else if (p < eps[1])
        ;
    else if (dt[1] <= 1.0e-15)
        return 0;
    else
    {
        dt[1] = 0.5 * dt[0];
        // cout << "umen'shenie shaga" << dt[1] << endl;
    }
    return 1;
}

double *VctrCpy(double *Vctr1, double *Vctr2, int n)
{
    for (size_t i = 0; i < n; i++)
    {
        Vctr2[i] = Vctr1[i];
    }

    return Vctr2;
}

int main()
{
    int iter = 0;
    FILE *file = fopen("result.txt", "w");

    double eps = 1e-7;
    double eps1 = 1e-12;
    double eps2 = 1e-6;
    bool can_step = false;
    double dt = 1e-8;
    double now_T = 0.000001;
    double T = 0.001;

    // �������
    // double a[n][n + 1];
    double **a = (double **)calloc(n, sizeof(double *));

    for (int i = 0; i < n; i++)
        a[i] = (double *)calloc(n, sizeof(double));

    double x[n];
    double xk[n], p[n];
    double dPhi[n];
    double Phi[n];

    double y[n];
    double *dfn1 = new double[n];
    double *dfn2 = new double[n];

    for (int i = 0; i < n; i++)
    {
        Phi[i] = 0.0;
    }

    for (int i = 0; i < n; i++)
    {
        dPhi[i] = 0.0;
    }

    // �������� ��������� �����������
    double I_n1_iter = 0.0;
    double E = 0.0;

    double I2 = 0.0;
    double U_n1_C3 = 0.0;
    double U_n1_C4 = 0.0;
    double U_n1_C5 = 0.0;

    // printf("=======================================================================\n");
    // printf("=======================================================================\n");
    // printf("=======================================================================\n");

    while (now_T < T)
    {
        double P = 1e-4;
        E = E_n * sin((2 * M_PI * now_T) / P);
        double I_E = E / R1;
        do
        {
            iter++;
            // {---------------------------------
            Phi_plus_dPhi(Phi, dPhi);
            for (int i = 0; i < n; i++)
                dPhi[i] = 0.0;

            I2 = It * (exp((Phi[2] - Phi[3]) / mFt) - 1.0);

            double I_C3 = (C3 / dt) * (Phi[3] - U_n1_C3);
            double I_C4 = (C4 / dt) * (Phi[1] - U_n1_C4);
            double I_C5 = (C5 / dt) * (Phi[2] - Phi[3] - U_n1_C5);

            double I_R2 = (Phi[1] - Phi[2]) / R2;
            double I_R3 = (Phi[2] - Phi[3]) / R3;
            {
                a[0][0] = (dt / L) + (1 / R1);
                a[0][1] = -(1 / R1);
                a[0][2] = 0.0;
                a[0][3] = 0.0;

                a[1][0] = -(1 / R1);
                a[1][1] = (1 / R1) + (C4 / dt) + (1 / R2);
                a[1][2] = -(1 / R2);
                a[1][3] = 0.0;

                a[2][0] = 0.0;
                a[2][1] = -(1 / R2);
                a[2][2] = (1 / R2) + (C5 / dt) + (1 / R3) + (It * exp((Phi[2] - Phi[3]) / mFt)) / mFt;
                a[2][3] = -(C5 / dt) - (1 / R3) - (It * exp((Phi[2] - Phi[3]) / mFt)) / mFt;

                a[3][0] = 0.0;
                a[3][1] = 0.0;
                a[3][2] = -(C5 / dt) - (1 / R3) - (It * exp((Phi[2] - Phi[3]) / mFt)) / mFt;
                a[3][3] = (C5 / dt) + (1 / R3) + (1 / R4) + (It * exp((Phi[2] - Phi[3]) / mFt)) / mFt;

                y[0] = -(-(I_n1_iter + (dt / L) * (-Phi[0])) + I_E + (Phi[0] - Phi[1]) / R1);
                y[1] = -(-I_E + I_C4 + I_R2 - (Phi[0] - Phi[1]) / R1);
                y[2] = -(-I_R2 + I_C5 + I_R3 + I2);
                y[3] = -(-I_C5 - I_R3 - I2 + Phi[3]/R4);
            }
            
            {
                // print_matr(a);
                // printf("Y\n");
                // for (int i = 0; i < n; i++)
                //     printf("%.20lf  ", y[i]);
                // printf("\n");

                // printf("DPHI\n");
                // for (int i = 0; i < n; i++)
                //     printf("%.20lf  ", dPhi[i]);
                // printf("\n");

                // printf("PHI\n");
                // for (int i = 0; i < n; i++)
                //     printf("%.20lf  ", Phi[i]);
                // printf("\n");

                MethodGauss(n, a, y, dPhi);
            }

            if (iter > 25)
            {
                cout << "ne soshlos' =(" << endl;
                return (0);
            }

        } while (sup_norm(dPhi) >= eps);
        // printf("=========================STEP====================================\n");

        iter = 0;

        Phi_plus_dPhi(Phi, dPhi);
        for (int i = 0; i < n; i++)
            dPhi[i] = 0.0;

        if (now_T > 0.000002)
        {
            // double k = pr2(dfn2, dfn1, dPhi, n);
            //cout<<pp<<endl;
            //  PrintVctr(df,n);

            // dtmod(k);
            double k = sup_norm(dPhi);

            if (k < eps2)
            {
                can_step = true;

                // if (k < eps1)
                //     dt *= 1.5;
            }
            else
            {
                can_step = false;
                dt = dt / 2.0;
            }
        }
        VctrCpy(dfn1, dfn2, n);
        VctrCpy(dPhi, dfn1, n);

        U_n1_C3 = Phi[3];
        U_n1_C4 = Phi[1];
        U_n1_C5 = Phi[2] - Phi[3];
        I_n1_iter = I_n1_iter + (dt / L) * (-Phi[0]);

        fprintf(file, "%.20lf    %.10lf     \n", now_T, Phi[3]);
        fflush(file);

        now_T += dt;
    }
    fclose(file);

    if (gnuplot)
    {
        FILE *gnu = popen("gnuplot -persist", "w");
        fprintf(gnu, "set xrange [0:0.001]\n");
        // fprintf(gnu, "set yrange [-0.1:0.1]\n");
        fprintf(gnu, "set grid xtics ytics\n");
        fprintf(gnu, "plot \"result.txt\" using 1:2 w li lw 2 lt rgb 'blue'\n");
        fflush(gnu);
    }
    return 0;
}

//============================================
double
norm(double x[n], double **a)
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
