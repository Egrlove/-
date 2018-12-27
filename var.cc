#include <math.h>

#include <stdio.h>
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

double sup_norm(double const *dV)
{
  double max = dV[0];
  for (int i = 0; i < n; i++)
    if (fabs(dV[i]) > max)
      max = fabs(dV[i]);
  return max;
}

int main()
{
  int iter = 0;
  FILE *file = fopen("result.txt", "w");

  // ��������
  double eps = 0.00001;

  double eps1 = 1e-9;
  double eps2 = 1e-4;
  bool can_step = false;
  double delta_T = 1e-6;
  double now_T = 0.0;
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

    do
    {
      iter++;
      {
        double P = 1e-4;
        I_n1 = I_n1 + (delta_T / L) * (-Phi[0]);
        E = E_n * sin((2 * 3.14 * now_T) / P);
        I2 = It * (exp(U_n1_C5 / mFt) - 1.0);

        a[0][4] = -(-(I_n1 + (delta_T / L) * Phi[0]) + E / R1);
        a[1][4] = -(-E / R1 + (C4 / delta_T * (Phi[1] - U_n1_C4)) + (Phi[1] - Phi[2]) / R2);
        a[2][4] = -(-(Phi[1] - Phi[2]) / R2 + (C5 / delta_T * (Phi[2] - Phi[3] - U_n1_C5)) + (Phi[2] - Phi[3]) / R3 + I2);
        a[3][4] = -(-(C5 / delta_T * (Phi[2] - Phi[3] - U_n1_C5)) - (Phi[2] - Phi[3]) / R3 - I2 + (C3 / delta_T * (Phi[3] - U_n1_C3)));
      }

      {
        a[0][0] = (delta_T / L) + (1 / R1);
        a[0][1] = -(1 / R1);
        a[0][2] = 0.0;
        a[0][3] = 0.0;

        a[1][0] = -(1 / R1);
        a[1][1] = (1 / R1) + (C4 / delta_T) + (1 / R2);
        a[1][2] = -(1 / R2);
        a[1][3] = 0.0;

        a[2][0] = 0.0;
        a[2][1] = -(1 / R2);
        a[2][2] = (1 / R2) + (C5 / delta_T) + (1 / R3);
        a[2][3] = -(C5 / delta_T) - (1 / R3);

        a[3][0] = 0.0;
        a[3][1] = 0.0;
        a[3][2] = -(C5 / delta_T) - (1 / R3) + (It * exp(U_n1_C5 / mFt)) / mFt;
        a[3][3] = (C5 / delta_T) + (1 / R3) + (C3 / delta_T) - (It * exp(U_n1_C5 / mFt) / mFt);
      }
      // ��������� ������ ������ ������

      print_matr(a);
      gauss_forward(a, n, n + 1);
      gauss_reverse(a, n, n + 1);
      // print_matr(a);

      for (int i = 0; i < n; i++)
      {
        // dPhi[i] = (a[i][n] - xk[i]);
        dPhi[i] = (a[i][n]);
      }

      if (sup_norm(dPhi) < eps1)
      {
        delta_T = 2.0 * delta_T;
        can_step = true;
      }
      else if (sup_norm(dPhi) < eps2)
      {
        can_step = true;
      }
      else
      {
        can_step = false;
        delta_T = delta_T / 2.0;
      }

      for (int i = 0; i < n; i++)
      {
        Phi[i] += dPhi[i];
        // xk[i] = 15;
      }

      if (iter > 25)
      {
        // cout<<"ne soshlos' =("<<endl;
        return (0);
      }

    } while (sup_norm(dPhi) >= eps);
    // } while (not can_step);
    iter = 0;
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
    now_T += delta_T;
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
