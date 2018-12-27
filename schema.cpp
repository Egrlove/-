#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <cstring>

#define R 1.0e3
#define L 2.53e-4
#define C 1.0e-6
#define C3 1.0e-6
#define It 1.0e-12
#define mYt 0.026
#define Re 0.0001
using namespace std;
double Rd[2] = {1.0e6, 1.0e6};
double Cd[2] = {2e-12, 2e-12};
double r[2] = {20.0, 20.0};
double dt[2] = {1.0e-8, 1.0e-8};
double In = 0.0;
double Ii = 0.0;
double Un[4] = {0.0};
double eps[2] = {2.5e-3, 1.0e-2};
struct XY
{
  double point[3];
  double time;
};
struct XY *Points = new struct XY[10000000];
void obnul(double *f, int n)
{
  for (size_t i = 0; i < n; i++)
  {
    f[i] = 0.0;
  }
}
double dIpodF(double f, int n)
{
  double x = It * exp(f / mYt) / mYt + 1 / Rd[n] + Cd[n] / dt[1];
  return x;
}
double Idiod(double f, int n)
{
  double x = It * (exp(f / mYt) - 1) + f / Rd[n] + Cd[n] * (f - Un[n + 1]) / dt[1];
  return x;
}
void InitYacc(double **Yacc, double *df)
{
  double Y01 = dt[1] / L + C / dt[1];
  double Y12 = 1 / r[1];
  double Yd1 = It * exp((df[4] - df[1]) / mYt) / mYt + 1 / Rd[0] + Cd[0] / dt[1];
  double Yd2 = It * exp((df[2] - df[3]) / mYt) / mYt + 1 / Rd[1] + Cd[1] / dt[1];
  double Yf3 = C3 / dt[1] + 1 / R;
  double Yf4 = 1 / r[0];
  Yacc[0][0] = Y01 - 1 / Re;
  Yacc[0][1] = -Y01;

  Yacc[1][0] = -Y01;
  Yacc[1][1] = Y01 + Y12 + Yd1;
  Yacc[1][2] = -Y12;
  Yacc[1][4] = -Yd1;

  Yacc[2][1] = -Y12;
  Yacc[2][2] = Y12 + Yd2;
  Yacc[2][3] = -Yd2;

  Yacc[3][2] = -Yd2;
  Yacc[3][3] = Yd2 + Yf3;
  Yacc[4][1] = -Yd1;
  Yacc[4][4] = Yd1 + Yf4;
  Yacc[0][2] = Yacc[0][3] = Yacc[0][4] = Yacc[1][3] = 0.0;
  Yacc[2][0] = Yacc[3][0] = Yacc[4][0] = Yacc[3][1] = 0.0;
  Yacc[2][4] = Yacc[3][4] = Yacc[4][2] = Yacc[4][3] = 0.0;
}
void InitTok(int t, double *df, double *I)
{
  double k = pow(10, -7);
  double Fd1 = It * (exp((df[4] - df[1]) / mYt) - 1) +
               (df[4] - df[1]) / Rd[0] +
               Cd[0] * (df[4] - df[1] - Un[1]) / dt[1];
  double Fd2 = It * (exp((df[2] - df[3]) / mYt) - 1) +
               (df[2] - df[3]) / Rd[1] +
               Cd[1] * (df[2] - df[3] - Un[2]) / dt[1];

  double tik = t * dt[1];
  I[0] = -(10 * sin(tik * 2 * M_PIl / 1e-4) / Re - df[0] / Re + ((df[0] - df[1]) * dt[1] / L) + In +
           C * (df[0] - df[1] - Un[0]) / dt[1]);
  //     cout<<10*sin(t*2*3.1415/10)/Re<<endl<<endl;
  I[1] = -(-In - (df[0] - df[1]) * dt[1] / L - C * (df[0] - df[1] - Un[0]) / dt[1] -
           Fd1 - (df[2] - df[1]) / r[1]);
  I[2] = -(Fd2 + (df[2] - df[1]) / r[1]);
  I[3] = -(-Fd2 + C3 * (df[3] - Un[3]) / dt[1] + df[3] / R);
  I[4] = -(Fd1 + df[4] / r[0]);
}
void MtrxCpy(double **Mtrx1, double **Mtrx2, int n)
{
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < n; j++)
    {
      Mtrx2[i][j] = Mtrx1[i][j];
    }
  }
}
double *VctrCpy(double *Vctr1, double *Vctr2, int n)
{
  for (size_t i = 0; i < n; i++)
  {
    Vctr2[i] = Vctr1[i];
    //  cout<< Vctr1[i]<<"  ";
  }

  return Vctr2;
}
void PrintMtrx(double **Mtrx, int n)
{
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < n; j++)
    {
      cout << Mtrx[i][j] << "     ";
    }
    cout << endl;
  }
}
void PrintVctr(double *Vctr, int n)
{
  for (size_t i = 0; i < n; i++)
  {
    cout << Vctr[i] << "  ";
  }
  cout << endl;
}
void PrintI(double *I)
{
  for (size_t i = 0; i < 5; i++)
  {
    cout << I[i] << "  ";
  }
  cout << endl;
}
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
double NormaSupr(double *f, int n)
{
  double q = fabs(f[0]);
  for (size_t i = 1; i < n; i++)
  {
    if (fabs(f[i]) > q)
      q = fabs(f[i]);
    i++;
  }
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
    cout << "uvelichenie shaga" << endl;
  }
  else if (p < eps[1])
    ;
  else if (dt[1] <= 1.0e-15)
    return 0;
  else
  {
    dt[1] = 0.5 * dt[0];
    cout << "umen'shenie shaga" << dt[1] << endl;
  }
  return 1;
}
void plusdf(double *df, double *dfper, int n)
{
  for (size_t i = 0; i < n; i++)
  {
    df[i] += dfper[i];
  }
}
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
int main(int argc, char const *argv[])
{
  FILE *gnuplt = popen("gnuplot -persist", "w");
  FILE *gnuplt1 = popen("gnuplot -persist", "w");
  int j = 0;
  int n = 5;
  int key = 0, keypr2 = 1;
  double Un1[5] = {0.0};
  double In1 = 0.0;
  double T = 1.0e-3 / dt[1];
  int t = 0;
  double pp;
  double q = pow(10, -6);
  double q1;
  double tim = dt[1] * t;
  double *I = new double[n];
  double *dfn1 = new double[n];
  double *dfn2 = new double[n];
  double *df = new double[n];
  double *dfper = new double[n];
  double **Yacc = new double *[n];
  for (size_t i = 0; i < n; i++)
  {
    Yacc[i] = new double[5];
  }
  double key1 = 0.0;
  while (t < T)
  {
    do
    {
      plusdf(df, dfper, n);
      obnul(dfper, n);
      InitYacc(Yacc, df);
      InitTok(t, df, I);
      MethodGauss(n, Yacc, I, dfper);
      q1 = NormaSupr(dfper, n);
      //  q1 = Norma2(dfper,n);
      if (key > 7)
      {
        cout << "ne soshlos' =(" << endl;
        break;
      }
      key++;
    } while (q1 > q);
    if (key > 7)
      break;
    key = 0;

    plusdf(df, dfper, n);
    obnul(dfper, n);
    if (tim >= key1)
    {
      //    cout<<tim<<endl;
      PrintVctr(df, n);
      std::cout << In << '\n';
      Points[j].point[0] = df[0];
      Points[j].point[1] = df[1];
      Points[j].point[2] = df[3];
      Points[j].time = tim;
      j++;
      key1 += 1e-6;
    }
    if (t % 1000000 == 0)
      std::cout << "solve " << t / 1000000 << " mln iter" << '\n';
    if (t > 3)
    {
      pp = pr2(dfn2, dfn1, df, n);
      //cout<<pp<<endl;
      //  PrintVctr(df,n);
      keypr2 = dtmod(pp);
      T = 1.0e-3 / dt[1];
      if (!keypr2)
      {
        cout << "ne soshlos' mdeee =(" << endl;
        break;
      }
    }
    VctrCpy(dfn1, dfn2, n);
    VctrCpy(df, dfn1, n);

    Un[0] = Un1[0];
    Un[1] = Un1[1];
    Un[2] = Un1[2];
    Un[3] = Un1[3];
    In = In;
    Un[0] = df[0] - df[1];
    Un[1] = df[4] - df[1];
    Un[2] = df[2] - df[3];
    Un[3] = df[3];
    In = In + (df[0] - df[1]) * dt[1] / L;

    tim += dt[1];

    dt[0] = dt[1];
    t++;
  }
  cout << "proizvedeno " << t << " iterac" << endl;
  cout << " Press 'F' to output plots" << endl;
  delete I;
  delete dfper;
  delete df;
  fprintf(gnuplt, "plot '-' using 1:2 with lines\n");
  fprintf(gnuplt1, "plot '-' using 1:2 with lines\n");
  int ii = 0;
  for (int i = 0; i < j; i++)
  {
    fprintf(gnuplt, "%lf %lf\n", Points[i].time, Points[i].point[0]);
    fprintf(gnuplt1, "%lf %lf\n", Points[i].time, Points[i].point[2]);
  }
  fprintf(gnuplt, "%s\n", "e");
  fprintf(gnuplt1, "%s\n", "e");
  scanf("%d\n", &t);
  return 0;
}
