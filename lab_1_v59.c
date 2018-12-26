#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int gnuplot = 1;
//=============================================

void gauss_forward(double **A, int n, int m) {
  for (int j = 0; j < n - 1; j++)
    for (int i = j + 1; i < n; i++) {
      double coeff = A[i][j] / A[j][j];

      for (int k = j; k < m; k++)
        A[i][k] -= A[j][k] * coeff;
    }
}

//=============================================

void gauss_reverse(double **A, int n, int m) {
  for (int j = n - 1; j > 0; j--)
    for (int i = j - 1; i > -1; i--) {
      double coeff = A[i][j] / A[j][j];
      A[i][m - 1] -= A[j][m - 1] * coeff;

      for (int k = j; k > j - i - 1; k--)
        A[i][k] -= A[j][k] * coeff;
    }
}

//===============================================
void assembling(double **finite_el, int order, double *b_el, double u_beg,
                double u_end, double **res, double *b, int res_size) {
  int step = order - 1;
  int end = res_size - order + 1;

  for (int i = 0; i < end; i += step) {
    for (int j = 0; j < order; j++) {
      for (int k = 0; k < order; k++) {
        res[i + j][i + k] += finite_el[j][k];
      }
      b[i + j] += b_el[j];
    }
  }
}

//=================================================

void output(int length, int num_of_nodes, double *U, const char *file_name) {
  double frac = 2. / 71.;
  double e1 = exp(3. * length * sqrt(frac));
  double e2 = exp(2. * 3. * length * sqrt(frac));
  double C2 = (10. * e1 - 4.) * e1 / (3. * (1. - e2));
  double C1 = (-10. + 4. * e1) / (3. * (1. - e2));

  FILE *file = fopen(file_name, "w");
  FILE *gnu = popen("gnuplot -persist", "w");

  double x = 0.;
  double step = (double)length / ((double)num_of_nodes - 1.);

  double er = 0;
  for (int i = 0; i < num_of_nodes; i++) {
    double y = C1 * exp(3. * x * sqrt(frac)) + C2 * exp(-3. * x * sqrt(frac)) +
               10. / 3.;
    fprintf(file, "%lf    %lf    %lf    %e\n", x, y, U[i], fabs(U[i] - y));
    if (fabs(U[i] - y) > er) {
      er = fabs(U[i] - y);
    }
    x += step;
  }
  fprintf(file, "MAX ERROR %e\n", er);

  fclose(file);

  if (gnuplot) {
    if (file_name == "res_lin.txt") {
      FILE *gnu = popen("gnuplot -persist", "w");
      fprintf(gnu, "set xrange [0:90]\n");
      fprintf(gnu, "set grid xtics ytics\n");
      fprintf(gnu, "plot \"res_lin.txt\" using 1:2 w li lw 2 lt rgb 'blue', "
                   "\"res_lin.txt\" using 1:3 w li lw 2 lt rgb 'red'\n");
      fflush(gnu);
    }
    if (file_name == "res_quad.txt") {
      FILE *gnu = popen("gnuplot -persist", "w");
      fprintf(gnu, "set xrange [0:90]\n");
      fprintf(gnu, "set grid xtics ytics\n");
      fprintf(gnu, "plot \"res_quad.txt\" using 1:2 w li lw 2 lt rgb 'blue', "
                   "\"res_quad.txt\" using 1:3 w li lw 2 lt rgb 'red'\n");
      fflush(gnu);
    }
  }
}

//=================================================

void calculate(int length, int num_of_el, double **finite_el, int order,
               double *b_el, const char *file_name) {
  int num_of_nodes = num_of_el * (order - 1) + 1;

  double *U = (double *)calloc(num_of_nodes, sizeof(double));
  U[0] = 0.0;
  U[num_of_nodes - 1] = 2.0;

  double **res = (double **)calloc(num_of_nodes, sizeof(double *));
  double *b = (double *)calloc(num_of_nodes, sizeof(double));

  for (int i = 0; i < num_of_nodes; i++)
    res[i] = (double *)calloc(num_of_nodes, sizeof(double));

  assembling(finite_el, order, b_el, U[0], U[-1], res, b, num_of_nodes);

  int aug_n = num_of_nodes - 2;
  int aug_m = num_of_nodes - 1;

  double **augmented_matrix = (double **)calloc(aug_n, sizeof(double *));

  for (int i = 0; i < aug_n; i++)
    augmented_matrix[i] = (double *)calloc(aug_m, sizeof(double));

  for (int i = 0; i < aug_n; i++) {
    for (int j = 0; j < aug_m - 1; j++)
      augmented_matrix[i][j] = res[i + 1][j + 1];

    augmented_matrix[i][aug_m - 1] = b[i + 1];
  }

  for (int i = 0; i < order - 1; i++) {
    augmented_matrix[i][aug_m - 1] -= U[0] * finite_el[i][0];
    augmented_matrix[aug_n - 1 - i][aug_m - 1] -=
        U[num_of_nodes - 1] * finite_el[order - 2 - i][order - 1];
  }

  gauss_forward(augmented_matrix, aug_n, aug_m);
  gauss_reverse(augmented_matrix, aug_n, aug_m);

  for (int i = 0; i < aug_n; i++)
    U[i + 1] = augmented_matrix[i][aug_m - 1] / augmented_matrix[i][i];

  output(length, num_of_nodes, U, file_name);
}

//=================================================

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    printf("Введите количество КЭ\n");
    exit(0);
  }

  int length = 90;
  int num_of_el = atoi(argv[1]);
  double L = (double)length / (double)num_of_el;

  //------------

  double a = 71. / L + 6. * L;
  double b = -71. / L + 3. * L;
  int order = 2;

  double **finite_el_lin = (double **)calloc(2, sizeof(double *));

  for (int i = 0; i < 2; i++)
    finite_el_lin[i] = (double *)calloc(2, sizeof(double));

  finite_el_lin[0][0] = finite_el_lin[1][1] = a;
  finite_el_lin[1][0] = finite_el_lin[0][1] = b;

  double b_el_lin[2] = {30. * L, 30. * L};

  char *file_name_lin = "res_lin.txt";

  calculate(length, num_of_el, finite_el_lin, order, b_el_lin, file_name_lin);

  for (int i = 0; i < 2; i++)
    free(finite_el_lin[i]);

  free(finite_el_lin);

  //-------------

  a = 71. * 7. / (3. * L) - (-18.) * 2. * L / 15.;
  b = -71. * 8. / (3. * L) - (-18.0) * L / 15.;
  double c = 71. / (3. * L) - (-1.) * (-18.) * L / 30.;
  double d = 71. * 16. / (3. * L) - (-18.) * 8. * L / 15.;
  order = 3;

  double **finite_el = (double **)calloc(3, sizeof(double *));

  for (int i = 0; i < 3; i++)
    finite_el[i] = (double *)calloc(3, sizeof(double));

  finite_el[0][0] = finite_el[2][2] = a;
  finite_el[0][1] = finite_el[1][0] = finite_el[1][2] = finite_el[2][1] = b;
  finite_el[0][2] = finite_el[2][0] = c;
  finite_el[1][1] = d;

  double b_el[3] = {60. * L / 6., 60. * 2. * L / 3., 60. * L / 6.};

  char *file_name_quad = "res_quad.txt";

  calculate(length, num_of_el, finite_el, order, b_el, file_name_quad);

  for (int i = 0; i < 3; i++)
    free(finite_el[i]);

  free(finite_el);

  return 0;
}
