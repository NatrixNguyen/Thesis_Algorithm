#define _CRT_SECURE_NO_DEPRECATE
//#include <conio.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define fopen_s(fp, fmt, mode)          *(fp)=fopen( (fmt), (mode))

double Pa[100], Pn[100];
double Ma[100][100], Mn[100][100];
double Na[100][100], Nn[100][100];
double Y_paa[100][100], Y_pan[100][100];
double Y_pra[100][100], Y_prn[100][100];
double V_paa[100][100], V_pan[100][100];
double V_pra[100][100], V_prn[100][100];
double Phi_a, Phi_n;
double PAS_a[500], PAS_n[500];
double PAS_amax, PAS_nmax;
double PAS_amin, PAS_nmin;
double dot, eucl;

double Mca[100][100], Mcn[100][100];
double Nca[100][100], Ncn[100][100];
double Yc_paa[100][100], Yc_pan[100][100];
double Yc_pra[100][100], Yc_prn[100][100];
double Vc_paa[100][100], Vc_pan[100][100];
double Vc_pra[100][100],  Vc_prn[100][100];
double Phic_a, Phic_n;
double PASc_a[500], PASc_n[500];

double S1[100][100], S2[100][100];
double mid1[100][100];
double mid2[100][100];
double c1[100][100];

double f, w, Ts;

int i, j, q, cnt, m, k, alpha, sample, count, length;

const double pi = 22.0 / 7.0;

void MatrixMultiply(double a[][100], double b[][100], int row1, int col1, int col2, double c[][100]) {

	int i, j, q;
	for (i = 0; i < row1; i++) {
		for (j = 0; j < col2; j++) {
			c[i][j] = 0;
			for (q = 0; q < col1; q++)
			{
				c[i][j] += (a[i][q] * b[q][j]);
			}
		}
	}
}

//Dot product of two vector
double DotProduct(double a[][100], double b[][100], int k) {
	int i;
	double c = 0;
	for (i = 0; i < k; i++) {
		c += a[i][0] * b[i][0];
	}
	return c;
}

//Euclidian length of vector
double Euclidian(double a[][100], double b[][100], int k) {
	int i;
	double c = 0;
	double d = 0;
	double e = 0;
	for (i = 0; i < k; i++) {
		c += pow(a[i][0], 2);
		d += pow(b[i][0], 2);
	}

	e = sqrt(c) * sqrt(d);
	return e;
}

/*For calculating Determinant of the Matrix */
double determinant(double a[][100], double k)
{
	double s = 1, det = 0, b[100][100];
	int i, j, m, n, c;
	if (k == 1)
	{
		return (a[0][0]);
	}
	else
	{
		det = 0;
		for (c = 0; c < k; c++)
		{
			m = 0;
			n = 0;
			for (i = 0; i < k; i++)
			{
				for (j = 0; j < k; j++)
				{
					b[i][j] = 0;
					if (i != 0 && j != c)
					{
						b[m][n] = a[i][j];
						if (n < (k - 2))
							n++;
						else
						{
							n = 0;
							m++;
						}
					}
				}
			}
			det = det + s * (a[0][c] * determinant(b, k - 1));
			s = -1 * s;
		}
	}

	return (det);
}

void InvertMatrix(double num[100][100], double f, double inverse[][100])
{
	double b[100][100], fac[100][100], d;
	int p, q, m, n, i, j;
	for (q = 0; q < f; q++)
	{
		for (p = 0; p < f; p++)
		{
			m = 0;
			n = 0;
			for (i = 0; i < f; i++)
			{
				for (j = 0; j < f; j++)
				{
					if (i != q && j != p)
					{
						b[m][n] = num[i][j];
						if (n < (f - 2))
							n++;
						else
						{
							n = 0;
							m++;
						}
					}
				}
			}
			fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
		}
	}

	for (i = 0; i < f; i++)
	{
		for (j = 0; j < f; j++)
		{
			b[i][j] = fac[j][i];
		}
	}
	d = determinant(num, f);
	for (i = 0; i < f; i++)
	{
		for (j = 0; j < f; j++)
		{
			inverse[i][j] = b[i][j] / d;
		}
	}

}

//input for trix H
void inputH(double H1[][100], double H2[][100], double current[100], int m) {
	double temp;
	for (i = 0; i < m - 4; i++) {

		for (j = 0; j < 3; j++) {
			H1[i][j] = current[i + j];
		}
	}
	for (i = 0; i < m - 4; i++) {
		temp = H1[i][0];
		H1[i][0] = H1[i][2];
		H1[i][2] = temp;
	}
	for (i = 0; i < m - 4; i++) {

		for (j = 0; j < 3; j++) {
			H2[j][i] = H1[i][j];
		}
	}
}

//print output (just for test the algorithm)
void output(double a[][100], int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			printf("%lf ", a[i][j]);
		}
		printf("\n");
	}
}

typedef struct
{
	double voltage_a;
	double voltage_n;

	double current_a;
	double current_n;

	//double voltage_dc;
	//double current_dc;
	//double voltage_fa;
	//double voltage_fn;

	//double current_fa;
	//double current_fn;
} Signal;


int main(void) {

	clock_t start, end;   // Khai báo biến thời gian
	double time_use;      // Thời gian sử dụng
	start = clock();     // Lấy thời gian trước khi thực hiện thuật toán

	FILE *file1;
	//system("pwd");
	file1 = fopen("/Users/namnguyen/Downloads/senior_project_NGUYEN_PHUONG_NAM/Coding/5_test_power_of_pv_and_load_1/Algorithm/Data/AC_AC_NG1.txt", "r");

	if (file1 == NULL) {
		printf("Error\n");
		//printf("Error: %d (%s)\n", errno, strerror(errno));
		//perror("Error opening file. \n");
		return 1;
	}

	volatile int samplefloat = 401;

	volatile Signal S[401];


	int read = 0;
	int records = 0;
	do {
		read = fscanf(file1, "%lf, %lf, %lf, %lf",
			&S[records].voltage_a,
			&S[records].voltage_n,

			&S[records].current_a,
			&S[records].current_n );


		if (read <= samplefloat) records++;

		if (read == samplefloat && (!feof(file1)))
		{
			printf("File format incorrect.\n");
			return 1;
		}
		if (ferror(file1))
		{
			printf("Error reading file.\n");
			return 1;
		}
	} while (!feof(file1));

	fclose(file1);


	printf("\n%d records read.\n\n", records);


	for (i = 0; i < 50; i++)
	{
		printf("%d %f %f %f %f",
			i,
			S[i].voltage_a,
			S[i].voltage_n,

			S[i].current_a,
			S[i].current_n

		);
		printf("\n");
	}

	f = 50;
	w = 2 * pi * 50;
	Ts = 0.000526;
	k = 38;
	alpha = 6;
	sample = k + alpha;
	count = 0; cnt = 0;

	for (i = 0; i < k; i++) {
		for (j = 0; j < 2; j++) {
			if (j == 0) {
				S1[i][j] = cos(w * (i + 1) * Ts);
			}
			if (j == 1) {
				S1[i][j] = sin(w * (i + 1) * Ts);
			}
			S2[j][i] = S1[i][j];
		}
	}
	MatrixMultiply(S2, S1, 2, k, 2, mid1);
	InvertMatrix(mid1, 2, mid2);
	MatrixMultiply(mid2, S2, 2, 2, k, c1);
	while (count < 4) {
		count++;
		if (count < sample) {
			Pa[count] = S[count].voltage_a;
			Pn[count] = S[count].voltage_n;
		}
		else {
			cnt++;
			for (i = 0; i < sample - 1; i++) {
				Pa[i] = Pa[i + 1];
				Pn[i] = Pn[i + 1];
			}

			Pa[sample - 1] = S[count].voltage_a;
			Pn[sample - 1] = S[count].voltage_n;

			for (i = 0; i < k; i++) {
				Ma[i][0] = Pa[i];
				Mn[i][0] = Pn[i];

				Na[i][0] = Pa[i + alpha];
				Nn[i][0] = Pn[i + alpha];
			}
		MatrixMultiply(c1, Ma, 2, k, 1, Y_paa);
		MatrixMultiply(c1, Mn, 2, k, 1, Y_pan);
		MatrixMultiply(c1, Mca, 2, k, 1, Yc_paa);
		MatrixMultiply(c1, Mcn, 2, k, 1, Yc_pan);
		MatrixMultiply(c1, Na, 2, k, 1, Y_pra);
		MatrixMultiply(c1, Nn, 2, k, 1, Y_prn);
		MatrixMultiply(c1, Nca, 2, k, 1, Yc_pra);
		MatrixMultiply(c1, Ncn, 2, k, 1, Yc_prn);
		MatrixMultiply(S1, Y_paa, k, 2, 1, V_paa);
		MatrixMultiply(S1, Y_pan, k, 2, 1, V_pan);
		MatrixMultiply(S1, Y_pra, k, 2, 1, V_pra);
		MatrixMultiply(S1, Y_prn, k, 2, 1, V_prn);
		MatrixMultiply(S1, Yc_paa, k, 2, 1, Vc_paa);
		MatrixMultiply(S1, Yc_pan, k, 2, 1, Vc_pan);
		MatrixMultiply(S1, Yc_pra, k, 2, 1, Vc_pra);
		MatrixMultiply(S1, Yc_prn, k, 2, 1, Vc_prn);

		dot = DotProduct(V_pra, V_paa, k);
		eucl = Euclidian(V_pra, V_paa, k);
		Phi_a = dot / eucl;
		Phi_a = acos(Phi_a);

		dot = DotProduct(V_prn, V_pan, k);
		eucl = Euclidian(V_prn, V_pan, k);
		Phi_n = dot / eucl;
		Phi_n = acos(Phi_n);

		dot = DotProduct(Vc_pra, Vc_paa, k);
		eucl = Euclidian(Vc_pra, Vc_paa, k);
		Phic_a = dot / eucl;
		Phic_a = acos(Phic_a);

		dot = DotProduct(Vc_prn, Vc_pan, k);
		eucl = Euclidian(Vc_prn, Vc_pan, k);
		Phic_n = dot / eucl;
		Phic_n = acos(Phic_n);

		PAS_a[count] = fabs((Phi_a * 180 / pi) - ((alpha * 360) / k));
		PAS_n[count] = fabs((Phi_n * 180 / pi) - ((alpha * 360) / k));

		PASc_a[count] = fabs((Phic_a * 180 / pi) - ((alpha * 360) / k));
		PASc_n[count] = fabs((Phic_n * 180 / pi) - ((alpha * 360) / k));

		PAS_amax = fmax(PAS_amax, PAS_a[count]);
		PAS_amin = fmin(PAS_amin, PAS_a[count]);

		PAS_nmax = fmax(PAS_nmax, PAS_n[count]);
		PAS_nmin = fmin(PAS_nmin, PAS_n[count]);
    }
    }
    /*for (int i = 0; i < sample; i++)
    {
        if (PAS_amax > 0)
        {
            printf("Fault Detected");
        }
        else
        {
            printf("Normal condition, keep execution");
        }
        return 1;
    }
	*/

	
	end = clock();  // lấy thời gian sau khi thực hiện 
	time_use = (double)(end - start) / CLOCKS_PER_SEC;    //Tính thời gian sử dụng

	printf("Time used = %lf\n", time_use);

	//FILE *file;
	//system("pwd");
	//fopen_s(&file, "/Users/namnguyen/Downloads/senior_project_NGUYEN_PHUONG_NAM/Coding/5_test_power_of_pv_and_load_1/Algorithm/Data/AC_AC_NG1_output.csv", "w");
	FILE *fp;
    fp = fopen("filenameeeee.txt", "w");
	if (fp == NULL) {
        printf("Unable to open the output file.\n");
        return 1;
    }
	for (int i = 0; i < 401; i++) {
	fprintf(fp, "%lf, %lf, %lf, %lf \n", 
			PAS_a[i], PAS_n[i], PASc_a[i], PASc_n[i]);
	printf("\n");		
	}

	fclose(fp);

	return 0;
    }
