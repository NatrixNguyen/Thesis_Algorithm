

/**
 * main.c
 */

 // C program for the above approach
#define _CRT_SECURE_NO_DEPRECATE
#include <conio.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>




 //declare for PAS
double Pa[100], Pb[100], Pc[100], Pn[100];
double Ma[100][100], Mb[100][100], Mc[100][100], Mn[100][100];
double Na[100][100], Nb[100][100], Nc[100][100], Nn[100][100];
double Y_paa[100][100], Y_pab[100][100], Y_pac[100][100], Y_pan[100][100];
double Y_pra[100][100], Y_prb[100][100], Y_prc[100][100], Y_prn[100][100];
double V_paa[100][100], V_pab[100][100], V_pac[100][100], V_pan[100][100];
double V_pra[100][100], V_prb[100][100], V_prc[100][100], V_prn[100][100];
double Phi_a, Phi_b, Phi_c, Phi_n;
double PAS_a[500], PAS_b[500], PAS_c[500], PAS_n[500];
double PAS_amax, PAS_bmax, PAS_cmax, PAS_nmax;
double PAS_amin, PAS_bmin, PAS_cmin, PAS_nmin;
double dot, eucl;


double Mca[100][100], Mcb[100][100], Mcc[100][100], Mcn[100][100];
double Nca[100][100], Ncb[100][100], Ncc[100][100], Ncn[100][100];
double Yc_paa[100][100], Yc_pab[100][100], Yc_pac[100][100], Yc_pan[100][100];
double Yc_pra[100][100], Yc_prb[100][100], Yc_prc[100][100], Yc_prn[100][100];
double Vc_paa[100][100], Vc_pab[100][100], Vc_pac[100][100], Vc_pan[100][100];
double Vc_pra[100][100], Vc_prb[100][100], Vc_prc[100][100], Vc_prn[100][100];
double Phic_a, Phic_b, Phic_c, Phic_n;
double PASc_a[500], PASc_b[500], PASc_c[500], PASc_n[500];



//declare for D parameter
double current_known_a[100], current_known_b[100], current_known_c[100], current_known_n[100];
double cur_a[100], cur_b[100], cur_c[100], cur_n[100];
double H1_a[100][100], H1_b[100][100], H1_c[100][100], H1_n[100][100];
double H2_a[100][100], H2_b[100][100], H2_c[100][100], H2_n[100][100];
double mid1_a[100][100], mid1_b[100][100], mid1_c[100][100], mid1_n[100][100];
double mid2_a[100][100], mid2_b[100][100], mid2_c[100][100], mid2_n[100][100];
double mid3_a[100][100], mid3_b[100][100], mid3_c[100][100], mid3_n[100][100];
double current1_a[100][100], current1_b[100][100], current1_c[100][100], current1_n[100][100];
double coefficient_a[100][100], coefficient_b[100][100], coefficient_c[100][100], coefficient_n[100][100];
double current_measured_a, current_measured_b, current_measured_c, current_measured_n;
double current_predicted_a, current_predicted_b, current_predicted_c, current_predicted_n;
double Da[500], Db[500], Dc[500], Dn[500];
double minDa, maxDa;

double volt_known_a[100], volt_known_b[100], volt_known_c[100], volt_known_n[100];
double H1c_a[100][100], H1c_b[100][100], H1c_c[100][100], H1c_n[100][100];
double H2c_a[100][100], H2c_b[100][100], H2c_c[100][100], H2c_n[100][100];
double mid1c_a[100][100], mid1c_b[100][100], mid1c_c[100][100], mid1c_n[100][100];
double mid2c_a[100][100], mid2c_b[100][100], mid2c_c[100][100], mid2c_n[100][100];
double mid3c_a[100][100], mid3c_b[100][100], mid3c_c[100][100], mid3c_n[100][100];
double volt1_a[100][100], volt1_b[100][100], volt1_c[100][100], volt1_n[100][100];
double coefficientv_a[100][100], coefficientv_b[100][100], coefficientv_c[100][100], coefficientv_n[100][100];
double volt_measured_a, volt_measured_b, volt_measured_c, volt_measured_n;
double volt_predicted_a, volt_predicted_b, volt_predicted_c, volt_predicted_n;
double Dca[500], Dcb[500], Dcc[500], Dcn[500];

//declare for other elementary
double S1[100][100], S2[100][100];
double mid1[100][100];
double mid2[100][100];
double c1[100][100];

double f, w, Ts;

int i, j, q, cnt, m, k, alpha, sample, count, length;;

const double pi = 22.0 / 7.0;





//Multiply two matrix
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
	double voltage_b;
	double voltage_c;
	double voltage_n;

	double current_a;
	double current_b;
	double current_c;
	double current_n;

	double voltage_dc;
	double current_dc;
	double voltage_fa;
	double voltage_fb;
	double voltage_fc;
	double voltage_fn;

	double current_fa;
	double current_fb;
	double current_fc;
	double current_fn;
} Signal;


int main(void) {

	clock_t start, end;   // Khai báo biến thời gian
	double time_use;      // Thời gian sử dụng
	start = clock();     // Lấy thời gian trước khi thực hiện thuật toán

	FILE* file1;
	file1 = fopen("D:\\PV\\AC faults Fluke\\Test\\AC\\A-G-7.csv", "r");

	if (file1 == NULL) {
		printf("Error opening file. \n");
		return 1;
	}

	volatile int samplefloat = 401;

	volatile Signal S[401];


	int read = 0;
	int records = 0;
	do {
		read = fscanf(file1, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf",
			&S[records].voltage_a,
			&S[records].voltage_b,
			&S[records].voltage_c,
			&S[records].voltage_n,

			&S[records].current_a,
			&S[records].current_b,
			&S[records].current_c,
			&S[records].current_n,

			&S[records].voltage_dc,
			&S[records].current_dc,

			&S[records].voltage_fa,
			&S[records].voltage_fb,
			&S[records].voltage_fc,
			&S[records].voltage_fn,

			&S[records].current_fa,
			&S[records].current_fb,
			&S[records].current_fc,
			&S[records].current_fn
		);





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
		printf("%d %f %f %f %f  %f %f %f %f",
			i,
			S[i].voltage_a,
			S[i].voltage_b,
			S[i].voltage_c,
			S[i].voltage_n,

			S[i].current_a,
			S[i].current_b,
			S[i].current_c,
			S[i].current_n

		);
		printf("\n");
	}
	

	
	//Input frequency of inverter for f
	f = 50;

	//Compute w from f
	w = 2 * pi * 50;

	//Choose Ts for sampling period
	Ts = 0.000526;

	//Choose the amount of samples in one sampling window
	k = 38;

	//Choose the length between two sampling window (0.1*sample<=alpha<=0.25*sample)
	alpha = 6;

	//The amount of samples (include samples in one window plus the length between two window)
	sample = k + alpha;

	count = 0; cnt = 0;


	/*Compute constant matrix S1, S2
		S1=[cos(wTs)	sin(wTs)
			cos(w2Ts)	sin(w2Ts)
					...
			cos(wkTs)	sin(wkTs)]

		S2=[cos(wTs)	cos(w2Ts)	...	cos(wkTs)
			sin(wTs)	sin(w2Ts)	...	sin(wkTs)]
	*/

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

	// Compute C1=(S2*S1)^-1*S2

	MatrixMultiply(S2, S1, 2, k, 2, mid1);
	printf("mid1 = \n");
	output(mid1, 2, 2);
	printf("\n");
	InvertMatrix(mid1, 2, mid2);
	printf("mid2 = \n");
	output(mid2, 2, 2);
	printf("\n");
	MatrixMultiply(mid2, S2, 2, 2, k, c1);
	printf("c1 = \n");
	output(c1, 2, k);
	printf("\n");

	while (count<400) {
		
		//Phase Angle Shift
		/* Write k voltage samples from each phase to Ma, Mb, Mc as past sampling window, respectively;
		 each sample of each phase is recorded after Ts; all 3 phase are taken sample in parallel
			Ma=[ia(Ts)	ia(2Ts)	...	ia(kTs)]
			Mb=[ib(Ts)	ib(2Ts)	...	ib(kTs)]
			Mc=[ic(Ts)	ic(2Ts)	...	ic(kTs)]
		   Write k voltage samples from each phase to Na, Nb, Nc as present sampling window, respectively;
		 each sample of each phase is recorded after Ts; all 3 phase are taken sample in parallel;
		 the first sample of this array is taking after "alpha" samples
			Na=[ia(Ts)	ia(2Ts)	...	ia(kTs)]
			Nb=[ib(Ts)	ib(2Ts)	...	ib(kTs)]
			Nc=[ic(Ts)	ic(2Ts)	...	ic(kTs)]
		*/
		count++;
		if (count < sample) {
			Pa[count] = S[count].voltage_a;
			Pb[count] = S[count].voltage_b;
			Pc[count] = S[count].voltage_c;
			Pn[count] = S[count].voltage_n;

			cur_a[count] = S[count].current_a;
			cur_b[count] = S[count].current_b;
			cur_c[count] = S[count].current_c;
			cur_n[count] = S[count].current_n;

		}
		else {
			cnt++;
			for (i = 0; i < sample - 1; i++) {
				Pa[i] = Pa[i + 1];
				Pb[i] = Pb[i + 1];
				Pc[i] = Pc[i + 1];
				Pn[i] = Pn[i + 1];

				cur_a[i] = cur_a[i + 1];
				cur_b[i] = cur_b[i + 1];
				cur_c[i] = cur_c[i + 1];
				cur_n[i] = cur_n[i + 1];

			}

			Pa[sample - 1] = S[count].voltage_a;
			Pb[sample - 1] = S[count].voltage_b;
			Pc[sample - 1] = S[count].voltage_c;
			Pn[sample - 1] = S[count].voltage_n;

			cur_a[sample - 1] = S[count].current_a;
			cur_b[sample - 1] = S[count].current_b;
			cur_c[sample - 1] = S[count].current_c;
			cur_n[sample - 1] = S[count].current_n;

			for (i = 0; i < k; i++) {
				Ma[i][0] = Pa[i];
				Mb[i][0] = Pb[i];
				Mc[i][0] = Pc[i];
				Mn[i][0] = Pn[i];

				Na[i][0] = Pa[i + alpha];
				Nb[i][0] = Pb[i + alpha];
				Nc[i][0] = Pc[i + alpha];
				Nn[i][0] = Pn[i + alpha];

				Mca[i][0] = cur_a[i];
				Mcb[i][0] = cur_b[i];
				Mcc[i][0] = cur_c[i];
				Mcn[i][0] = cur_n[i];

				Nca[i][0] = cur_a[i + alpha];
				Ncb[i][0] = cur_b[i + alpha];
				Ncc[i][0] = cur_c[i + alpha];
				Ncn[i][0] = cur_n[i + alpha];
			}
			

		/*Compute Y_paa, Y_pab, Y_pac from S1, S2 and Ma, Mb, Mc
			Y_paa=c1*Ma
			Y_pab=c1*Mb
			Y_pac=c1*Mc
		*/
		//Matrix caculation

		MatrixMultiply(c1, Ma, 2, k, 1, Y_paa);
		/*printf("Ya = \n");
		output(Y_paa, 2, 1);
		printf("\n");*/
		MatrixMultiply(c1, Mb, 2, k, 1, Y_pab);
		MatrixMultiply(c1, Mc, 2, k, 1, Y_pac);
		MatrixMultiply(c1, Mn, 2, k, 1, Y_pan);

		MatrixMultiply(c1, Mca, 2, k, 1, Yc_paa);
		/*printf("Ya = \n");
		output(Yc_paa, 2, 1);
		printf("\n");*/
		MatrixMultiply(c1, Mcb, 2, k, 1, Yc_pab);
		MatrixMultiply(c1, Mcc, 2, k, 1, Yc_pac);
		MatrixMultiply(c1, Mcn, 2, k, 1, Yc_pan);

		/*Compute Y_pra, Y_prb, Y_prc from S1, S2 and Na, Nb, Nc
			Y_paa=c1*Na
			Y_pab=c1*Nb
			Y_pac=c1*Nc
		*/
		//Matrix caculation

		MatrixMultiply(c1, Na, 2, k, 1, Y_pra);
		MatrixMultiply(c1, Nb, 2, k, 1, Y_prb);
		MatrixMultiply(c1, Nc, 2, k, 1, Y_prc);
		MatrixMultiply(c1, Nn, 2, k, 1, Y_prn);

		MatrixMultiply(c1, Nca, 2, k, 1, Yc_pra);
		MatrixMultiply(c1, Ncb, 2, k, 1, Yc_prb);
		MatrixMultiply(c1, Ncc, 2, k, 1, Yc_prc);
		MatrixMultiply(c1, Ncn, 2, k, 1, Yc_prn);

		/*Compute V_paa, V_pab, V_pac, which are the estimated fundamental frequency component of Ma, Mb, Mc respectively
			V_paa=S1*Y_paa
			V_pab=S1*Y_pab
			V_pac=S1*Y_pac
		  Compute V_pra, V_prb, V_prc, which are the estimated fundamental frequency component of Na, Nb, Nc respectively
			V_pra=S1*Y_pra
			V_prb=S1*Y_prb
			V_prc=S1*Y_prc
		*/
		//Matrix caculation

		MatrixMultiply(S1, Y_paa, k, 2, 1, V_paa);
		MatrixMultiply(S1, Y_pab, k, 2, 1, V_pab);
		MatrixMultiply(S1, Y_pac, k, 2, 1, V_pac);
		MatrixMultiply(S1, Y_pan, k, 2, 1, V_pan);

		MatrixMultiply(S1, Y_pra, k, 2, 1, V_pra);
		MatrixMultiply(S1, Y_prb, k, 2, 1, V_prb);
		MatrixMultiply(S1, Y_prc, k, 2, 1, V_prc);
		MatrixMultiply(S1, Y_prn, k, 2, 1, V_prn);


		MatrixMultiply(S1, Yc_paa, k, 2, 1, Vc_paa);
		MatrixMultiply(S1, Yc_pab, k, 2, 1, Vc_pab);
		MatrixMultiply(S1, Yc_pac, k, 2, 1, Vc_pac);
		MatrixMultiply(S1, Yc_pan, k, 2, 1, Vc_pan);

		MatrixMultiply(S1, Yc_pra, k, 2, 1, Vc_pra);
		MatrixMultiply(S1, Yc_prb, k, 2, 1, Vc_prb);
		MatrixMultiply(S1, Yc_prc, k, 2, 1, Vc_prc);
		MatrixMultiply(S1, Yc_prn, k, 2, 1, Vc_prn);


		/*Find the phase difference between the two vector in radian
			Phi=arcos((V_pr.V_pa)/(abs(V_pr)*abs(V_pa))
		*/
		dot = DotProduct(V_pra, V_paa, k);
		//printf("dota = %lf \n", dot);
		eucl = Euclidian(V_pra, V_paa, k);
		//printf("Euclidiana = %lf\n", eucl);
		Phi_a = dot / eucl;
		Phi_a = acos(Phi_a);
		//printf("Phi_a = %lf\n", Phi_a * 180 / pi);

		dot = DotProduct(V_prb, V_pab, k);
		//printf("dotb = %lf \n", dot);
		eucl = Euclidian(V_prb, V_pab, k);
		//printf("Euclidianb = %lf\n", eucl);
		Phi_b = dot / eucl;
		Phi_b = acos(Phi_b);
		//printf("Phi_b = %lf\n", Phi_b * 180 / pi);

		dot = DotProduct(V_prc, V_pac, k);
		//printf("dotc = %lf \n", dot);
		eucl = Euclidian(V_prc, V_pac, k);
		//printf("Euclidian = %lf\n", eucl);
		Phi_c = dot / eucl;
		Phi_c = acos(Phi_c);
		//printf("Phi_a = %lf\n", Phi_c * 180 / pi);

		dot = DotProduct(V_prn, V_pan, k);
		//printf("dotc = %lf \n", dot);
		eucl = Euclidian(V_prn, V_pan, k);
		//printf("Euclidian = %lf\n", eucl);
		Phi_n = dot / eucl;
		Phi_n = acos(Phi_n);
		//printf("Phi_a = %lf\n", Phi_n * 180 / pi);



		dot = DotProduct(Vc_pra, Vc_paa, k);
		//printf("dota = %lf \n", dot);
		eucl = Euclidian(Vc_pra, Vc_paa, k);
		//printf("Euclidiana = %lf\n", eucl);
		Phic_a = dot / eucl;
		Phic_a = acos(Phic_a);
		//printf("Phi_a = %lf\n", Phic_a * 180 / pi);

		dot = DotProduct(Vc_prb, Vc_pab, k);
		//printf("dotb = %lf \n", dot);
		eucl = Euclidian(Vc_prb, Vc_pab, k);
		//printf("Euclidianb = %lf\n", eucl);
		Phic_b = dot / eucl;
		Phic_b = acos(Phic_b);
		//printf("Phi_b = %lf\n", Phic_b * 180 / pi);

		dot = DotProduct(Vc_prc, Vc_pac, k);
		//printf("dotc = %lf \n", dot);
		eucl = Euclidian(Vc_prc, Vc_pac, k);
		//printf("Euclidian = %lf\n", eucl);
		Phic_c = dot / eucl;
		Phic_c = acos(Phic_c);
		//printf("Phi_a = %lf\n", Phic_c * 180 / pi);

		dot = DotProduct(Vc_prn, Vc_pan, k);
		//printf("dotc = %lf \n", dot);
		eucl = Euclidian(Vc_prn, Vc_pan, k);
		//printf("Euclidian = %lf\n", eucl);
		Phic_n = dot / eucl;
		Phic_n = acos(Phic_n);
		//printf("Phi_a = %lf\n", Phic_n * 180 / pi);

		//Find PAS
		PAS_a[count] = fabs((Phi_a * 180 / pi) - ((alpha * 360) / k));
		PAS_b[count] = fabs((Phi_b * 180 / pi) - ((alpha * 360) / k));
		PAS_c[count] = fabs((Phi_c * 180 / pi) - ((alpha * 360) / k));
		PAS_n[count] = fabs((Phi_n * 180 / pi) - ((alpha * 360) / k));

		/*printf("PAS_a[%d] = %lf\n", count, PAS_a[count]);
		printf("PAS_b[%d] = %lf\n", count, PAS_b[count]);
		printf("PAS_c[%d] = %lf\n", count, PAS_c[count]);
		printf("PAS_n[%d] = %lf\n", count, PAS_n[count]);
		printf("\n");*/


		PASc_a[count] = fabs((Phic_a * 180 / pi) - ((alpha * 360) / k));
		PASc_b[count] = fabs((Phic_b * 180 / pi) - ((alpha * 360) / k));
		PASc_c[count] = fabs((Phic_c * 180 / pi) - ((alpha * 360) / k));
		PASc_n[count] = fabs((Phic_n * 180 / pi) - ((alpha * 360) / k));

		/*printf("PASc_a[%d] = %lf\n", count, PASc_a[count]);
		printf("PASc_b[%d] = %lf\n", count, PASc_b[count]);
		printf("PASc_c[%d] = %lf\n", count, PASc_c[count]);
		printf("PASc_n[%d] = %lf\n", count, PASc_n[count]);
		printf("\n");*/


		PAS_amax = fmax(PAS_amax, PAS_a[count]);
		PAS_amin = fmin(PAS_amin, PAS_a[count]);

		PAS_bmax = fmax(PAS_bmax, PAS_b[count]);
		PAS_bmin = fmin(PAS_bmin, PAS_b[count]);

		PAS_cmax = fmax(PAS_cmax, PAS_c[count]);
		PAS_cmin = fmin(PAS_cmin, PAS_c[count]);

		PAS_nmax = fmax(PAS_nmax, PAS_n[count]);
		PAS_nmin = fmin(PAS_nmin, PAS_n[count]);



		//D parameter
		m = 33;
		//copy m known samples to to array current_known
		for (i = 0; i < m; i++) {
			current_known_a[i] = cur_a[i];
			current_known_b[i] = cur_b[i];
			current_known_c[i] = cur_c[i];
			current_known_n[i] = cur_n[i];

			volt_known_a[i] = Pa[i];
			volt_known_b[i] = Pb[i];
			volt_known_c[i] = Pc[i];
			volt_known_n[i] = Pn[i];
		}

		/*printf("current_a = \n");
		output(current_a, k, 1);
		printf("\n");

		printf("current_known_a = \n");
		output(current_known_a, k, 1);
		printf("\n");*/

		//input for matrix H
		inputH(H1_a, H2_a, current_known_a, m);

		/*printf("H1= \n");
		output(H1_a, m, m);
		printf("\n");

		printf("H2= \n");
		output(H2_a, m, m);
		printf("\n");*/

		inputH(H1_b, H2_b, current_known_b, m);
		inputH(H1_c, H2_c, current_known_c, m);
		inputH(H1_n, H2_n, current_known_n, m);

		inputH(H1c_a, H2c_a, volt_known_a, m);
		inputH(H1c_b, H2c_b, volt_known_b, m);
		inputH(H1c_c, H2c_c, volt_known_c, m);
		inputH(H1c_n, H2c_n, volt_known_n, m);

		/* Calculate a=(H1*H2)^-1*H2*I */
		MatrixMultiply(H2_a, H1_a, 3, m - 4, 3, mid1_a);

		/*printf("mid1_a= \n");
		output(mid1_a, m, m);
		printf("\n");*/

		InvertMatrix(mid1_a, 3, mid2_a);

		/*printf("mid2_a= \n");
		output(mid2_a, m, m);
		printf("\n");*/

		MatrixMultiply(mid2_a, H2_a, 3, 3, m - 4, mid3_a);

		/*printf("mid3_a= \n");
		output(mid3_a, m, m);
		printf("\n");*/

		MatrixMultiply(H2_b, H1_b, 3, m - 4, 3, mid1_b);
		InvertMatrix(mid1_b, 3, mid2_b);
		MatrixMultiply(mid2_b, H2_b, 3, 3, m - 4, mid3_b);

		MatrixMultiply(H2_c, H1_c, 3, m - 4, 3, mid1_c);
		InvertMatrix(mid1_c, 3, mid2_c);
		MatrixMultiply(mid2_c, H2_c, 3, 3, m - 4, mid3_c);

		MatrixMultiply(H2_n, H1_n, 3, m - 4, 3, mid1_n);
		InvertMatrix(mid1_n, 3, mid2_n);
		MatrixMultiply(mid2_n, H2_n, 3, 3, m - 4, mid3_n);


		MatrixMultiply(H2c_a, H1c_a, 3, m - 4, 3, mid1c_a);
		InvertMatrix(mid1c_a, 3, mid2c_a);
		MatrixMultiply(mid2c_a, H2c_a, 3, 3, m - 4, mid3c_a);

		MatrixMultiply(H2c_b, H1c_b, 3, m - 4, 3, mid1c_b);
		InvertMatrix(mid1c_b, 3, mid2c_b);
		MatrixMultiply(mid2c_b, H2c_b, 3, 3, m - 4, mid3c_b);

		MatrixMultiply(H2c_c, H1c_c, 3, m - 4, 3, mid1c_c);
		InvertMatrix(mid1c_c, 3, mid2c_c);
		MatrixMultiply(mid2c_c, H2c_c, 3, 3, m - 4, mid3c_c);

		MatrixMultiply(H2c_n, H1c_n, 3, m - 4, 3, mid1c_n);
		InvertMatrix(mid1c_n, 3, mid2c_n);
		MatrixMultiply(mid2c_n, H2c_n, 3, 3, m - 4, mid3c_n);

		for (i = 0; i < m - 4; i++) {
			current1_a[i][0] = cur_a[i + 3];
			current1_b[i][0] = cur_b[i + 3];
			current1_c[i][0] = cur_c[i + 3];
			current1_n[i][0] = cur_n[i + 3];

			volt1_a[i][0] = Pa[i + 3];
			volt1_b[i][0] = Pb[i + 3];
			volt1_c[i][0] = Pc[i + 3];
			volt1_n[i][0] = Pn[i + 3];
		}
		MatrixMultiply(mid3_a, current1_a, 3, m - 4, 1, coefficient_a);

		/*printf("coefficient_a= \n");
		output(coefficient_a, m, m);
		printf("\n");*/

		MatrixMultiply(mid3_b, current1_b, 3, m - 4, 1, coefficient_b);
		MatrixMultiply(mid3_c, current1_c, 3, m - 4, 1, coefficient_c);
		MatrixMultiply(mid3_n, current1_n, 3, m - 4, 1, coefficient_n);


		MatrixMultiply(mid3c_a, volt1_a, 3, m - 4, 1, coefficientv_a);
		MatrixMultiply(mid3c_b, volt1_b, 3, m - 4, 1, coefficientv_b);
		MatrixMultiply(mid3c_c, volt1_c, 3, m - 4, 1, coefficientv_c);
		MatrixMultiply(mid3c_n, volt1_n, 3, m - 4, 1, coefficientv_n);
		/*Find predicted samples
			I(m)=a1*I(m-1)+a2*I(m-2)+a3*I(m-3)
			for third-order linear predictor*/
		for (i = m; i < k; i++) {
			current_known_a[i] = current_known_a[i - 1] * coefficient_a[0][0] + current_known_a[i - 2] * coefficient_a[1][0] + current_known_a[i - 3] * coefficient_a[2][0];
			current_known_b[i] = current_known_b[i - 1] * coefficient_b[0][0] + current_known_b[i - 2] * coefficient_b[1][0] + current_known_b[i - 3] * coefficient_b[2][0];
			current_known_c[i] = current_known_c[i - 1] * coefficient_c[0][0] + current_known_c[i - 2] * coefficient_c[1][0] + current_known_c[i - 3] * coefficient_c[2][0];
			current_known_n[i] = current_known_n[i - 1] * coefficient_n[0][0] + current_known_n[i - 2] * coefficient_n[1][0] + current_known_n[i - 3] * coefficient_n[2][0];

			volt_known_a[i] = volt_known_a[i - 1] * coefficientv_a[0][0] + volt_known_a[i - 2] * coefficientv_a[1][0] + volt_known_a[i - 3] * coefficientv_a[2][0];
			volt_known_b[i] = volt_known_b[i - 1] * coefficientv_b[0][0] + volt_known_b[i - 2] * coefficientv_b[1][0] + volt_known_b[i - 3] * coefficientv_b[2][0];
			volt_known_c[i] = volt_known_c[i - 1] * coefficientv_c[0][0] + volt_known_c[i - 2] * coefficientv_c[1][0] + volt_known_c[i - 3] * coefficientv_c[2][0];
			volt_known_n[i] = volt_known_n[i - 1] * coefficientv_n[0][0] + volt_known_n[i - 2] * coefficientv_n[1][0] + volt_known_n[i - 3] * coefficientv_n[2][0];
		}

		/*Find D parameter
			D=(I_measured-I_predicted) (for last 100% samples of 1 cycle)*/
		current_measured_a = 0;
		current_predicted_a = 0;

		current_measured_b = 0;
		current_predicted_b = 0;

		current_measured_c = 0;
		current_predicted_c = 0;

		current_measured_n = 0;
		current_predicted_n = 0;


		volt_measured_a = 0;
		volt_predicted_a = 0;

		volt_measured_b = 0;
		volt_predicted_b = 0;

		volt_measured_c = 0;
		volt_predicted_c = 0;

		volt_measured_n = 0;
		volt_predicted_n = 0;

		for (i = m; i < k; i++) {
			current_measured_a += cur_a[i];
			current_predicted_a += current_known_a[i];

			current_measured_b += cur_b[i];
			current_predicted_b += current_known_b[i];

			current_measured_c += cur_c[i];
			current_predicted_c += current_known_c[i];

			current_measured_n += cur_n[i];
			current_predicted_n += current_known_n[i];


			volt_measured_a += Pa[i];
			volt_predicted_a += volt_known_a[i];

			volt_measured_b += Pb[i];
			volt_predicted_b += volt_known_b[i];

			volt_measured_c += Pc[i];
			volt_predicted_c += volt_known_c[i];

			volt_measured_n += Pn[i];
			volt_predicted_n += volt_known_n[i];
		}

		Da[count] = fabs(current_measured_a - current_predicted_a);
		Db[count] = fabs(current_measured_b - current_predicted_b);
		Dc[count] = fabs(current_measured_c - current_predicted_c);
		Dn[count] = fabs(current_measured_n - current_predicted_n);

		Dca[count] = fabs(volt_measured_a - volt_predicted_a);
		Dcb[count] = fabs(volt_measured_b - volt_predicted_b);
		Dcc[count] = fabs(volt_measured_c - volt_predicted_c);
		Dcn[count] = fabs(volt_measured_n - volt_predicted_n);

		
		/*Declare if found occur
		if ((PAS_a > 1) || (PAS_b > 1) || (PAS_c > 1)) {
			printf("Faults");
			return 0;
		}
		else { return 1; }*/
	}
}
/*printf("PAS_amax = %lf\n", PAS_amax);
printf("PAS_amin = %lf\n", PAS_amin);

printf("PAS_bmax = %lf\n", PAS_bmax);
printf("PAS_bmin = %lf\n", PAS_bmin);

printf("PAS_cmax = %lf\n", PAS_cmax);
printf("PAS_cmin = %lf\n", PAS_cmin);

printf("count = %d, cnt = %d \n",count,cnt);*/


end = clock();  // lấy thời gian sau khi thực hiện 
time_use = (double)(end - start) / CLOCKS_PER_SEC;    //Tính thời gian sử dụng

printf("Time used = %lf\n", time_use);

FILE* fp = NULL;

fopen_s(&fp, "D:\\PV\\AC faults Fluke\\Test\\AC\\P20\\P20-AG7.csv", "w");

// Ghi dữ liệu theo định dạng chỉ định vào file
for (int i = 0; i < 401; i++) {
	fprintf(fp, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf \n", 
			PAS_a[i], PAS_b[i], PAS_c[i], PAS_n[i], Da[i], Db[i], Dc[i], Dn[i],
			PASc_a[i], PASc_b[i], PASc_c[i], PASc_n[i], Dca[i], Dcb[i], Dcc[i], Dcn[i]);
}

fclose(fp);

}





