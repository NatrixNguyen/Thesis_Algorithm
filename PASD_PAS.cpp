﻿#define _CRT_SECURE_NO_DEPRECATE
//#include <conio.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>




 //declare for PAS
float Pa[100], Pb[100], Pc[100], Pn[100];
float Ma[100][100], Mb[100][100], Mc[100][100], Mn[100][100];
float Na[100][100], Nb[100][100], Nc[100][100], Nn[100][100];
float Y_paa[100][100], Y_pab[100][100], Y_pac[100][100], Y_pan[100][100];
float Y_pra[100][100], Y_prb[100][100], Y_prc[100][100], Y_prn[100][100];
float V_paa[100][100], V_pab[100][100], V_pac[100][100], V_pan[100][100];
float V_pra[100][100], V_prb[100][100], V_prc[100][100], V_prn[100][100];
float Phi_a, Phi_b, Phi_c, Phi_n;
float PAS_a[500], PAS_b[500], PAS_c[500], PAS_n[500];
float PAS_amax, PAS_bmax, PAS_cmax, PAS_nmax;
float PAS_amin, PAS_bmin, PAS_cmin, PAS_nmin;
float dot, eucl;


float Mca[100][100], Mcb[100][100], Mcc[100][100], Mcn[100][100];
float Nca[100][100], Ncb[100][100], Ncc[100][100], Ncn[100][100];
float Yc_paa[100][100], Yc_pab[100][100], Yc_pac[100][100], Yc_pan[100][100];
float Yc_pra[100][100], Yc_prb[100][100], Yc_prc[100][100], Yc_prn[100][100];
float Vc_paa[100][100], Vc_pab[100][100], Vc_pac[100][100], Vc_pan[100][100];
float Vc_pra[100][100], Vc_prb[100][100], Vc_prc[100][100], Vc_prn[100][100];
float Phic_a, Phic_b, Phic_c, Phic_n;
float PASc_a[500], PASc_b[500], PASc_c[500], PASc_n[500];



//declare for D parameter
float cur_a[100], cur_b[100], cur_c[100], cur_n[100];



//declare for other elementary
float S1[100][100], S2[100][100];
float mid1[100][100];
float mid2[100][100];
float c1[100][100];

float f, w, Ts;

int i, j, q, cnt, m, k, alpha, sample, count, length;;

const float pi = 22.0 / 7.0;





//Multiply two matrix
void MatrixMultiply(float a[][100], float b[][100], int row1, int col1, int col2, float c[][100]) {

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
float DotProduct(float a[][100], float b[][100], int k) {
	int i;
	float c = 0;
	for (i = 0; i < k; i++) {
		c += a[i][0] * b[i][0];
	}
	return c;
}

//Euclidian length of vector
float Euclidian(float a[][100], float b[][100], int k) {
	int i;
	float c = 0;
	float d = 0;
	float e = 0;
	for (i = 0; i < k; i++) {
		c += pow(a[i][0], 2);
		d += pow(b[i][0], 2);
	}

	e = sqrt(c) * sqrt(d);
	return e;
}

/*For calculating Determinant of the Matrix */
float determinant(float a[][100], float k)
{
	float s = 1, det = 0, b[100][100];
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

void InvertMatrix(float num[100][100], float f, float inverse[][100])
{
	float b[100][100], fac[100][100], d;
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
void inputH(float H1[][100], float H2[][100], float current[100], int m) {
	float temp;
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

typedef struct
{
	float voltage_a;
	//float voltage_b;
	//float voltage_c;
	float voltage_n;

	float current_a;
	//float current_b;
	//float current_c;
	float current_n;

	//float voltage_dc;
	//float current_dc;
	//float voltage_fa;
	//float voltage_fb;
	//float voltage_fc;
	//float voltage_fn;

	//float current_fa;
	//float current_fb;
	//float current_fc;
	//float current_fn;
} Signal;


int main(void) {

	clock_t start, end;   // Khai báo biến thời gian
	float time_use;      // Thời gian sử dụng
	start = clock();     // Lấy thời gian trước khi thực hiện thuật toán

	FILE* file1;
	file1 = fopen("data_test.txt", "r");

	if (file1 == NULL) {
		printf("Error opening file. \n");
		return 1;
	}

	volatile int samplefloat = 401;

	volatile Signal S[401];


	int read = 0;
	int records = 0;
	do {
		read = fscanf(file1, "%6f %6f %6f %6f",
			&S[records].voltage_a,
			//&S[records].voltage_b,
			//&S[records].voltage_c,
			&S[records].voltage_n,

			&S[records].current_a,
			//&S[records].current_b,
			//&S[records].current_c,
			&S[records].current_n

			//&S[records].voltage_dc,
			//&S[records].current_dc,

			//&S[records].voltage_fa,
			//&S[records].voltage_fb,
			//&S[records].voltage_fc,
			//&S[records].voltage_fn,

			//&S[records].current_fa,
			//&S[records].current_fb,
			//&S[records].current_fc,
			//&S[records].current_fn
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


	for (i = 0; i < 400; i++)
	{
		printf("%d %6f %6f %6f %6f ",
			i,
			S[i].voltage_a,
			//S[i].voltage_b,
			//S[i].voltage_c,
			S[i].voltage_n,

			S[i].current_a,
			//S[i].current_b,
			//S[i].current_c,
			S[i].current_n

		);
		printf("\n");
	}
	

	
	//Input frequency of inverter for f
	f = 50;

	//Compute w from f
	w = 2 * pi * 50;

	//Choose Ts for sampling period
	Ts = 0.00526;

	//Choose the amount of samples in one sampling window
	k = 10;

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


	InvertMatrix(mid1, 2, mid2);


	MatrixMultiply(mid2, S2, 2, 2, k, c1);



	while (count<400) {

		count++;
		if (count < sample) {
			Pa[count] = S[count].voltage_a;
			Pn[count] = S[count].voltage_n;

			cur_a[count] = S[count].current_a;
			cur_n[count] = S[count].current_n;

		}
		else {
			cnt++;
			for (i = 0; i < sample - 1; i++) {
				Pa[i] = Pa[i + 1];
				Pn[i] = Pn[i + 1];

				cur_a[i] = cur_a[i + 1];
				cur_n[i] = cur_n[i + 1];

			}

			Pa[sample - 1] = S[count].voltage_a;
			Pn[sample - 1] = S[count].voltage_n;

			cur_a[sample - 1] = S[count].current_a;
			cur_n[sample - 1] = S[count].current_n;

			for (i = 0; i < k; i++) {
				Ma[i][0] = Pa[i];
				Mn[i][0] = Pn[i];

				Na[i][0] = Pa[i + alpha];
				Nn[i][0] = Pn[i + alpha];

				Mca[i][0] = cur_a[i];
				Mcn[i][0] = cur_n[i];

				Nca[i][0] = cur_a[i + alpha];
				Ncn[i][0] = cur_n[i + alpha];
			}
			
		//Matrix caculation

		MatrixMultiply(c1, Ma, 2, k, 1, Y_paa);
		MatrixMultiply(c1, Mn, 2, k, 1, Y_pan);

		MatrixMultiply(c1, Mca, 2, k, 1, Yc_paa);
		MatrixMultiply(c1, Mcn, 2, k, 1, Yc_pan);

		//Matrix caculation

		MatrixMultiply(c1, Na, 2, k, 1, Y_pra);
		MatrixMultiply(c1, Nn, 2, k, 1, Y_prn);

		MatrixMultiply(c1, Nca, 2, k, 1, Yc_pra);
		MatrixMultiply(c1, Ncn, 2, k, 1, Yc_prn);

		//Matrix caculation

		MatrixMultiply(S1, Y_paa, k, 2, 1, V_paa);

		MatrixMultiply(S1, Y_pan, k, 2, 1, V_pan);

		MatrixMultiply(S1, Y_pra, k, 2, 1, V_pra);

		MatrixMultiply(S1, Y_prn, k, 2, 1, V_prn);


		MatrixMultiply(S1, Yc_paa, k, 2, 1, Vc_paa);

		MatrixMultiply(S1, Yc_pan, k, 2, 1, Vc_pan);

		MatrixMultiply(S1, Yc_pra, k, 2, 1, Vc_pra);

		MatrixMultiply(S1, Yc_prn, k, 2, 1, Vc_prn);


		/*Find the phase difference between the two vector in radian
			Phi=arcos((V_pr.V_pa)/(abs(V_pr)*abs(V_pa))
		*/
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


		//Find PAS
		PAS_a[count] = fabs((Phi_a * 180 / pi) - ((alpha * 360) / k));

		PAS_n[count] = fabs((Phi_n * 180 / pi) - ((alpha * 360) / k));

		PASc_a[count] = fabs((Phic_a * 180 / pi) - ((alpha * 360) / k));

		PASc_n[count] = fabs((Phic_n * 180 / pi) - ((alpha * 360) / k));


}
}

end = clock();  // lấy thời gian sau khi thực hiện 
time_use = (float)(end - start) / CLOCKS_PER_SEC;    //Tính thời gian sử dụng

printf("Time used = %6f\n", time_use);

FILE* fp = NULL;

fp = fopen("Output_PAS.txt", "w+");

// Ghi dữ liệu theo định dạng chỉ định vào file
for (int i = 0; i < 401; i++) {
	fprintf(fp, "%6f, %6f, %6f, %6f \n", 
			PAS_a[i], PAS_n[i],
			PASc_a[i], PASc_n[i]);
}

fclose(fp);
return 0;
}





