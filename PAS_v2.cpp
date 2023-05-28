

/**
 * main.c
 */

#include <stdio.h>
#include <math.h>




 //declare 
double Ma[100][100], Mb[100][100], Mc[100][100];
double Na[100][100], Nb[100][100], Nc[100][100];
double Y_paa[100][100], Y_pab[100][100], Y_pac[100][100];
double Y_pra[100][100], Y_prb[100][100], Y_prc[100][100];
double V_paa[100][100], V_pab[100][100], V_pac[100][100];
double V_pra[100][100], V_prb[100][100], V_prc[100][100];
double Phi_a, Phi_b, Phi_c;
double PAS_a, PAS_b, PAS_c;
double dot, eucl;
int i, j, q, cnt, m, k, alpha, sample;


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
		c += (pow(a[i][0],2));
		d += (pow(b[i][0],2));
	}
	e = sqrt(c*d);
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


//print output (just for test the algorithm)
void output(double a[][100], int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			printf("%lf ", a[i][j]);
		}
		printf("\n");
	}
}

int main(void) {
	
	//Input frequency of inverter for f
	double f = 50;

	//Compute w from f
	double w = 2 * pi * 50;

	//Choose Ts for sampling period
	double Ts = 0.01;

	//Choose the amount of samples in one sampling window
	k = 20;

	//Choose the length between two sampling window (0.1*sample<=alpha<=0.25*sample)
	alpha = 5;

	//The amount of samples (include samples in one window plus the length between two window)
	sample = k + alpha;

	/*Compute constant matrix S1, S2
		S1=[cos(wTs)	sin(wTs)
			cos(w2Ts)	sin(w2Ts)
					...
			cos(wkTs)	sin(wkTs)]

		S2=[cos(wTs)	cos(w2Ts)	...	cos(wkTs)
			sin(wTs)	sin(w2Ts)	...	sin(wkTs)]
	*/

	double S1[100][100], S2[100][100];

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
	printf("S1= \n");
	output(S1, k, 2);
	printf("\n");

	printf("S2= \n");
	output(S2, 2, k);
	printf("\n");

	// Compute C1=(S2*S1)^-1*S2
	double mid1[100][100];
	double mid2[100][100];
	double c1[100][100];
	MatrixMultiply(S2, S1, 2, k, 2, mid1);

	printf("mid1= \n");
	output(mid1, 2, 2);
	printf("\n");

	InvertMatrix(mid1, 2, mid2);
	printf("mid2= \n");
	output(mid2, 2, 2);
	printf("\n");

	MatrixMultiply(mid2, S2, 2, 2, k, c1);
	printf("c1= \n");
	output(c1, 2, k);
	printf("\n");


	while (1) {
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
		for (int i = 0; i < k; i++) {
			Ma[i][0] = 5*cos(w*(i+8)*Ts+100);
			//Mb[i][0] = (double)i + 1;
			//Mc[i][0] = (double)i - 1;
		}
		printf("Ma = \n");
		output(Ma, k, 2);
		printf("\n");

		for (int i = 0; i < k; i++) {
			Na[i][0] = 5*cos(w*(i+alpha+8)*Ts+100);
			//Nb[i][0] = (double)i + 2;
			//Nc[i][0] = (double)i + 3;
		}
		printf("Na = \n");
		output(Na, k, 2);
		printf("\n");
		/*Compute Y_paa, Y_pab, Y_pac from S1, S2 and Ma, Mb, Mc
			Y_paa=c1*Ma
			Y_pab=c1*Mb
			Y_pac=c1*Mc
		*/
		//Matrix caculation

		MatrixMultiply(c1, Ma, 2, k, 1, Y_paa);
		//MatrixMultiply(c1, Mb, 2, k, 1, Y_pab);
		//MatrixMultiply(c1, Mc, 2, k, 1, Y_pac);
		printf("Y_paa = \n");
		output(Y_paa, 2, 1);
		printf("\n");

		/*Compute Y_pra, Y_prb, Y_prc from S1, S2 and Na, Nb, Nc
			Y_paa=c1*Na
			Y_pab=c1*Nb
			Y_pac=c1*Nc
		*/
		//Matrix caculation

		MatrixMultiply(c1, Na, 2, k, 1, Y_pra);
		//MatrixMultiply(c1, Nb, 2, k, 1, Y_prb);
		//MatrixMultiply(c1, Nc, 2, k, 1, Y_prc);

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
		//MatrixMultiply(S1, Y_pab, 2, k, 1, V_pab);
		//MatrixMultiply(S1, Y_pac, 2, k, 1, V_pac);
		printf("V_paa = \n");
		output(V_paa, k, 1);
		printf("\n");

		MatrixMultiply(S1, Y_pra, k, 2, 1, V_pra);
		//MatrixMultiply(S1, Y_prb, 2, k, 1, V_prb);
		//MatrixMultiply(S1, Y_prc, 2, k, 1, V_prc);
		printf("V_pra = \n");
		output(V_pra, k, 1);
		printf("\n");


		/*Find the phase difference between the two vector in radian
			Phi=arcos((V_pr.V_pa)/(abs(V_pr)*abs(V_pa))
		*/
		//Matrix caculation
		
		dot=DotProduct(V_pra, V_paa, k);
		
		printf("%lf \n", dot);

		eucl=Euclidian(V_pra, V_paa, k);
		printf("Euclidian = %lf\n", eucl);
		Phi_a = dot / eucl;
		Phi_a = acos(Phi_a);
		//Phi_b = acos(DotProduct(V_prb, V_pab, k) / (Euclidian(V_prb,k) * Euclidian(V_pab,k)));
		//Phi_c = acos(DotProduct(V_prc, V_pac, k) / (Euclidian(V_prc,k) * Euclidian(V_pac,k)));

		printf("Phi_a = %lf\n", Phi_a*180/pi);

		//printf("Phi_b = %lf\n", Phi_b);
		//printf("Phi_c = %lf\n", Phi_c);
		//Find PAS
		PAS_a = fabs((Phi_a * 180 / pi) - ((alpha * 360) / k));
		//PAS_b = fabs(Phi_b * 180 / pi - (alpha * 360) / k);
		//PAS_c = fabs(Phi_c * 180 / pi - (alpha * 360) / k);

		printf("PAS_a = %lf\n", PAS_a);
		//printf("PAS_b = %lf\n", PAS_b);
		//printf("PAS_c = %lf\n", PAS_c);

		break;

		/*Declare if found occur
		if ((PAS_a > 1) || (PAS_b > 1) || (PAS_c > 1)) {
			printf("Faults");
			return 0;
		}
		else { return 1; }*/
	}
}




