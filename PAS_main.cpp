#define _CRT_SECURE_NO_DEPRECATE
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#define pi M_PI


/*Declaration for input of Voltage*/
double Pa[200], Pn[200];
double Ma[200][200], Mn[200][200];
double Na[200][200], Nn[200][200];
double Y_paa[200][200], Y_pan[200][200];
double Y_pra[200][200], Y_prn[200][200];
double V_paa[200][200], V_pan[200][200];
double V_pra[200][200], V_prn[200][200];
double Phi_a, Phi_n;
double PAS_a[500], PAS_n[500];
double PAS_amax, PAS_nmax;
double PAS_amin, PAS_nmin;
double dot, eucl;

/*Declaration for input of Current*/
double Mca[200][200], Mcn[200][200];
double Nca[200][200], Ncn[200][200];
double Yc_paa[200][200], Yc_pan[200][200];
double Yc_pra[200][200], Yc_prn[200][200];
double Vc_paa[200][200], Vc_pan[200][200];
double Vc_pra[200][200], Vc_prn[200][200];
double Phic_a, Phic_n;
double PASc_a[500], PASc_n[500];

// declare for D parameter
double cur_a[200], cur_b[200], cur_c[200], cur_n[200];

// declare for other elementary
double S1[200][200], S2[200][200];
double mid1[200][200];
double mid2[200][200];
double c1[200][200];

double determinant(double [][25], double);

void cofactor(double [][25], double);

void transpose(double [][25], double [][25], double);

double f, w, Ts;

int i, j, q, cnt, m, k, alpha, sample, count, length;

// Multiply two matrix
void MatrixMultiply(double a[][200], double b[][200], int row1, int col1, int col2, double c[][200])
{

	int i, j, q;
	for (i = 0; i < row1; i++)
	{
		for (j = 0; j < col2; j++)
		{
			c[i][j] = 0;
			for (q = 0; q < col1; q++)
			{
				c[i][j] += (a[i][q] * b[q][j]);
			}
		}
	}
}

// Dot product of two vector
double DotProduct(double a[][200], double b[][200], int k)
{
	int i;
	double c = 0;
	for (i = 0; i < k; i++)
	{
		c += a[i][0] * b[i][0];
	}
	return c;
}

// Euclidian length of vector
double Euclidian(double a[][200], double b[][200], int k)
{
	int i;
	double c = 0;
	double d = 0;
	double e = 0;
	for (i = 0; i < k; i++)
	{
		c += pow(a[i][0], 2);
		d += pow(b[i][0], 2);
	}

	e = sqrt(c) * sqrt(d);
	return e;
}

/*For calculating Determinant of the Matrix */
/*
double determinant(double a[][200], double k)
{
	double s = 1, det = 0, b[200][200];
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
*/
/*
void InvertMatrix(double num[200][200], double f, double inverse[][200])
{
	double b[200][200], fac[200][200], d;
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
*/

    void InvertMatrix()
    {
      double a[25][25], index, d;
      int i, j;
      index = k;
      for (i = 0;i < index; i++)
        {
         for (j = 0;j < index; j++)
           {
            a[i][j]= mid1[i][j];
            }
        }
      d = determinant(a, index);
      if (d == 0)
       printf("\nInverse of Entered Matrix is not possible\n");
      else
       cofactor(a, index);
    }
    /*For calculating Determinant of the Matrix */
    double determinant(double a[25][25], double index)
    {
      double s = 1, det = 0, b[25][25];
      int i, j, m, n, c;
      if (index == 1)
        {
         return (a[0][0]);
        }
      else
        {
         det = 0;
         for (c = 0; c < index; c++)
           {
            m = 0;
            n = 0;
            for (i = 0;i < index; i++)
              {
                for (j = 0 ;j < index; j++)
                  {
                    b[i][j] = 0;
                    if (i != 0 && j != c)
                     {
                       b[m][n] = a[i][j];
                       if (n < (index - 2))
                        n++;
                       else
                        {
                         n = 0;
                         m++;
                         }
                       }
                   }
                 }
              det = det + s * (a[0][c] * determinant(b, index - 1));
              s = -1 * s;
              }
        }
        return (det);
    }
    void cofactor(double num[25][25], double f)
    {
     double b[25][25], fac[25][25];
     int p, q, m, n, i, j;
     for (q = 0;q < f; q++)
     {
       for (p = 0;p < f; p++)
        {
         m = 0;
         n = 0;
         for (i = 0;i < f; i++)
         {
           for (j = 0;j < f; j++)
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
      transpose(num, fac, f);
    }
    /*Finding transpose of matrix*/ 
    void transpose(double num[25][25], double fac[25][25], double r)
    {
      int i, j;
      double b[25][25], inverse[25][25], d;
      for (i = 0;i < r; i++)
        {
         for (j = 0;j < r; j++)
           {
             b[i][j] = fac[j][i];
            }
        }
      d = determinant(num, r);
      for (i = 0;i < r; i++)
        {
         for (j = 0;j < r; j++)
           {
            inverse[i][j] = b[i][j] / d;
            }
        }
     
       for (i = 0;i < r; i++)
        {
         for (j = 0;j < r; j++)
           {
            mid2[i][j] = inverse[i][j];
            }
         }
    }
// input for trix H
// Try thissssssssssssssssss................................................
void inputH(double H1[][200], double H2[][200], double current[200], int m)
{
	double temp;
	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < m - 4; i++)
		{
			H1[i][j] = current[i + j];
		}
	}

	for (int i = 0; i < m - 4; i++)
	{
		temp = H1[i][0];
		H1[i][0] = H1[i][2];
		H1[i][2] = temp;
	}

	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < m - 4; i++)
		{
			H2[j][i] = H1[i][j];
		}
	}
}
//////////////////////////////////////////////////

typedef struct
{
	double voltage_a;
	double voltage_n;

	double current_a;
	double current_n;
} Signal;

int main(void)
{

	clock_t start, end; // Khai báo biến thời gian
	double time_use;		// Thời gian sử dụng
	start = clock();	// Lấy thời gian trước khi thực hiện thuật toán

	FILE *file1;
	file1 = fopen("data_test_2.txt", "r");

	if (file1 == NULL)
	{
		printf("Error opening file. \n");
		return 1;
	}

	volatile int sampledouble = 201;

	volatile Signal S[201];

	int read = 0;
	int records = 0;
    
	do
	{
		//read = fscanf(file1, "%lf %lf %lf %lf", &S[records].voltage_a, &S[records].voltage_n, &S[records].current_a, &S[records].current_n);
        for (i = 0; i < 200; i++) {
            fscanf(file1, "%lf", &S[i].voltage_a);
            fscanf(file1, "%lf", &S[i].voltage_n);
            fscanf(file1, "%lf", &S[i].current_a);
            fscanf(file1, "%lf", &S[i].current_n);
        }

		printf("\n");
		if (read <= sampledouble)
			records++;

		if (read == sampledouble && (!feof(file1)))
		{
			printf("File format incorrect.\n");
			return 1;
		}
		if (ferror(file1))
		{
			printf("Error reading file.\n");
			return 1;
		}
        
	} while (feof(file1));

	fclose(file1);

	printf("\n%d records read.\n\n", records);

	for (i = 0; i < 200; i++)
	{
		printf("%d %lf %lf %lf %lf ",
			   i,
			   S[i].voltage_a,
			   S[i].voltage_n,

			   S[i].current_a,
			   S[i].current_n

		);
		printf("\n");
	}
    
	// Input frequency of inverter for f
	f = 50;

	// Compute w from f
	w = 2 * pi * 50;

	// Choose Ts for sampling period
	Ts = 0.00526;

	// Choose the amount of samples in one sampling window
	k = 10;

	// Choose the length between two sampling window (0.1*sample<=alpha<=0.25*sample)
	alpha = 6;

	// The amount of samples (include samples in one window plus the length between two window)
	sample = k + alpha;

	count = 0;
	cnt = 0;

    for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < k; j++)
		{
			if (i == 0)
			{
				S1[i][j] = cos(w * (j + 1) * Ts);
			}
			else if (i == 1)
			{
				S1[i][j] = sin(w * (j + 1) * Ts);
			}
			S2[j][i] = S1[i][j];
		}
    }
    //fixed --> Wrong in the dimension of matrix S1 and S2

    printf("Matrix S1:\n");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < k; j++) {
            printf("%lf ", S1[i][j]);
        }
        printf("\n");
    }

	printf("Matrix S2:\n");
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < 2; j++) {
            printf("%lf ", S2[i][j]);
        }
        printf("\n");
    }

    // Compute C1=(S2*S1)^-1*S2

	MatrixMultiply(S2, S1, k, 2, k, mid1);

    	printf("Matrix mid1:\n");
        for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            printf("%lf ", mid1[i][j]);
        }
        printf("\n");
    	}

	InvertMatrix();
	//InvertMatrix(mid1, k, mid2);

	printf("Matrix mid2:\n");
        for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            printf("%lf ", mid2[i][j]);
        }
        printf("\n");
    	}

	MatrixMultiply(mid2, S2, 2, 2, k, c1);
    	//while (count < 400)
	while (count < 200)
	{

		count++;
		if (count < sample)
		{
			Pa[count] = S[count].voltage_a;
			Pn[count] = S[count].voltage_n;

			cur_a[count] = S[count].current_a;
			cur_n[count] = S[count].current_n;
		}
		else
		{
			cnt++;
			for (i = 0; i < sample - 1; i++)
			{
				Pa[i] = Pa[i + 1];
				Pn[i] = Pn[i + 1];

				cur_a[i] = cur_a[i + 1];
				cur_n[i] = cur_n[i + 1];
			}

			Pa[sample - 1] = S[count].voltage_a;
			Pn[sample - 1] = S[count].voltage_n;

			cur_a[sample - 1] = S[count].current_a;
			cur_n[sample - 1] = S[count].current_n;

			for (i = 0; i < k; i++)
			{
				Ma[i][0] = Pa[i];
				Mn[i][0] = Pn[i];

				Na[i][0] = Pa[i + alpha];
				Nn[i][0] = Pn[i + alpha];

				Mca[i][0] = cur_a[i];
				Mcn[i][0] = cur_n[i];

				Nca[i][0] = cur_a[i + alpha];
				Ncn[i][0] = cur_n[i + alpha];
			}

			// Matrix caculation

			MatrixMultiply(c1, Ma, 2, k, 1, Y_paa);
			MatrixMultiply(c1, Mn, 2, k, 1, Y_pan);

			MatrixMultiply(c1, Mca, 2, k, 1, Yc_paa);
			MatrixMultiply(c1, Mcn, 2, k, 1, Yc_pan);

			// Matrix caculation

			MatrixMultiply(c1, Na, 2, k, 1, Y_pra);
			MatrixMultiply(c1, Nn, 2, k, 1, Y_prn);

			MatrixMultiply(c1, Nca, 2, k, 1, Yc_pra);
			MatrixMultiply(c1, Ncn, 2, k, 1, Yc_prn);

			// Matrix caculation

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

			// Find PAS
			PAS_a[count] = fabs((Phi_a * 180 / pi) - ((alpha * 360) / k));

			PAS_n[count] = fabs((Phi_n * 180 / pi) - ((alpha * 360) / k));

			PASc_a[count] = fabs((Phic_a * 180 / pi) - ((alpha * 360) / k));

			PASc_n[count] = fabs((Phic_n * 180 / pi) - ((alpha * 360) / k));
		}
	}

    end = clock();									  // lấy thời gian sau khi thực hiện
	time_use = (double)(end - start) / CLOCKS_PER_SEC; // Tính thời gian sử dụng

	printf("Time used = %lf\n", time_use);

	FILE *fp = NULL;

	fp = fopen("PAS_data_output_3.txt", "w+");

	// Ghi dữ liệu theo định dạng chỉ định vào file
	for (int i = 0; i < 201; i++)
	{
		fprintf(fp, "%lf, %lf, %lf, %lf \n",
				PAS_a[i], PAS_n[i],
				PASc_a[i], PASc_n[i]);
	}

	fclose(fp);
	return 0;
}