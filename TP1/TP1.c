# include <stdio.h>
# include <stdlib.h>
# include <math.h>


/*
TP1 - Solving Linear Systems: { Descent Method, Ascent Method, Gauss Elimination, LU Decomposition, Cholesky Decomposition }

Autor: ABBASSI Yassir
*/

///////////////Prototypes///////////////

void input_v( int, float []);
void display_v( int , float []);

void input_m( int n, float [][n]);
void display_m( int n, float [][n]);

float *descent( int n, float [][n], float []);
float *ascent( int n, float [][n], float []);

void row_reduction ( int n, float [][n], float []);
void lu( int n, float [][n], float [][n]);
void cholesky( int n, float [][n], float [][n], float [][n]);

// Extras
void copy(int n, float [][n], float [][n], float [], float[]);

///////////////////////////////////////

void input_v(int n, float V[])
{
	printf("\n");
	for(int i=0; i < n; i++)
	{
		printf("\t[%d] = ", i);
		scanf("%f", &V[i]);
	}
}
	
void input_m(int n, float M[][n])
{
	printf("\n");
	for(int i=0; i < n; i++)
	{
		for( int j = 0; j < n; j++)
		{
			printf("\t[%d][%d] = ", i, j);
			scanf("%f", &M[i][j]);
		}
	}
}


void display_v( int n, float V[])
{
	printf("\n");
	for(int i=0; i < n; i++)
	{
		printf("\t%f\n",V[i]);
	}
}

void display_m( int n, float M[][n])
{
	printf("\n");
	for(int i=0; i<n; i++)
	{
		for(int j=0; j < n; j++)
		{
			printf("\t%f",M[i][j]);
		}
		printf("\n");
	}
}



float *descent(int n, float A[][n], float B[])
{
	float s;
	float *x;
	
	x = malloc (sizeof(float) * n);
	
	x[0] = B[0] / A[0][0];

	for(int i = 1; i < n; i++)
	{
		s = 0;
		for(int j = 0; j < i; j++)
		{
			s = s + A[i][j] * x[j];
		}
		x[i] = (B[i]-s) / A[i][i];
	}
	return x;
}



float *ascent(int n, float A[][n], float B[])
{
	float s;
	float *x;
	
	x = malloc (sizeof(float) * n);
	
	x[n-1] = B[n-1] / A[n-1][n-1];
	
	for(int i = n-2; i > -1 ; i--)
	{
		s = 0;
		
		for(int j = i+1; j < n; j++)
		{
			s = s + A[i][j] * x[j];
		}
		x[i] = (B[i]-s) / A[i][i];
	}
	return x;
}



void row_reduction (int n, float A[][n], float B[])
{
	float r;
	for(int k = 0; k < n-1; k++)
	{
		for(int i = k+1; i < n; i++)
		{
			r = A[i][k] / A[k][k];
			
			for(int j=0; j < n; j++)
			{
				A[i][j] = A[i][j] - (r * A[k][j]);
			}
			B[i] = B[i] - (r * B[k]);
		}
	}
}



void lu(int n, float U[][n], float L[][n])
{
	for(int i=0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			if( i == j)
			{
				L[i][j] = 1;
				break;
			}
			L[i][j] = 0;
		}
	}
	
	for(int k = 0; k < n-1; k++)
	{
		for(int i = k+1; i < n; i++)
		{
			L[i][k] = U[i][k] / U[k][k];
			
			for(int j=0; j < n; j++)
			{
				U[i][j] = U[i][j] - (L[i][k] * U[k][j]);
			}
		}
	}
}



void cholesky(int n, float A[][n], float C[][n], float CT[][n])
{
	float sigma = 0;
	
	C[0][0] = sqrt(A[0][0]);
	CT[0][0] = C[0][0];
	
	for(int i = 1; i < n; i++)
	{
		C[i][0] = A[i][0] / C[0][0];
		CT[0][i] = C[i][0];
	}
	
	for(int j = 1; j < n; j++)
	{
		
		for(int i = 0; i < j-1; i++)
		{
			C[i][j] = 0;
			CT[j][i] = 0;
		}
		
		sigma = 0;
		for(int i = 0; i < j-1; i++)
		{
			sigma = sigma + C[j][i] * C[j][i];
		}
		
		C[j][j] = sqrt ( A[j][j] - sigma );
		CT[j][j] = C[j][j];
		
		
		for(int i = j+1; i < n; i++)
		{
			sigma = 0;
			for(int t = 0; t < j-1; t++)
			{
				sigma = sigma + (C[i][t] * C[j][t]);
			}
			C[i][j] = (1/C[j][j]) * ( A[i][j] - sigma );
			CT[j][i] = C[i][j];
		}
	}
}



void copy(int n, float MAT[][n], float MAT_sub[][n], float VECT[], float VECT_sub[])
{
	for(int i = 0; i < n; i++)
	{
		VECT_sub[i] = VECT[i];
		for(int j = 0; j < n; j++)
		{
			MAT_sub[i][j] = MAT[i][j];
		}
	}
}



int main()
{
	printf("\n\t\tLinear System Solving Program\n\n");
	
	int n;
	printf("Enter the Dimension of the Linear System: ");
	scanf("%d", &n);
	
	// Declration of input & output variables
	float MAT[n][n], VECT[n];
	float *x = (float*) malloc(n * sizeof(float));
	
	// Setting matrix & vector coefficients
	printf("\nSet the coefficients of the Matrix");
	input_m ( n, MAT );
	
	printf("\nSet the coefficients of the Vector");
	input_v ( n,  VECT );
	
	// Declaration of copy variables
	float A[n][n], B[n];
	copy(n, MAT, A, VECT, B);
	
	
	// Setting Switch loop environment
	int option;
	
	printf("\nSelect your option:\n\t(1) Row Reduction\t(2) LU Decomposition\n\t(3) Cholesky Method\n\tYour Answer: ");
	
	scanf("%d", &option);
	
	switch(option)
	{
		case 1:
		{
			printf("\nYour Matrix");
			display_m(n ,A);
			
			printf("\nYour Vector");
			display_v(n ,B);
			
			printf("\nApplaying Gaussian Elimination...");
			row_reduction(n , A, B);
			
			printf("\nSimplified Matrix");
			display_m(n ,A);
			
			printf("\nResulting Vector");
			display_v(n ,B);
			
			// Solution
			x = ascent(n, A, B);
			printf("\nA X = B solution using only Gauss Elimination:");
			display_v(n ,x);
			break;
		}
		
		
		case 2:
		{
			float L[n][n];
			float *Y = (float*) malloc(n * sizeof(float));
			
			printf("\nYour Matrix");
			display_m(n ,A);
			
			printf("\nYour Vector");
			display_v(n ,B);
			
			printf("\nApplaying LU Decomposition program...");
			lu(n , A, L);
			
			printf("\nLower Matrix:");
			display_m(n ,L);
			
			printf("\nUpper Matrix");
			display_m(n ,A);
			
			// Solution
			Y = descent( n , L , B );
			x = ascent( n , A , Y);
			printf("\nA X = B solution using LU Decomposition Method:");
			display_v(n ,x);
			break;
		}
		
		
		case 3:
		{
			float C[n][n], CT[n][n];
			float *Y = (float*) malloc(n * sizeof(float));
			
			printf("\nYour Matrix");
			display_m(n ,A);
			
			printf("\nYour Vector");
			display_v(n ,B);
			
			printf("\nApplaying Cholesky Decomposition program...");
			cholesky(n, A, C, CT);
			
			printf("\nCholesky Matrix:");
			display_m(n ,C);
			
			printf("\nCholesky Transpose Matrix");
			display_m(n ,CT);
			
			// Solution
			Y = descent( n , C , B );
			x = ascent( n , CT , Y);
			printf("\nA X = B solution using Cholesky Decomposition Method:");
			display_v(n ,x);
			break;
		}
		
		default: ;
			printf("\nUnregistered input, PROGRAM WILL EXIT\n");
	}
	return 0;
}
