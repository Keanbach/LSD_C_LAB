# include <stdio.h>
# include <stdlib.h>
# include <math.h>


//////// Prototypes ////////

float left_rect ( float(*)(float), float, float, int );
float right_rect ( float(*)(float) , float , float , int );
float trapz ( float(*)(float) , float , float , int );
float simpson ( float(*)(float), float , float , int );
float fct1( float);
float fct2( float);

////////////////////////////

// Left Rectangle Method
float left_rect ( float (*FUNCT)(float), float a, float b, int n)
{
	float f = 0;
	float x = a;
	float h = (b-a) / n;

	for(int i = 0; i < n-1; i++)
	{
		f = f + FUNCT(x);
		x = x + h;
	}
	f = f + FUNCT(x);

	return h*f;
}


// Right Rectangle Method
float right_rect ( float (*FUNCT)(float), float a, float b, int n)
{
	float f = 0;
	float x = a;
	float h = (b-a) / n;

	for(int i = 1; i < n; i++)
	{
		f = f + FUNCT(x);
		x = x + h;
	}
	f = f + FUNCT(x);

	return h*f;
}


// Trapezoidal Rule
float trapz ( float (*FUNCT)(float), float a, float b, int n)
{
    float h = (b-a) / n;
    float sigma = 0;
    float x = a;

    for(int i = 1; i < n-1 ; i++)
    {
        sigma = sigma + FUNCT(x);
        x = x + h;
    }
    sigma = sigma + FUNCT(x);

    x = ( FUNCT(a) + FUNCT(b)) / 2;
    
    return (sigma + x) * h;
}


// Simpson Method
float simpson ( float (*FUNCT)(float), float a, float b, int n)
{
    float h = (b-a) / n;
    float sigma = 0;
    float x = a;

    for(int i = 0; i < n-1 ; i++)
    {
        sigma = sigma + FUNCT(x) + FUNCT(x+h) + 4 * FUNCT(x + (h/2));
        x = x + h;
    }
    sigma = sigma + FUNCT(x) + FUNCT(x+h) + 4 * FUNCT(x + (h/2));
    
    return (h/6) * sigma;
}



// Integrable functions
float fct1( float t)
{
	return exp(-(t*t)/2);
}

float fct2( float t)
{
	return t;
}



int main()
{
	float a = 0;
	float b = 1;
	int n = 500;
	
	printf("\n%f\n", left_rect( fct1, a, b, n));
	printf("\n%f\n", right_rect( fct1, a, b, n));
	printf("\n%f\n", trapz( fct1, a, b, n));
	printf("\n%f\n", simpson( fct1, a, b, n));
	
	return 0;
}







