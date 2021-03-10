#include <stdio.h>
#include <stdint.h>
#include <math.h>

#define MAX_FLOAT 3.4028234664E+38
// Function declartions
double secant(double(*)(double),double, double, int);
// Function denfitions

/*
 * Name: secant(double(*f)(double),double x1, double x2, int MAX_ITER) 
 * Description: 
 * Parameters:
 * Return:  
 */ 
double secant(double(*f)(double),double x1, double x2, int MAX_ITER)
{
	
	double maxError = 1E-8;
	double err = MAX_FLOAT;
	double f1 = (f(x1) > 0.0)?(f(x1)): -f(x1);
	double f2 = (f(x2) > 0.0)?(f(x2)): -f(x2);
	// pick most recent guess to have the smallest function eval
	if(f(x1) < f(x2)){
		x2 = x2 - x1;
		x1 = x2 + x1;
		x2 = x1 - x2;
	}
	// Iterates untill 1E-8 error is reached or max iterations is exceeded
	for(int i = 0;i < MAX_ITER;++i){
		double temp = x2;
		x2 = x2 - (f(x2)*(x2-x1))/(f(x2)-f(x1));
		x1 = temp;
		err = (x1-x2>0.0)?(x1-x2): -(x1-x2);
		//printf("%f\n",f(i));
		if(err < maxError){return x2;} 
	}
	printf("%s %d %s\n","ERROR: Failed to converge within", MAX_ITER,"iterations");
       	return MAX_FLOAT;	
}
