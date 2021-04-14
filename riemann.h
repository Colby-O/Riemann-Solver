#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "root_finder.h"
#include "flow_variables.h"

#define MAX_FLOAT 3.4028234664E+38

#define MAX(a,b) ((a > b) ? a:b)

double two_shock(double P)
{
	double C1 = P/g_PL;
	double C3 = P/g_PR;
        double LHS = (g_cL*(C1-1))/(g_gamma*sqrt(g_alpha*(1+(g_beta*C1))));
	double RHS = (g_cR*(C3-1))/(g_gamma*sqrt(g_alpha*(1+(g_beta*C3))));			
	return g_vL - LHS - g_vR - RHS;
}
double RF_shock(double x)
{
	double C1 = x/g_PL;
	double C3 = x/g_PR;
	return g_vL + (g_cL/(g_alpha*g_gamma))*(1-pow(C1,g_alpha)) - g_vR - (g_cR*(C3-1))/(g_gamma*sqrt(g_alpha*(1+(g_beta*C3))));
}
double shock_RF(double x)
{
        double C1 = x/g_PL;
	double C3 = x/g_PR;
	return g_vL - (g_cL*(C1-1))/(g_gamma*sqrt(g_alpha*(1+(g_beta*C1)))) - g_vR + (g_cR/(g_alpha*g_gamma))*(1-pow(C3,g_alpha));
}
double* set_elements(double* arr,double p1,double P1,double v1,double p2,double P2,double v2,double VL,double VR)
{
	arr = (double*) malloc(8*sizeof(double));
	// checks if there is  enough HEAP memory
	if (arr == NULL)
	{
		printf("ERROR: OUT OF MEMORY!!!\n");
		exit(1);
	}
	arr[0] = p1;
        arr[1] = P1;
        arr[2] = v1;
        arr[3] = p2;
        arr[4] = P2;
        arr[5] = v2;
        arr[6] = VL;
        arr[7] = VR;
	return arr;

}
double* riemannSolver(void)
{
	double* intermediate_states = NULL;
	double P1,p1,p2,v1,VL,VR,C1,C3;
	// Checks for 1,3-RF.
	P1 = pow(((g_cL + g_cR + (g_gamma*g_alpha*(g_vL-g_vR)))/((g_cL*pow(g_PL,-g_alpha)) + (g_cR*pow(g_PR,-g_alpha)))),pow(g_alpha,-1));
	if(g_PL > P1 && g_PR > P1)
	{
		C1 = P1/g_PL;
		C3 = P1/g_PR;
		p1 = g_pL*pow(C1,1/g_gamma);
		p2 = g_pR*pow(C3,1/g_gamma);
		v1 = g_vR - (g_cR/(g_alpha*g_gamma))*(1-pow(C3,g_alpha));
		VL = MAX_FLOAT;
	        VR = MAX_FLOAT; 	
		g_states[0] = 'R';
		g_states[1] = 'R';
		return set_elements(intermediate_states,p1,P1,v1,p2,P1,v1,MAX_FLOAT,MAX_FLOAT); 

	       	
	}
	else
	{
		// Checks for 1,3-shock
		P1 = secant(two_shock,0.0,MAX(g_PL,g_PR),100);
		if(g_PL < P1 && g_PR < P1)
		{
			C1 = P1/g_PL;
			C3 = P1/g_PR;
			p1 = g_pL * ((1+ g_beta*C1)/(g_beta + C1));
			p2 = g_pR * ((1+ g_beta*C3)/(g_beta + C3));
			v1 = g_vL -g_cL * (C1 - 1)/(g_gamma*sqrt(g_alpha*(1+(g_beta*C1))));
			VL = g_vL - g_cL*sqrt(g_alpha*(1+g_beta*C1));
			VR = g_vR + g_cR*sqrt(g_alpha*(1+g_beta*C3));
			g_states[0] = 'S';
			g_states[1] = 'S';
			return set_elements(intermediate_states,p1,P1,v1,p2,P1,v1,VL,VR);
		}
		else
		{
			// checks for 1-RF, 3-shock
			P1 = secant(RF_shock,0.0,MAX(g_PL,g_PR),100);
			if(g_PL > P1 && P1 > g_PR)
			{
				C1 = P1/g_PL;
				C3 = P1/g_PR;
				p1 = g_pL*pow(C1,1/g_gamma);
				p2 = g_pR * ((1+ g_beta*C3)/(g_beta + C3));
				v1 = g_vL + (g_cL/(g_alpha*g_gamma))*(1-pow(C1,g_alpha));
				VR = g_vR + g_cR*sqrt(g_alpha*(1+g_beta*C3));
				VL = MAX_FLOAT; 
				g_states[0] = 'R';
				g_states[1] = 'S';
				return set_elements(intermediate_states,p1,P1,v1,p2,P1,v1,MAX_FLOAT,VR);

			}
			else
			{
                        	// checks for 1-shock, 3-RF
				P1 = secant(shock_RF,0.0,MAX(g_PL,g_PR),100);
				if(g_PL < P1 && P1 < g_PR)
				{
					C1 = P1/g_PL;
					C3 = P1/g_PR;
					p1 = g_pL * ((1+ g_beta*C1)/(g_beta + C1));
					p2 =  g_pR*pow(C3,1/g_gamma);
					v1 = g_vR - (g_cR/(g_alpha*g_gamma))*(1-pow(C3,g_alpha));
					VL = g_vL - g_cL*sqrt(g_alpha*(1+g_beta*C1));
					VR = MAX_FLOAT; 
					g_states[0] = 'S';
					g_states[1] = 'R';
					return set_elements(intermediate_states,p1,P1,v1,p2,P1,v1,VL,MAX_FLOAT);
				}
			}
		}
	}
	printf("ERROR: NO FEASIBLE STATE FOUND!!!");
	exit(2);
}
