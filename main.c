#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "riemann.h"

#define NUM_POINTS 206
#define NUM_COMMANDS 14

#define MAX(a,b) ((a > b) ? a:b)
#define MIN(a,b) ((a < b) ? a:b)

// computes p, P , v, and c for a 1-family Refraction Fan (RF)
double* one_family(double u)
{
	double p = g_pL * pow(((2*g_cL + (g_vL - u)*(g_gamma - 1))/(g_cL*(g_gamma + 1))),2/(g_gamma-1));
	double P =  g_PL * pow(((2*g_cL + (g_vL - u)*(g_gamma - 1))/(g_cL*(g_gamma + 1))),(2*g_gamma)/(g_gamma-1));
	double v = ((2*(g_cL + u) + g_vL*(g_gamma - 1))/(g_gamma + 1));
	double c = ((2*g_cL + (g_vL - u)*(g_gamma - 1))/(g_gamma + 1));
	double* vars = (double*)malloc(4*sizeof(double));
	vars[0] = p;
	vars[1] = P;
	vars[2] = v;
	vars[3] = c;
	return vars;
}
// computes p, P, v, c for a 3-family Refraction Fan (RF)
double* three_family(double u)
{
	double p = g_pR * pow(((2*g_cR + (u - g_vR)*(g_gamma - 1))/(g_cR*(g_gamma + 1))),2/(g_gamma-1));
	double P = g_PR * pow(((2*g_cR + (u - g_vR)*(g_gamma - 1))/(g_cR*(g_gamma + 1))),(2*g_gamma)/(g_gamma-1));
	double v = ((2*(u - g_cR) + g_vR*(g_gamma - 1))/(g_gamma + 1));
	double c = ((2*g_cR + (u - g_vR)*(g_gamma - 1))/(g_gamma + 1));
	double* vars = (double*)malloc(4*sizeof(double));
	vars[0] = p;
	vars[1] = P;
	vars[2] = v;
	vars[3] = c;
	return vars;
}

int main(int argc, char**argv)
{
	// Checks for vaild number of arugments
	if(argc != 7)
	{
		printf("ERROR: NOT ENOUGH ARUGMENTS!!!\n");
		exit(3);
	}
	
	// sets arugments to corrsponding left and right states
	sscanf(argv[1],"%lf", &g_pL);
	sscanf(argv[2],"%lf",&g_PL);
	sscanf(argv[3],"%lf",&g_vL);

	sscanf(argv[4],"%lf",&g_pR);
	sscanf(argv[5],"%lf",&g_PR);
	sscanf(argv[6],"%lf",&g_vR);

	// Checks for vaild densities and pressures
	if(g_pL < 0 || g_PL < 0 || g_pR < 0 || g_PR < 0)
	{
		printf("ERROR: INVAID ARUGMENTS (DENSITY AND PRESSURE CANNOT BE NEAGTIVE)!!!\n");
		exit(3);
	}

	// gets solution to the reimeann problem.
	double* intermediate_states = riemannSolver();

	// allocate HEAP memory for x/t, density, pressure, velocity and sound speed data points for plots
	double* xt_vals = (double*)malloc(206*sizeof(double));
	double* p_vals = (double*)malloc(206*sizeof(double));
	double* P_vals = (double*)malloc(206*sizeof(double));
	double* v_vals = (double*)malloc(206*sizeof(double));
	double* c_vals = (double*)malloc(206*sizeof(double));

	// checks if HEAP memory is available
	if (xt_vals == NULL || p_vals == NULL || P_vals == NULL || v_vals == NULL || c_vals ==NULL)
	{
		printf("ERROR: OUT OF MEMORY");
		exit(1);
	}

	// initializes each element of the arrays.
	for (int i = 0; i < 100; ++i)
	{
		xt_vals[i] = 0;
		p_vals[i] = 0;
		P_vals[i] = 0;
		v_vals[i] = 0;
		c_vals[i] = 0;

	}
	// Sets points for the contact discontinuity. 
	xt_vals[102] = intermediate_states[2];
	xt_vals[103] = intermediate_states[2];

	p_vals[102] = intermediate_states[0];
	p_vals[103] = intermediate_states[3];

	P_vals[102] = intermediate_states[1];
	P_vals[103] = intermediate_states[4];

	v_vals[102] = intermediate_states[2];
	v_vals[103] = intermediate_states[5];

	c_vals[102] = sqrt(g_gamma*intermediate_states[1]/intermediate_states[0]);
	c_vals[103] = sqrt(g_gamma*intermediate_states[4]/intermediate_states[3]);

	if (g_states[0] == 'R')
	{
		double lb = g_vL - g_cL; // lower x/t value for the 1-RF
		double ub = intermediate_states[2] - sqrt(g_gamma*intermediate_states[1]/intermediate_states[0]); // upper x/t value for the 1-RF
		double h = (ub - lb)/100; // step size to computes plot points

		xt_vals[0] = lb - 1.0; // first plot point
		p_vals[0] = g_pL; // left Density
		P_vals[0] = g_PL; // left Pressure
		v_vals[0] = g_vL; // left velocity
		c_vals[0] = g_cL; // left sound speed

		for (int i = 0; i < 101; ++i)
		{
			xt_vals[i+1] = lb + h*i;
			double*  RF_vals = one_family(xt_vals[i+1]);
			p_vals[i+1] = RF_vals[0];
			P_vals[i+1] = RF_vals[1];
			v_vals[i+1] = RF_vals[2];
			c_vals[i+1] = RF_vals[3];
			free(RF_vals);
		}
	}
	else 
	{
		// shocks position
		double xt = (intermediate_states[6] > 0)?(intermediate_states[6]): intermediate_states[6];

		xt_vals[0] = xt - 1.0; // first plot point
		p_vals[0] = g_pL; // left Density
		P_vals[0] = g_PL; // left Pressure
		v_vals[0] = g_vL; // left velocity
		c_vals[0] = g_cL; // left sound speed
		
		xt_vals[1] = xt; // second points
		p_vals[1] = g_pL;
		P_vals[1] = g_PL;
		v_vals[1] = g_vL;
		c_vals[1] = g_cL;

		for (int i = 1; i < 101; ++i)
		{
			xt_vals[i+1] = xt;
			p_vals[i+1] = intermediate_states[0];
			P_vals[i+1] = intermediate_states[1];
			v_vals[i+1] = intermediate_states[2];
			c_vals[i+1] = sqrt(g_gamma*intermediate_states[1]/intermediate_states[0]);
		}
	
	}

	if (g_states[1] == 'R')
	{

		double lb = intermediate_states[2] + sqrt(g_gamma*intermediate_states[4]/intermediate_states[3]); // lower x/t value for the 3-RF
		double ub = g_vR + g_cR; // upper x/t value for the 3-RF
		double h = (ub - lb)/100; // step size to compute plot points

		for (int i = 0; i < 101; ++i)
		{
			xt_vals[i+104] = lb + h*i;
			double* RF_vals = three_family(xt_vals[i+104]);
			p_vals[i+104] = RF_vals[0];
			P_vals[i+104] = RF_vals[1];
			v_vals[i+104] = RF_vals[2];
			c_vals[i+104] = RF_vals[3];
			free(RF_vals);

		}
		// sets the final values
		xt_vals[205] = xt_vals[204] + 1.0;
		p_vals[205] = g_pR;
		P_vals[205] = g_PR;
		v_vals[205] = g_vR;
		c_vals[205] = g_cR;
	}
	else 
	{
		// shocks position
		double xt = (intermediate_states[7] > 0)?(intermediate_states[7]): intermediate_states[7];

		for (int i = 0; i < 100; ++i)
		{
			xt_vals[i+104] = xt;
			p_vals[i+104] = intermediate_states[3];
			P_vals[i+104] = intermediate_states[4];
			v_vals[i+104] = intermediate_states[5];
			c_vals[i+104] = sqrt(g_gamma*intermediate_states[4]/intermediate_states[3]);
		}
		// sets the final values
		xt_vals[204] = xt;
		xt_vals[205] = xt + 1.0;

		p_vals[204] = g_pR;
		p_vals[205] = g_pR;

		P_vals[204] = g_PR;
		P_vals[205] = g_PR;

		v_vals[204] = g_vR;
		v_vals[205] = g_vR;

		c_vals[204] = g_cR;
		c_vals[205] = g_cR;
	}
	
	// prints the values of each state
	printf("Left State: \nρₗ: %lf  Pₗ: %lf  vₗ: %lf\n",g_pL,g_PL,g_vL);
	printf("Right State: \nρᵣ: %lf  Pᵣ: %lf  vᵣ: %lf\n",g_pR,g_PR,g_vR);
	printf("Intermediate State 1: \nρ₁: %lf  P₁: %lf  v₁: %lf\n",intermediate_states[0],intermediate_states[1],intermediate_states[2]);
	printf("Intermediate State 2: \nρ₂: %lf  P₂: %lf  v₂: %lf\n",intermediate_states[3],intermediate_states[4],intermediate_states[5]);
	printf("Shock Speeds: \nV₁: %lf  V₂: %lf\n",intermediate_states[6],intermediate_states[7]);

	// gets the max and min values of p,P,v,and c for plottig
	double p_max = MAX(MAX(g_pL,g_pR),MAX(intermediate_states[0],intermediate_states[3]));	
	double p_min = MIN(MIN(g_pL,g_pR),MIN(intermediate_states[0],intermediate_states[3]));

	double P_max = MAX(MAX(g_PL,g_PR),MAX(intermediate_states[1],intermediate_states[4]));
	double P_min = MIN(MIN(g_PL,g_PR),MIN(intermediate_states[1],intermediate_states[4]));

	double v_max = MAX(MAX(g_vL,g_vR),MAX(intermediate_states[2],intermediate_states[5]));
	double v_min = MIN(MIN(g_vL,g_vR),MIN(intermediate_states[2],intermediate_states[5]));

	double c_max = MAX(MAX(g_cL,g_cR),MAX(sqrt(g_gamma*intermediate_states[1]/intermediate_states[0]),sqrt(g_gamma*intermediate_states[4]/intermediate_states[3])));
	double c_min = MIN(MIN(g_cL,g_cR),MIN(sqrt(g_gamma*intermediate_states[1]/intermediate_states[0]),sqrt(g_gamma*intermediate_states[4]/intermediate_states[3])));

	// sets up GNUplots  plotting commands  
	char p_plot[100],P_plot[100],v_plot[100],c_plot[100];

	double factor = 0.1; // factor that adds extra bit onto the y-axis from the max/min values

	sprintf(p_plot,"plot [] [%lf:%lf] 'data.txt' u 1:2 w linespoints pt 0 notitle", p_min - p_min*factor,p_max + p_max*factor);
	sprintf(P_plot,"plot [] [%lf:%lf] 'data.txt' u 1:3 w linespoints pt 0 notitle", P_min - P_min*factor,P_max + P_max*factor);
	sprintf(c_plot,"plot [] [%lf:%lf] 'data.txt' u 1:4 w linespoints pt 0 notitle", c_min - c_min*factor,c_max + c_max*factor);
	sprintf(v_plot,"plot [] [%lf:%lf] 'data.txt' u 1:5 w linespoints pt 0 notitle", v_min - fabs(v_min*factor),v_max + fabs(v_max*factor));

	// Plots the data
        char* commandsForGnuplot[] = {"set title \"Density Plot\"","set xlabel \"x/t\"","set ylabel \"ρ\"","set multiplot layout 2,2",p_plot,"set title \"Pressure Plot\"","set ylabel \"P\"",P_plot,"set title \"Sound Speed Plot\"","set ylabel \"cₛ\"",c_plot,"set title \"Velocity Plot\"","set ylabel \"v\"",v_plot};
        FILE* temp = fopen("data.txt", "w");
        FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
  
        for (int i=0; i < NUM_POINTS; i++)
        {
      		 fprintf(temp, "%lf %lf %lf %lf %lf \n", xt_vals[i],p_vals[i], P_vals[i],c_vals[i],v_vals[i]); //Write the data to a temporary file
      		 fflush(gnuplotPipe);
        }

        for (int i=0; i < NUM_COMMANDS; i++)
        {
   		 fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
		 fflush(gnuplotPipe);
        }
	
	free(intermediate_states);
	return 0;
}
