#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "root_finder.h"
#include "flow_variables.h"


double f(double x){
	double C1 = x/g_PL;
	double C3 = x/g_PR;
        double part_one = (g_cL*(C1-1))/(g_gamma*sqrt(g_alpha*(1+(g_beta*C1))));
	double part_two = (g_cR*(C3-1))/(g_gamma*sqrt(g_alpha*(1+(g_beta*C3))));			
	return g_vL - part_one - g_vR - part_two;
}
double f1(double x){
	double C1 = x/g_PL;
	double C3 = x/g_PR;
	return g_vL + (g_cL/(g_alpha*g_gamma))*(1-pow(C1,g_alpha)) - g_vR - (g_cR*(C3-1))/(g_gamma*sqrt(g_alpha*(1+(g_beta*C3))));
}
double f2(double x){
        double C1 = x/g_PL;
	double C3 = x/g_PR;
	return g_vL - (g_cL*(C1-1))/(g_gamma*sqrt(g_alpha*(1+(g_beta*C1)))) - g_vR + (g_cR/(g_alpha*g_gamma))*(1-pow(C3,g_alpha));
}
double* set_elements(double* arr,double p1,double P1,double v1,double p2,double P2,double v2,double VL,double VR){
	arr = (double*) malloc(8);
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
double* riemannSolver(/*double g_pL, double g_PL, double g_vL, double g_pR, double g_PR, double g_vR, double g_gamma*/){
	double* intermediate_states = NULL;
	g_states = (char*)malloc(3);
	g_states[1] = 'C';
	double P1,P2,p1,p2,v1,v2,VL,VR,C1,C3;
	P1 = pow(((g_cL + g_cR + (g_gamma*g_alpha*(g_vL-g_vR)))/((g_cL*pow(g_PL,-g_alpha)) + (g_cR*pow(g_PR,-g_alpha)))),1/g_alpha);
	if(g_PL > P1 && g_PR > P1){
		C1 = P1/g_PL;
		C3 = P1/g_PR;
		p1 = g_pL*pow(C1,1/g_gamma);
		printf("%f\n",pow(C1,1/g_gamma));
		p2 = g_pR*pow(C3,1/g_gamma);
		v1 = g_vR - (g_cR/(g_alpha*g_gamma))*(1-pow(C3,g_alpha));
		v2 = v1;
		P2 = P1;
		printf("%s%f%c%f%c%f%s%f%c%f%c%f%c\n","(p1,P1,v1) (p2,P2,v2): (",p1,',',P1,',',v1,") (",p2,',',P2,',',v2,')');
		g_states[0] = 'R';
		g_states[2] = 'R';
		return set_elements(intermediate_states,p1,P1,v1,p2,P2,v2,MAX_FLOAT,MAX_FLOAT); 

	       	
	}
	else{
		P1 = secant(f,0.0,(((g_PL)<(g_PR))?(g_PL):(g_PR)),100);
		if(g_PL < P1 && g_PR < P1){
			C1 = P1/g_PL;
			C3 = P1/g_PR;
			p1 = g_pL * ((1+ g_beta*C1)/(g_beta + C1));
			p2 = g_pR * ((1+ g_beta*C3)/(g_beta + C3));
			v1 = g_vL -g_cL * (C1 - 1)/(g_gamma*sqrt(g_alpha*(1+(g_beta*C1))));
			VL = g_vL - g_cL*sqrt(g_alpha*(1+g_beta*C1));
			VR = g_vR + g_cR*sqrt(g_alpha*(1+g_beta*C3));
			v2 = v1;
			P2 = P1;
			printf("%s%f%c%f%c%f%s%f%c%f%c%f%c\n","(p1,P1,v1) (p2,P2,v2): (",p1,',',P1,',',v1,") (",p2,',',P2,',',v2,')');
			g_states[0] = 'S';
			g_states[2] = 'S';
			return set_elements(intermediate_states,p1,P1,v1,p2,P2,v2,VL,VR);
		}
		else{
			
			P1 = secant(f1,0.0,(((g_PL)<(g_PR))?(g_PL):(g_PR)),100);
			if(g_PL > P1 && P1 > g_PR){
				C1 = P1/g_PL;
				C3 = P1/g_PR;
				p1 = g_pL*pow(C1,1/g_gamma);
				p2 = g_pR * ((1+ g_beta*C3)/(g_beta + C3));
				v1 = g_vL + (g_cL/(g_alpha*g_gamma))*(1-pow(C1,g_alpha));
				v2 = v1;
				P2 = P1;
				VR = g_vR + g_cR*sqrt(g_alpha*(1+g_beta*C3));
				printf("%s%f%c%f%c%f%s%f%c%f%c%f%c\n","(p1,P1,v1) (p2,P2,v2): (",p1,',',P1,',',v1,") (",p2,',',P2,',',v2,')');
				g_states[0] = 'R';
				g_states[2] = 'S';
				return set_elements(intermediate_states,p1,P1,v1,p2,P2,v2,MAX_FLOAT,VR);

			}
			else{
                        
				P1 = secant(f2,0.0,(((g_PL)<(g_PR))?(g_PL):(g_PR)),100);
				//printf("%f\n",P1);
				if(g_PL < P1 && P1 < g_PR){
					C1 = P1/g_PL;
					C3 = P1/g_PR;
					p1 = g_pL * ((1+ g_beta*C1)/(g_beta + C1));
					p2 =  g_pR*pow(C3,1/g_gamma);
					v1 = g_vR - (g_cR/(g_alpha*g_gamma))*(1-pow(C3,g_alpha));
					v2 = v1;
					P2 = P1;
					VL = g_vL - g_cL*sqrt(g_alpha*(1+g_beta*C1));
					printf("%s%f%c%f%c%f%s%f%c%f%c%f%c\n","(p1,P1,v1) (p2,P2,v2): (",p1,',',P1,',',v1,") (",p2,',',P2,',',v2,')');
					g_states[0] = 'S';
					g_states[2] = 'R';
					return set_elements(intermediate_states,p1,P1,v1,p2,P2,v2,VL,MAX_FLOAT);
				}
			}
		}
	}
	return intermediate_states;
}
