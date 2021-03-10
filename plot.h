#include <stdio.h>
#include <math.h>
#include "cursor.h"

#ifndef _FLOW_VARIABLES_H_
#define _FLOW_VARIAVLES_H_
#endif

#define MIN(a,b) ((a < b) ? a : b)
#define MAX(a,b) ((a > b) ? a : b) 

void init_grid(double*);
void plot_line();
void plot_jump();
void plot_lienar();
void plot_shock();
void plot_RF();
void plot_CD();
void plot_riemann(double*);

double map(double x, double in_min, double in_max, double out_min, double out_max){
	return (x-in_min)*(out_max-out_min)/(in_max-in_min) + out_min;
}

void init_grids(double* vars){
	cursor_go_home_clear();
	int h = get_height()/2 - 1;
	int w = get_width()/2 - 10;
	int xzoom = 1;
	int yzoom = 1;
	int x_margin = 10;
	int y_margin = 3;
	int title_h = 2;
	double p_max = MAX(MAX(g_pL,g_pR),MAX(vars[0],vars[3]));
	double p_min = MIN(MIN(g_pL,g_pR),MIN(vars[0],vars[3]));
	double P_max = MAX(MAX(g_PL,g_PR),MAX(vars[1],vars[4]));
	double P_min = MIN(MIN(g_PL,g_PR),MIN(vars[1],vars[4]));
	double v_max = MAX(MAX(g_vL,g_vR),MAX(vars[2],vars[5]));
	//printf("%f",MAX(MAX(g_vL,g_vR),MAX(vars[2],vars[5])));
	double v_min = MIN(MIN(g_vL,g_vR),MIN(vars[2],vars[5]));
	double c_max = MAX(MAX(sqrt(g_gamma*g_PL/g_pL),sqrt(g_gamma*g_PR/g_pR)),MAX(sqrt(g_gamma*vars[1]/vars[0]),sqrt(g_gamma*vars[4]/vars[3])));
        double c_min = MIN(MIN(sqrt(g_gamma*g_PL/g_pL),sqrt(g_gamma*g_PR/g_pR)),MIN(sqrt(g_gamma*vars[1]/vars[0]),sqrt(g_gamma*vars[4]/vars[3])));
	//printf("%f	%f",c_min,c_max);
        // Plot titles	
	cursor_goto_pos(title_h,w/2);
	printf("Density Solution");
        cursor_goto_pos(title_h,w/2+w+x_margin);
	printf("Pressure Solution");
        cursor_goto_pos(h+title_h,w/2);
	printf("Speed Of Sound Solution");
        cursor_goto_pos(h+title_h,w/2+w+x_margin);
	printf("Velocity Solution");
	// xlabels
	cursor_goto_pos(h+h+4,w/2+x_margin);
	printf("x/t");
	cursor_goto_pos(h+h+4,w/2+w+2*x_margin);
	printf("x/t");
	// ylabels
        cursor_goto_pos(h/2,0);
	printf("ρ");
	cursor_goto_pos(h/2,w+x_margin);
	printf("P");
        cursor_goto_pos(h/2+h,0);
	printf("c");
	cursor_goto_pos(h/2+h,w+x_margin);
	printf("v");
	for(int i = h - 1; y_margin < i;--i){
		cursor_goto_pos(i,x_margin);
		printf("|");
		cursor_goto_pos(i,w+8+x_margin);
		printf("|");
                cursor_goto_pos(i+h,w+8+x_margin);
		printf("|");
                cursor_goto_pos(i+h,x_margin);
		printf("|");
		if(i%4 == 0){
			// NOTE: upper bounds on the map are wrong wont work all the time
			double temp = map(h-i-2,-1,h-y_margin-3,p_min,p_max);
                	cursor_goto_pos(i,x_margin-5);
	        	printf("%.2f-",temp);
			temp = map(h-i-2,-1,h-y_margin-3,P_min,P_max);
                        cursor_goto_pos(i,w+8+x_margin-5);
			printf("%.2f-",temp);
			// FIXED both bound they will never work!!!
			temp = map(h-i-2,-1,h-y_margin-3,c_min,c_max);
			cursor_goto_pos(i+h,x_margin-5);
			printf("%.2f-",temp);
			// SAME FIX!!!
                        temp = map(h-i-2,-1,h-y_margin-3,v_min,v_max);
			cursor_goto_pos(i+h,w+8+x_margin-5);
			printf("%.2f-",temp);
		}
	}
	
	for(int i = 0; i < w;++i){
		cursor_goto_pos(h,i+x_margin+1);
		printf("-");
                cursor_goto_pos(h,i+w+9+x_margin);
		printf("-");
                cursor_goto_pos(h+h,i+x_margin+1);
		printf("-");
                cursor_goto_pos(h+h,i+w+9+x_margin);
		printf("-");
		if(i%5 == 0){
			double temp = map(i+1,1,w,-2,2);
			cursor_goto_pos(h+h+1,i+x_margin+1);
			printf("|");
                	cursor_goto_pos(h+h+2,i+x_margin+1);
			printf("%.1f",temp);

			cursor_goto_pos(h+h+1,i+2*x_margin+w-1);
			printf("|");
			cursor_goto_pos(h+h+2,i+2*x_margin+w-1);
			printf("%.1f",temp);
		}
	}

	
	/*for(int i = h-1; 1 < i;--i){
		for(int j = 1; j < w;++j)
		if(i == h - j){
			cursor_goto_pos(i+2,j+6);
			printf("%c",'@');
		}
	
	}
	
	for(int i = 42; i < 80; ++i){
		cursor_goto_pos(h-42,i+8);
		printf("%c",'@');
	}
        for(int i = h-42; i < h; ++i){
		cursor_goto_pos(i+1,88);
		printf("%c",'@');
	} */
		
}

/*void plot_riemann(double* intermediate_variables){
	double l_state,1_state,2_state,r_state;
	switch(g_state[0]){
		case 'S':
			
			break;
		case 'R':

			break;
		default:
			printf("ERROR: Invaild State!!!");
	}
	switch(g_states[2]){
		case 'S':

			break;

		case 'R':
			
			break;
		default:
			printf("ERROR: Invaild State!!!");
	}	
}*/

