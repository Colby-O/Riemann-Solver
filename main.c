#include "Riemann.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "plot.h"
//#include "cursor.h"
float func(float x){
	return sin(x); 
}

int main(void){	
	double* test = riemannSolver();
	/*
        printf("State 1: %c\n", g_states[0]);
	printf("State 2: %c\n", g_states[1]);
	printf("State 3: %c\n", g_states[2]);
        for(int i = 0; i < 8; i++){
		if(test[i] != MAX_FLOAT){
        		printf("%f ", test[i]);
		}
		else{
			printf("%s ", "N/A");
		}
	}
	printf("\n");
	free(test);*/	
	init_grids(test);
	//printf("Height: %d\n",get_height());
	//printf("Width: %d\n",get_width());
}
