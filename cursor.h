#include <stdio.h>
#include <sys/ioctl.h>
#include <unistd.h>


int get_height(void){
	struct winsize w;
	ioctl(STDOUT_FILENO,TIOCGWINSZ,&w);
	return w.ws_row;
}

int get_width(void){
        struct winsize w;
	ioctl(STDOUT_FILENO,TIOCGWINSZ,&w);
	return w.ws_col;
}

void cursor_go_home_clear(void){
	printf("\033[2J\033[H");	
}

void cursor_clear_screen(void){
	printf("\033[2J");
}

void cursor_goto_home(void){
	printf("\033[H");
}

void cursor_move_up(int x){
	printf("\033[xA");
}

void cursor_move_down(int x){
	 printf("\033[xB");
}

void cursor_move_right(int x){
	 printf("\033[xC");
}

void cursor_move_left(int x){
	 printf("\033[xD");
}

void cursor_goto_pos(int x,int y){
	printf("%c[%d;%df", 0x1B, x, y);
	fflush(stdout);
}

