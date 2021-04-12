#define g_cL (double)(sqrt(g_gamma*g_PL/g_pL)) 
#define g_cR (double)(sqrt(g_gamma * g_PR/g_pR)) 
#define g_beta (double)((g_gamma + 1)/(g_gamma - 1))
#define g_alpha (double)((g_gamma - 1)/(2*g_gamma))

static const double g_gamma = 5./3.;

// Left States
double g_pL;
double g_PL;
double g_vL;

// Right States
double g_pR;
double g_PR;
double g_vR;

char g_states[2] = {' ',' '};
