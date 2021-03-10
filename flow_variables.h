#define g_cL sqrt(g_gamma*g_PL/g_pL) // NOTE: Only works for an adiabatic gas FIX LATER!!!!!!!!!
#define g_cR sqrt(g_gamma * g_PR/g_pR) // NOTE: Only works for an adiabatic gas FIX LATER!!!!!!!!
#define g_beta (g_gamma + 1)/(g_gamma - 1)
#define g_alpha (g_gamma - 1)/(2*g_gamma)

static const double g_gamma = 5./3.;

static const double g_pL = 3.0;
static const double g_PL = 3.0;
static const double g_vL = 0.0;

static const double g_pR = 1.0;
static const double g_PR = 1.0;
static const double g_vR = 0.0;

char* g_states = NULL;
