#ifndef config_h
#define config_h
#include <stdbool.h> // for bool

// define class name and unique id
#define MODEL_IDENTIFIER MassSpringOscillator
#define INSTANTIATION_TOKEN "{9B6E0A3D-27C5-4F1B-8E7A-1D4C6F08B932}"

#define MODEL_EXCHANGE

#define MAX_CONTINUOUS_STATES 2

#define SET_FLOAT64
#define GET_BOOLEAN
#define SET_BOOLEAN

#define FIXED_SOLVER_STEP 1e-2
#define DEFAULT_STOP_TIME 20

typedef enum {
    vr_time, vr_x, vr_der_x, vr_v, vr_der_v, vr_F,
    vr_m, vr_k, vr_E, vr_res1, vr_ode_dae
} ValueReference;

typedef struct {
    double x;
    double der_x;
    double v;
    double der_v;
    double F;
    double m;
    double k;
    double E;
    double res1;
    bool ode_dae;
} ModelData;

#endif /* config_h */
