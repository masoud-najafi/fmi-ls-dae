#include "config.h"
#include "model.h"


void setStartValues(ModelInstance *comp) {
    M(x)     = 1.0;
    M(der_x) = 0;
    M(v)     = 0;
    M(der_v) = 0;
    M(F)     = 0;
    M(m)     = 1.0;
    M(k)     = 1.0;
    M(E)     = 0;
    M(res1)  = 0;
}

Status calculateValues(ModelInstance *comp) {

    if (!M(ode_dae)) {
        // ODE mode: the residual is solved analytically for the spring force
        M(F) = -M(k) * M(x);
    }

    M(der_x) = M(v);
    M(der_v) = M(F) / M(m);

    // residual of the constitutive equation of the spring
    M(res1) = M(F) + M(k) * M(x);

    // invariant: total mechanical energy (kinetic + potential)
    M(E) = 0.5 * M(m) * M(v) * M(v) + 0.5 * M(k) * M(x) * M(x);

    return OK;
}

Status getBoolean(ModelInstance* comp, ValueReference vr, bool values[], size_t nValues, size_t* index) {

    ASSERT_NVALUES(1);

    switch (vr) {
        case vr_ode_dae: values[(*index)++] = M(ode_dae); return OK;
        default:
            logError(comp, "GetBoolean is not allowed for value reference %u.", vr);
            return Error;
    }
}

Status setBoolean(ModelInstance* comp, ValueReference vr, const bool values[], size_t nValues, size_t* index) {

    if (comp->state != ConfigurationMode && comp->state != ReconfigurationMode) {
        logError(comp, "Structural variables can only be set in Configuration Mode or Reconfiguration Mode.");
        return Error;
    }

    ASSERT_NVALUES(1);

    switch (vr) {
        case vr_ode_dae: M(ode_dae) = values[(*index)++]; return OK;
        default:
            logError(comp, "SetBoolean is not allowed for value reference %u.", vr);
            return Error;
    }
}

Status getFloat64(ModelInstance* comp, ValueReference vr, double values[], size_t nValues, size_t* index) {

    ASSERT_NVALUES(1);

    calculateValues(comp);

    switch (vr) {
        case vr_time:  values[(*index)++] = comp->time; return OK;
        case vr_x:     values[(*index)++] = M(x);       return OK;
        case vr_der_x: values[(*index)++] = M(der_x);   return OK;
        case vr_v:     values[(*index)++] = M(v);       return OK;
        case vr_der_v: values[(*index)++] = M(der_v);   return OK;
        case vr_F:     values[(*index)++] = M(F);       return OK;
        case vr_m:     values[(*index)++] = M(m);       return OK;
        case vr_k:     values[(*index)++] = M(k);       return OK;
        case vr_E:     values[(*index)++] = M(E);       return OK;
        case vr_res1:  values[(*index)++] = M(res1);    return OK;
        default:
            logError(comp, "Get Float64 is not allowed for value reference %u.", vr);
            return Error;
    }
}

Status setFloat64(ModelInstance* comp, ValueReference vr, const double values[], size_t nValues, size_t* index) {

    ASSERT_NVALUES(1);

    switch (vr) {
        case vr_x: M(x) = values[(*index)++]; return OK;
        case vr_v: M(v) = values[(*index)++]; return OK;
        case vr_F: M(F) = values[(*index)++]; return OK;
        case vr_m: M(m) = values[(*index)++]; return OK;
        case vr_k: M(k) = values[(*index)++]; return OK;
        default:
            logError(comp, "Set Float64 is not allowed for value reference %u.", vr);
            return Error;
    }
}

size_t getNumberOfContinuousStates(ModelInstance* comp) {
    UNUSED(comp);
    return 2;
}

Status getContinuousStates(ModelInstance *comp, double x[], size_t nx) {
    UNUSED(nx);
    x[0] = M(x);
    x[1] = M(v);
    return OK;
}

Status setContinuousStates(ModelInstance *comp, const double x[], size_t nx) {
    UNUSED(nx);
    M(x) = x[0];
    M(v) = x[1];
    calculateValues(comp);
    return OK;
}

Status getDerivatives(ModelInstance *comp, double dx[], size_t nx) {
    UNUSED(nx);
    calculateValues(comp);
    dx[0] = M(der_x);
    dx[1] = M(der_v);
    return OK;
}

Status eventUpdate(ModelInstance *comp) {

    comp->valuesOfContinuousStatesChanged   = false;
    comp->nominalsOfContinuousStatesChanged = false;
    comp->terminateSimulation               = false;
    comp->nextEventTimeDefined              = false;

    return OK;
}
