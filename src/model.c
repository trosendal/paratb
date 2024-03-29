#include <R_ext/Rdynload.h>
#include "SimInf.h"

/**
* Make sure the necessary macros are defined so that the compiler can
* replace them when compiling the model.  'SIMINF_MODEL_RUN' defines
* the function name of the function that will be called from R to run
* a trajectory of the model.  'SIMINF_R_INIT' is the name of the
* function that R will call when this model is loaded into
* R. 'SIMINF_FORCE_SYMBOLS' defines whether R allows the entry point
* for the run function to be searched for as a character string.  If
* this file is compiled from SimInf (when calling run), the macros are
* defined by SimInf before calling 'R CMD SHLIB'.  If this file is
* compiled as part of a package, then the definitions are set in the
* variable 'PKG_CPPFLAGS' in 'src/Makevars' and 'src/Makevars.in'.
*/
#if !defined(SIMINF_MODEL_RUN)
#  error Definition for 'SIMINF_MODEL_RUN' is missing.
#endif
#if !defined(SIMINF_R_INIT)
#  error Definition for 'SIMINF_R_INIT' is missing.
#endif
#if !defined(SIMINF_FORCE_SYMBOLS)
#  error Definition for 'SIMINF_FORCE_SYMBOLS' is missing.
#endif
#define SIMINF_STR(name) #name
#define SIMINF_CALLDEF(name, n) {SIMINF_STR(name), (DL_FUNC) &name, n}

/* Compartments */
enum{S0, S1, T1H, T1L, S2, T2H, T2L, E2H, E2L, S3, E3H, E3L, L3H, L3L, H3, BM, BMp, BMn};

enum{P1, GAMMA_E, GAMMA_L, GAMMA_H, BETA, BETA_A, PHI, ETA, SIGMA_H, SIGMA_L, NU, MU_B, RHO_1, MU_1,
     RHO_2, MU_2, MU_3, MU_4, INTERCEPT, COEF, DELTA2, SP2, IND_SE};

/* Calculate the total number of adults in the population */
static int nadult(const int *u)
{
     return u[S3]+u[E3H]+u[E3L]+u[L3H]+u[L3L]+u[H3];
}

/* Calculate the population size */
static int npop(const int *u)
{
     return u[S1]+u[T1H]+u[T1L]+u[S2]+u[T2H]+u[T2L]+u[E2H]+u[E2L]+u[S3]+u[E3H]+u[E3L]+u[L3H]+u[L3L]+u[H3];
}

/* Calculate the number infectious in the population */
static int ipop(const int *u)
{
     return u[T1H]+u[T1L]+u[T2H]+u[T2L]+u[L3L]+u[L3H]+u[H3];
}

/* Calculate the number of infectious adults */
static int infectious_adults(const int *u)
{
     return u[L3L]+u[L3H]+u[H3];
}

/* Calculate the number of infected adults */
static int infected_adults(const int *u)
{
    return infectious_adults(u) + u[E3L]+u[E3H];
}

static double fraction_infected_adults(const int *u)
{
    const double I3 = infected_adults(u);
    const double N3 = nadult(u);
    return I3 / N3;
}

/* The calculation of the bulk milk test sensitivity dependant on the
 * number of animals positive contributing to the milk. */
static double bulk_milk_sensitivity(const int *u, const double* gdata)
{
    const double fi = fraction_infected_adults(u) * gdata[IND_SE];
    const double k = gdata[COEF];
    const double m = gdata[INTERCEPT];
    const double y = k * -log2(fi) + m;
    const double ey = exp(y);

    return ey / (1.0 + ey);
}

/* Rate of animals born into the S0 dummy compartment */
static double trFun1(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     const int N3 = nadult(u);
     return gdata[MU_B]*N3;
}

/* Rate that animals are 'born' from S0 to S1 */
static double trFun2(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     const int N3 = nadult(u);

     if (N3 < 1)
          return u[S0]*gdata[P1];

     return u[S0]*gdata[P1]*((u[S3]+(1.0-gdata[GAMMA_E])*(u[E3H]+u[E3L])+(1.0-gdata[GAMMA_L])*(u[L3H]+u[L3L])+(1-gdata[GAMMA_H])*u[H3])/((double)N3));
}

/* Rate that animals are 'born' from S0 to T1H */
static double trFun3(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     const int N3 = nadult(u);

     if (N3 < 1)
          return 0.0;

     return u[S0]*gdata[P1]*((gdata[GAMMA_E]*(u[E3H]+u[E3L])+gdata[GAMMA_L]*(u[L3H]+u[L3L])+gdata[GAMMA_H]*u[H3])/((double)N3))*gdata[ETA];
}

/* Rate that animals are 'born' from S0 to T1L */
static double trFun4(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     const int N3 = nadult(u);
     if (N3 < 1)
          return 0.0;
     return u[S0]*gdata[P1]*((gdata[GAMMA_E]*(u[E3H]+u[E3L])+gdata[GAMMA_L]*(u[L3H]+u[L3L])+gdata[GAMMA_H]*u[H3])/N3)*(1-gdata[ETA]);
}

/* Horizontal transmission from S1 to T1H */
static double trFun5(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     const int N = npop(u);
     const int I = ipop(u);
     if (N < 1)
          return 0.0;
     return u[S1]*(I)*(gdata[BETA]/((double)N))*gdata[ETA];
}

/* Horizontal transmission from S1 to T1L */
static double trFun6(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     const int N = npop(u);
     const int I = ipop(u);
     if (N < 1)
          return 0.0;
     return u[S1]*I*(gdata[BETA]/N)*(1-gdata[ETA]);
}

/* Ageing of calves from S1 to S2 */
static double trFun7(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[RHO_1]*u[S1];
}

/* Ageing of calves from T1H to T2H */
static double trFun8(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[RHO_1]*u[T1H];
}

/* Ageing of calves from T1L to T2L */
static double trFun9(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[RHO_1]*u[T1L];
}

/* Death in S1 */
static double trFun10(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[S1]*gdata[MU_1];
}

/* Death in T1H */
static double trFun11(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[T1H]*gdata[MU_1];
}

/* Death in T1L */
static double trFun12(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[T1L]*gdata[MU_1];
}

/* Horizontal transmission in heifers to T2H */
static double trFun13(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     const int N = npop(u);
     if (N < 1)
          return 0.0;
     return u[S2]*(u[T2H]+u[T2L])*(gdata[BETA_A]/N)*gdata[ETA];
}

/* Horizontal transmission in heifers to T2L */
static double trFun14(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     const int N = npop(u);
     if (N < 1)
          return 0.0;
     return u[S2]*(u[T2H]+u[T2L])*(gdata[BETA_A]/N)*(1-gdata[ETA]);
}

/* Heifers going from T2H to E2H */
static double trFun15(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[T2H]*gdata[PHI];
}
/* Heifers going from T2L to E2L */
static double trFun16(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[T2L]*gdata[PHI];
}

/* Ageing from from S2 to S3 */
static double trFun17(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[RHO_2]*u[S2];
}

/* Ageing from from E2H to E3H */
static double trFun18(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[E2H]*gdata[RHO_2];
}

/* Ageing from from E2L to E3L */
static double trFun19(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[E2L]*gdata[RHO_2];
}

/* Ageing from from T2H to T3H */
static double trFun20(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[T2H]*gdata[RHO_2];
}

/* Ageing from from T2L to T3L */
static double trFun21(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[T2L]*gdata[RHO_2];
}

/* Death in S2 */
static double trFun22(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[S2]*gdata[MU_2];
}

/* Death in T2H */
static double trFun23(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[T2H]*gdata[MU_2];
}

/* Death in T2L */
static double trFun24(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[T2L]*gdata[MU_2];
}

/* Death in E2H */
static double trFun25(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[E2H]*gdata[MU_2];
}

/* Death in E2L */
static double trFun26(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return u[E2L]*gdata[MU_2];
}

/* Horizontal transmission in S3 to E3L */
static double trFun27(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     const int N = npop(u);
     if (N < 1)
          return 0.0;
     return u[S3]*(u[L3H]+u[L3L]+u[H3])*(gdata[BETA_A]/N);
}

/* Rate at which adults move from latent to shedding - E3H --> L3H*/
static double trFun28(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[SIGMA_H]*u[E3H];
}

/* Rate at which adults move from latent to shedding - E3L --> L3L*/
static double trFun29(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[SIGMA_L]*u[E3L];
}

/* Rate which adults move from L3H to H3 */
static double trFun30(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[NU]*u[L3H];
}

/* Death in S3 */
static double trFun31(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[MU_3]*u[S3];
}

/* Death in E3H */
static double trFun32(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[MU_3]*u[E3H];
}

/* Death in E3L */
static double trFun33(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[MU_3]*u[E3L];
}

/* Death in L3H */
static double trFun34(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[MU_3]*u[L3H];
}

/* Death in L3L */
static double trFun35(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[MU_3]*u[L3L];
}

/* Death in H3 */
static double trFun36(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[MU_3]*u[H3];
}

/* Recovery in L3H. This is a model component that mimics active
 * selection against animals with clinical signs associated with para
 * TB (poor milk yield bad body condition etc.) Not just a death rate
 * because we wish to use real scheduled events which would be
 * disrupted by unplanned removal of animals*/
static double trFun37(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[MU_4]*u[L3H];
}

/* Recovery in L3L. Analogous justification to above*/
static double trFun38(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[MU_4]*u[L3L];
}

/* Recovery in H3. Analogous justification to above */
static double trFun39(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     return gdata[MU_4]*u[H3];
}

/* Rate at which bulk milk testing is scheduled by placing
 * a unit in the dummy BM compartment*/
static double trFun43(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     if (nadult(u) < 1)
          return 0.0;
     return gdata[DELTA2];
}

/* Rate at which the scheduled bulk milk tests become positive populating the dummy BMp compartment*/
static double trFun44(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     if (nadult(u) < 1)
          return 0.0;
     if (infected_adults(u) < 1)
         return u[BM] * gdata[P1] * (1.0 - gdata[SP2]);
     return u[BM] * gdata[P1] * bulk_milk_sensitivity(u, gdata);
}

/* Rate at which the scheduled bulk milk tests become negative populating the dummy BMn compartment*/
static double trFun45(
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     double t)
{
     if (nadult(u) < 1)
          return 0.0;
     if (infected_adults(u) < 1)
         return u[BM] * gdata[P1] * gdata[SP2];
     return u[BM] * gdata[P1] * (1.0 - bulk_milk_sensitivity(u, gdata));
}

static int ptsFun(
     double *v_new,
     const int *u,
     const double *v,
     const double *ldata,
     const double *gdata,
     int node,
     double t)
{
     return 0;
}

static SEXP SIMINF_MODEL_RUN(SEXP model, SEXP solver)
{
    static SEXP(*SimInf_run)(SEXP, SEXP, TRFun*, PTSFun) = NULL;
    TRFun tr_fun[] = {&trFun1, &trFun2, &trFun3, &trFun4, &trFun5, &trFun6,
                       &trFun7, &trFun8, &trFun9, &trFun10, &trFun11, &trFun12,
                       &trFun13, &trFun14, &trFun15, &trFun16, &trFun17, &trFun18,
                       &trFun19, &trFun20, &trFun21, &trFun22, &trFun23, &trFun24,
                       &trFun25, &trFun26, &trFun27, &trFun28, &trFun29, &trFun30,
                       &trFun31, &trFun32, &trFun33, &trFun34, &trFun35, &trFun36,
                       &trFun37, &trFun38, &trFun39, &trFun43, &trFun44, &trFun45};

    if (!SimInf_run) {
        SimInf_run = (SEXP(*)(SEXP, SEXP, TRFun*, PTSFun))
            R_GetCCallable("SimInf", "SimInf_run");

        if (!SimInf_run) {
            Rf_error("Cannot find function 'SimInf_run'.");
        }
    }

    return SimInf_run(model, solver, tr_fun, &ptsFun);
}

static const R_CallMethodDef callMethods[] =
{
    SIMINF_CALLDEF(SIMINF_MODEL_RUN, 2),
    {NULL, NULL, 0}
};

void SIMINF_R_INIT(DllInfo *info)
{
     R_registerRoutines(info, NULL, callMethods, NULL, NULL);
     R_useDynamicSymbols(info, FALSE);
     R_forceSymbols(info, SIMINF_FORCE_SYMBOLS);
}
