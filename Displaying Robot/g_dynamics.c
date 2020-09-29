/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) g_dynamics_ ## ID
#endif

#include <math.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_f1 CASADI_PREFIX(f1)
#define casadi_f2 CASADI_PREFIX(f2)
#define casadi_f3 CASADI_PREFIX(f3)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

static const casadi_int casadi_s0[7] = {3, 1, 0, 3, 0, 1, 2};
static const casadi_int casadi_s1[15] = {3, 3, 0, 3, 6, 9, 0, 1, 2, 0, 1, 2, 0, 1, 2};
static const casadi_int casadi_s2[11] = {3, 3, 0, 1, 3, 5, 0, 0, 1, 1, 2};

/* g_dynamics_JSIM:(q[3])->(H[3x3]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a3, a4, a5, a6, a7, a8, a9;
  a0=20.;
  a1=2.9999999999999999e-01;
  a2=arg[0]? arg[0][0] : 0;
  a3=cos(a2);
  a3=(a1*a3);
  a4=(a0*a3);
  a4=(a4*a3);
  a3=sin(a2);
  a3=(a1*a3);
  a5=(a0*a3);
  a5=(a5*a3);
  a4=(a4+a5);
  a5=1.;
  a4=(a4+a5);
  a3=5.9999999999999998e-01;
  a6=cos(a2);
  a6=(a3*a6);
  a7=(a0*a6);
  a8=(a7*a6);
  a9=sin(a2);
  a9=(a3*a9);
  a10=(a0*a9);
  a11=(a10*a9);
  a8=(a8+a11);
  a4=(a4+a8);
  a8=cos(a2);
  a8=(a3*a8);
  a11=(a0*a8);
  a12=(a11*a8);
  a2=sin(a2);
  a2=(a3*a2);
  a13=(a0*a2);
  a14=(a13*a2);
  a12=(a12+a14);
  a4=(a4+a12);
  if (res[0]!=0) res[0][0]=a4;
  a4=arg[0]? arg[0][1] : 0;
  a12=cos(a4);
  a12=(a1*a12);
  a14=(a0*a12);
  a6=(a14*a6);
  a15=sin(a4);
  a15=(a1*a15);
  a16=(a0*a15);
  a9=(a16*a9);
  a6=(a6+a9);
  a9=cos(a4);
  a9=(a3*a9);
  a17=(a0*a9);
  a18=(a17*a8);
  a4=sin(a4);
  a3=(a3*a4);
  a4=(a0*a3);
  a19=(a4*a2);
  a18=(a18+a19);
  a6=(a6+a18);
  if (res[0]!=0) res[0][1]=a6;
  a6=arg[0]? arg[0][2] : 0;
  a18=cos(a6);
  a18=(a1*a18);
  a19=(a0*a18);
  a8=(a19*a8);
  a6=sin(a6);
  a1=(a1*a6);
  a0=(a0*a1);
  a2=(a0*a2);
  a8=(a8+a2);
  if (res[0]!=0) res[0][2]=a8;
  a7=(a7*a12);
  a10=(a10*a15);
  a7=(a7+a10);
  a10=(a11*a9);
  a8=(a13*a3);
  a10=(a10+a8);
  a7=(a7+a10);
  if (res[0]!=0) res[0][3]=a7;
  a14=(a14*a12);
  a16=(a16*a15);
  a14=(a14+a16);
  a14=(a14+a5);
  a16=(a17*a9);
  a15=(a4*a3);
  a16=(a16+a15);
  a14=(a14+a16);
  if (res[0]!=0) res[0][4]=a14;
  a9=(a19*a9);
  a3=(a0*a3);
  a9=(a9+a3);
  if (res[0]!=0) res[0][5]=a9;
  a11=(a11*a18);
  a13=(a13*a1);
  a11=(a11+a13);
  if (res[0]!=0) res[0][6]=a11;
  a17=(a17*a18);
  a4=(a4*a1);
  a17=(a17+a4);
  if (res[0]!=0) res[0][7]=a17;
  a19=(a19*a18);
  a0=(a0*a1);
  a19=(a19+a0);
  a19=(a19+a5);
  if (res[0]!=0) res[0][8]=a19;
  return 0;
}

CASADI_SYMBOL_EXPORT int g_dynamics_JSIM(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int g_dynamics_JSIM_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int g_dynamics_JSIM_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_dynamics_JSIM_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int g_dynamics_JSIM_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_dynamics_JSIM_release(int mem) {
}

CASADI_SYMBOL_EXPORT void g_dynamics_JSIM_incref(void) {
}

CASADI_SYMBOL_EXPORT void g_dynamics_JSIM_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int g_dynamics_JSIM_n_in(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_int g_dynamics_JSIM_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real g_dynamics_JSIM_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_dynamics_JSIM_name_in(casadi_int i){
  switch (i) {
    case 0: return "q";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_dynamics_JSIM_name_out(casadi_int i){
  switch (i) {
    case 0: return "H";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_dynamics_JSIM_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_dynamics_JSIM_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int g_dynamics_JSIM_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 1;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

/* g_dynamics_RHS:(q[3],v[3],u[3])->(RHS[3]) */
static int casadi_f1(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a4, a5, a6, a7, a8, a9;
  a0=arg[2]? arg[2][0] : 0;
  a1=arg[2]? arg[2][1] : 0;
  a0=(a0-a1);
  a2=5.0000000000000000e-01;
  a3=2.9999999999999999e-01;
  a4=arg[0]? arg[0][0] : 0;
  a5=sin(a4);
  a5=(a3*a5);
  a6=20.;
  a7=cos(a4);
  a7=(a3*a7);
  a8=(a6*a7);
  a8=(a5*a8);
  a9=(a6*a5);
  a9=(a9*a7);
  a8=(a8+a9);
  a9=cos(a4);
  a9=(a3*a9);
  a7=sin(a4);
  a7=(a3*a7);
  a10=(a6*a7);
  a10=(a9*a10);
  a9=(a6*a9);
  a9=(a9*a7);
  a10=(a10+a9);
  a8=(a8-a10);
  a10=5.9999999999999998e-01;
  a9=sin(a4);
  a9=(a10*a9);
  a7=cos(a4);
  a7=(a10*a7);
  a11=(a6*a7);
  a12=(a9*a11);
  a13=(a6*a9);
  a14=(a13*a7);
  a12=(a12+a14);
  a14=cos(a4);
  a14=(a10*a14);
  a15=sin(a4);
  a15=(a10*a15);
  a16=(a6*a15);
  a17=(a14*a16);
  a18=(a6*a14);
  a19=(a18*a15);
  a17=(a17+a19);
  a12=(a12-a17);
  a8=(a8+a12);
  a12=sin(a4);
  a12=(a10*a12);
  a17=cos(a4);
  a17=(a10*a17);
  a19=(a6*a17);
  a20=(a12*a19);
  a21=(a6*a12);
  a22=(a21*a17);
  a20=(a20+a22);
  a22=cos(a4);
  a22=(a10*a22);
  a4=sin(a4);
  a4=(a10*a4);
  a23=(a6*a4);
  a24=(a22*a23);
  a25=(a6*a22);
  a26=(a25*a4);
  a24=(a24+a26);
  a20=(a20-a24);
  a8=(a8+a20);
  a20=arg[1]? arg[1][0] : 0;
  a8=(a8*a20);
  a8=(a2*a8);
  a8=(a8*a20);
  a24=arg[0]? arg[0][1] : 0;
  a26=sin(a24);
  a26=(a3*a26);
  a11=(a26*a11);
  a27=cos(a24);
  a27=(a3*a27);
  a16=(a27*a16);
  a11=(a11-a16);
  a16=sin(a24);
  a16=(a10*a16);
  a28=(a16*a19);
  a29=cos(a24);
  a29=(a10*a29);
  a30=(a29*a23);
  a28=(a28-a30);
  a11=(a11+a28);
  a11=(a11*a20);
  a28=cos(a24);
  a28=(a3*a28);
  a13=(a13*a28);
  a30=sin(a24);
  a30=(a3*a30);
  a18=(a18*a30);
  a13=(a13-a18);
  a18=cos(a24);
  a18=(a10*a18);
  a31=(a21*a18);
  a24=sin(a24);
  a10=(a10*a24);
  a24=(a25*a10);
  a31=(a31-a24);
  a13=(a13+a31);
  a31=arg[1]? arg[1][1] : 0;
  a13=(a13*a31);
  a11=(a11+a13);
  a11=(a2*a11);
  a11=(a11*a31);
  a8=(a8+a11);
  a11=arg[0]? arg[0][2] : 0;
  a13=sin(a11);
  a13=(a3*a13);
  a19=(a13*a19);
  a24=cos(a11);
  a24=(a3*a24);
  a23=(a24*a23);
  a19=(a19-a23);
  a19=(a19*a20);
  a23=cos(a11);
  a23=(a3*a23);
  a21=(a21*a23);
  a11=sin(a11);
  a3=(a3*a11);
  a25=(a25*a3);
  a21=(a21-a25);
  a25=arg[1]? arg[1][2] : 0;
  a21=(a21*a25);
  a19=(a19+a21);
  a19=(a2*a19);
  a19=(a19*a25);
  a8=(a8+a19);
  a19=-9.8000000000000007e+00;
  a5=(a6*a5);
  a5=(a19*a5);
  a21=(a6*a9);
  a21=(a19*a21);
  a5=(a5+a21);
  a21=(a6*a12);
  a21=(a19*a21);
  a5=(a5+a21);
  a8=(a8+a5);
  a8=(a8+a20);
  a0=(a0-a8);
  if (res[0]!=0) res[0][0]=a0;
  a0=arg[2]? arg[2][2] : 0;
  a1=(a1-a0);
  a8=(a6*a26);
  a7=(a8*a7);
  a5=(a6*a27);
  a15=(a5*a15);
  a7=(a7-a15);
  a15=(a6*a16);
  a21=(a15*a17);
  a11=(a6*a29);
  a32=(a11*a4);
  a21=(a21-a32);
  a7=(a7+a21);
  a7=(a7*a20);
  a21=(a6*a28);
  a9=(a9*a21);
  a32=(a6*a30);
  a14=(a14*a32);
  a9=(a9-a14);
  a14=(a6*a18);
  a33=(a12*a14);
  a34=(a6*a10);
  a35=(a22*a34);
  a33=(a33-a35);
  a9=(a9+a33);
  a9=(a9*a31);
  a7=(a7+a9);
  a7=(a2*a7);
  a7=(a7*a20);
  a21=(a26*a21);
  a8=(a8*a28);
  a21=(a21+a8);
  a27=(a27*a32);
  a5=(a5*a30);
  a27=(a27+a5);
  a21=(a21-a27);
  a27=(a16*a14);
  a5=(a15*a18);
  a27=(a27+a5);
  a5=(a29*a34);
  a30=(a11*a10);
  a5=(a5+a30);
  a27=(a27-a5);
  a21=(a21+a27);
  a21=(a21*a31);
  a21=(a2*a21);
  a21=(a21*a31);
  a7=(a7+a21);
  a14=(a13*a14);
  a34=(a24*a34);
  a14=(a14-a34);
  a14=(a14*a31);
  a15=(a15*a23);
  a11=(a11*a3);
  a15=(a15-a11);
  a15=(a15*a25);
  a14=(a14+a15);
  a14=(a2*a14);
  a14=(a14*a25);
  a7=(a7+a14);
  a26=(a6*a26);
  a26=(a19*a26);
  a14=(a6*a16);
  a14=(a19*a14);
  a26=(a26+a14);
  a7=(a7+a26);
  a7=(a7+a31);
  a1=(a1-a7);
  if (res[0]!=0) res[0][1]=a1;
  a1=(a6*a13);
  a17=(a1*a17);
  a7=(a6*a24);
  a4=(a7*a4);
  a17=(a17-a4);
  a17=(a17*a20);
  a4=(a6*a23);
  a12=(a12*a4);
  a26=(a6*a3);
  a22=(a22*a26);
  a12=(a12-a22);
  a12=(a12*a25);
  a17=(a17+a12);
  a17=(a2*a17);
  a17=(a17*a20);
  a18=(a1*a18);
  a10=(a7*a10);
  a18=(a18-a10);
  a18=(a18*a31);
  a16=(a16*a4);
  a29=(a29*a26);
  a16=(a16-a29);
  a16=(a16*a25);
  a18=(a18+a16);
  a18=(a2*a18);
  a18=(a18*a31);
  a17=(a17+a18);
  a4=(a13*a4);
  a1=(a1*a23);
  a4=(a4+a1);
  a24=(a24*a26);
  a7=(a7*a3);
  a24=(a24+a7);
  a4=(a4-a24);
  a4=(a4*a25);
  a2=(a2*a4);
  a2=(a2*a25);
  a17=(a17+a2);
  a6=(a6*a13);
  a19=(a19*a6);
  a17=(a17+a19);
  a17=(a17+a25);
  a0=(a0-a17);
  if (res[0]!=0) res[0][2]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int g_dynamics_RHS(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f1(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int g_dynamics_RHS_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int g_dynamics_RHS_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_dynamics_RHS_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int g_dynamics_RHS_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_dynamics_RHS_release(int mem) {
}

CASADI_SYMBOL_EXPORT void g_dynamics_RHS_incref(void) {
}

CASADI_SYMBOL_EXPORT void g_dynamics_RHS_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int g_dynamics_RHS_n_in(void) { return 3;}

CASADI_SYMBOL_EXPORT casadi_int g_dynamics_RHS_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real g_dynamics_RHS_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_dynamics_RHS_name_in(casadi_int i){
  switch (i) {
    case 0: return "q";
    case 1: return "v";
    case 2: return "u";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_dynamics_RHS_name_out(casadi_int i){
  switch (i) {
    case 0: return "RHS";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_dynamics_RHS_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_dynamics_RHS_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int g_dynamics_RHS_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 3;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

/* g_dynamics_ControlMap:(q[3])->(B[3x3,5nz]) */
static int casadi_f2(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1;
  a0=1.;
  if (res[0]!=0) res[0][0]=a0;
  a1=-1.;
  if (res[0]!=0) res[0][1]=a1;
  if (res[0]!=0) res[0][2]=a0;
  if (res[0]!=0) res[0][3]=a1;
  if (res[0]!=0) res[0][4]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int g_dynamics_ControlMap(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f2(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int g_dynamics_ControlMap_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int g_dynamics_ControlMap_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_dynamics_ControlMap_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int g_dynamics_ControlMap_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_dynamics_ControlMap_release(int mem) {
}

CASADI_SYMBOL_EXPORT void g_dynamics_ControlMap_incref(void) {
}

CASADI_SYMBOL_EXPORT void g_dynamics_ControlMap_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int g_dynamics_ControlMap_n_in(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_int g_dynamics_ControlMap_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real g_dynamics_ControlMap_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_dynamics_ControlMap_name_in(casadi_int i){
  switch (i) {
    case 0: return "q";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_dynamics_ControlMap_name_out(casadi_int i){
  switch (i) {
    case 0: return "B";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_dynamics_ControlMap_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_dynamics_ControlMap_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int g_dynamics_ControlMap_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 1;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

/* g_control_ForcesForComputedTorqueController:(q[3],v[3])->(c[3]) */
static int casadi_f3(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a4, a5, a6, a7, a8, a9;
  a0=5.0000000000000000e-01;
  a1=2.9999999999999999e-01;
  a2=arg[0]? arg[0][0] : 0;
  a3=sin(a2);
  a3=(a1*a3);
  a4=20.;
  a5=cos(a2);
  a5=(a1*a5);
  a6=(a4*a5);
  a6=(a3*a6);
  a7=(a4*a3);
  a7=(a7*a5);
  a6=(a6+a7);
  a7=cos(a2);
  a7=(a1*a7);
  a5=sin(a2);
  a5=(a1*a5);
  a8=(a4*a5);
  a8=(a7*a8);
  a7=(a4*a7);
  a7=(a7*a5);
  a8=(a8+a7);
  a6=(a6-a8);
  a8=5.9999999999999998e-01;
  a7=sin(a2);
  a7=(a8*a7);
  a5=cos(a2);
  a5=(a8*a5);
  a9=(a4*a5);
  a10=(a7*a9);
  a11=(a4*a7);
  a12=(a11*a5);
  a10=(a10+a12);
  a12=cos(a2);
  a12=(a8*a12);
  a13=sin(a2);
  a13=(a8*a13);
  a14=(a4*a13);
  a15=(a12*a14);
  a16=(a4*a12);
  a17=(a16*a13);
  a15=(a15+a17);
  a10=(a10-a15);
  a6=(a6+a10);
  a10=sin(a2);
  a10=(a8*a10);
  a15=cos(a2);
  a15=(a8*a15);
  a17=(a4*a15);
  a18=(a10*a17);
  a19=(a4*a10);
  a20=(a19*a15);
  a18=(a18+a20);
  a20=cos(a2);
  a20=(a8*a20);
  a2=sin(a2);
  a2=(a8*a2);
  a21=(a4*a2);
  a22=(a20*a21);
  a23=(a4*a20);
  a24=(a23*a2);
  a22=(a22+a24);
  a18=(a18-a22);
  a6=(a6+a18);
  a18=arg[1]? arg[1][0] : 0;
  a6=(a6*a18);
  a6=(a0*a6);
  a6=(a6*a18);
  a22=arg[0]? arg[0][1] : 0;
  a24=sin(a22);
  a24=(a1*a24);
  a9=(a24*a9);
  a25=cos(a22);
  a25=(a1*a25);
  a14=(a25*a14);
  a9=(a9-a14);
  a14=sin(a22);
  a14=(a8*a14);
  a26=(a14*a17);
  a27=cos(a22);
  a27=(a8*a27);
  a28=(a27*a21);
  a26=(a26-a28);
  a9=(a9+a26);
  a9=(a9*a18);
  a26=cos(a22);
  a26=(a1*a26);
  a11=(a11*a26);
  a28=sin(a22);
  a28=(a1*a28);
  a16=(a16*a28);
  a11=(a11-a16);
  a16=cos(a22);
  a16=(a8*a16);
  a29=(a19*a16);
  a22=sin(a22);
  a8=(a8*a22);
  a22=(a23*a8);
  a29=(a29-a22);
  a11=(a11+a29);
  a29=arg[1]? arg[1][1] : 0;
  a11=(a11*a29);
  a9=(a9+a11);
  a9=(a0*a9);
  a9=(a9*a29);
  a6=(a6+a9);
  a9=arg[0]? arg[0][2] : 0;
  a11=sin(a9);
  a11=(a1*a11);
  a17=(a11*a17);
  a22=cos(a9);
  a22=(a1*a22);
  a21=(a22*a21);
  a17=(a17-a21);
  a17=(a17*a18);
  a21=cos(a9);
  a21=(a1*a21);
  a19=(a19*a21);
  a9=sin(a9);
  a1=(a1*a9);
  a23=(a23*a1);
  a19=(a19-a23);
  a23=arg[1]? arg[1][2] : 0;
  a19=(a19*a23);
  a17=(a17+a19);
  a17=(a0*a17);
  a17=(a17*a23);
  a6=(a6+a17);
  a17=-9.8000000000000007e+00;
  a3=(a4*a3);
  a3=(a17*a3);
  a19=(a4*a7);
  a19=(a17*a19);
  a3=(a3+a19);
  a19=(a4*a10);
  a19=(a17*a19);
  a3=(a3+a19);
  a6=(a6+a3);
  a6=(a6+a18);
  if (res[0]!=0) res[0][0]=a6;
  a6=(a4*a24);
  a5=(a6*a5);
  a3=(a4*a25);
  a13=(a3*a13);
  a5=(a5-a13);
  a13=(a4*a14);
  a19=(a13*a15);
  a9=(a4*a27);
  a30=(a9*a2);
  a19=(a19-a30);
  a5=(a5+a19);
  a5=(a5*a18);
  a19=(a4*a26);
  a7=(a7*a19);
  a30=(a4*a28);
  a12=(a12*a30);
  a7=(a7-a12);
  a12=(a4*a16);
  a31=(a10*a12);
  a32=(a4*a8);
  a33=(a20*a32);
  a31=(a31-a33);
  a7=(a7+a31);
  a7=(a7*a29);
  a5=(a5+a7);
  a5=(a0*a5);
  a5=(a5*a18);
  a19=(a24*a19);
  a6=(a6*a26);
  a19=(a19+a6);
  a25=(a25*a30);
  a3=(a3*a28);
  a25=(a25+a3);
  a19=(a19-a25);
  a25=(a14*a12);
  a3=(a13*a16);
  a25=(a25+a3);
  a3=(a27*a32);
  a28=(a9*a8);
  a3=(a3+a28);
  a25=(a25-a3);
  a19=(a19+a25);
  a19=(a19*a29);
  a19=(a0*a19);
  a19=(a19*a29);
  a5=(a5+a19);
  a12=(a11*a12);
  a32=(a22*a32);
  a12=(a12-a32);
  a12=(a12*a29);
  a13=(a13*a21);
  a9=(a9*a1);
  a13=(a13-a9);
  a13=(a13*a23);
  a12=(a12+a13);
  a12=(a0*a12);
  a12=(a12*a23);
  a5=(a5+a12);
  a24=(a4*a24);
  a24=(a17*a24);
  a12=(a4*a14);
  a12=(a17*a12);
  a24=(a24+a12);
  a5=(a5+a24);
  a5=(a5+a29);
  if (res[0]!=0) res[0][1]=a5;
  a5=(a4*a11);
  a15=(a5*a15);
  a24=(a4*a22);
  a2=(a24*a2);
  a15=(a15-a2);
  a15=(a15*a18);
  a2=(a4*a21);
  a10=(a10*a2);
  a12=(a4*a1);
  a20=(a20*a12);
  a10=(a10-a20);
  a10=(a10*a23);
  a15=(a15+a10);
  a15=(a0*a15);
  a15=(a15*a18);
  a16=(a5*a16);
  a8=(a24*a8);
  a16=(a16-a8);
  a16=(a16*a29);
  a14=(a14*a2);
  a27=(a27*a12);
  a14=(a14-a27);
  a14=(a14*a23);
  a16=(a16+a14);
  a16=(a0*a16);
  a16=(a16*a29);
  a15=(a15+a16);
  a2=(a11*a2);
  a5=(a5*a21);
  a2=(a2+a5);
  a22=(a22*a12);
  a24=(a24*a1);
  a22=(a22+a24);
  a2=(a2-a22);
  a2=(a2*a23);
  a0=(a0*a2);
  a0=(a0*a23);
  a15=(a15+a0);
  a4=(a4*a11);
  a17=(a17*a4);
  a15=(a15+a17);
  a15=(a15+a23);
  if (res[0]!=0) res[0][2]=a15;
  return 0;
}

CASADI_SYMBOL_EXPORT int g_control_ForcesForComputedTorqueController(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f3(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int g_control_ForcesForComputedTorqueController_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int g_control_ForcesForComputedTorqueController_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_control_ForcesForComputedTorqueController_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int g_control_ForcesForComputedTorqueController_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_control_ForcesForComputedTorqueController_release(int mem) {
}

CASADI_SYMBOL_EXPORT void g_control_ForcesForComputedTorqueController_incref(void) {
}

CASADI_SYMBOL_EXPORT void g_control_ForcesForComputedTorqueController_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int g_control_ForcesForComputedTorqueController_n_in(void) { return 2;}

CASADI_SYMBOL_EXPORT casadi_int g_control_ForcesForComputedTorqueController_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real g_control_ForcesForComputedTorqueController_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_control_ForcesForComputedTorqueController_name_in(casadi_int i){
  switch (i) {
    case 0: return "q";
    case 1: return "v";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_control_ForcesForComputedTorqueController_name_out(casadi_int i){
  switch (i) {
    case 0: return "c";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_control_ForcesForComputedTorqueController_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_control_ForcesForComputedTorqueController_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int g_control_ForcesForComputedTorqueController_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 2;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
