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
  #define CASADI_PREFIX(ID) g_InverseKinematics_ ## ID
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
static const casadi_int casadi_s1[8] = {4, 1, 0, 4, 0, 1, 2, 3};
static const casadi_int casadi_s2[15] = {4, 3, 0, 3, 6, 9, 1, 2, 3, 1, 2, 3, 1, 2, 3};

/* g_InverseKinematics_Task:(q[3])->(Task[4]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a2, a3, a4, a5, a6, a7, a8, a9;
  a0=0.;
  if (res[0]!=0) res[0][0]=a0;
  a0=10.;
  a1=5.0000000000000000e-01;
  a2=arg[0]? arg[0][0] : 0;
  a3=sin(a2);
  a4=(a1*a3);
  a5=cos(a2);
  a6=arg[0]? arg[0][1] : 0;
  a7=sin(a6);
  a8=(a5*a7);
  a9=cos(a6);
  a10=(a3*a9);
  a8=(a8+a10);
  a10=(a1*a8);
  a10=(a4+a10);
  a11=2.5000000000000000e-01;
  a12=cos(a6);
  a5=(a5*a12);
  a6=sin(a6);
  a13=(a3*a6);
  a5=(a5-a13);
  a13=arg[0]? arg[0][2] : 0;
  a14=sin(a13);
  a5=(a5*a14);
  a13=cos(a13);
  a15=(a8*a13);
  a5=(a5+a15);
  a15=(a11*a5);
  a15=(a10+a15);
  a15=(a0*a15);
  a8=(a11*a8);
  a4=(a4+a8);
  a4=(a0*a4);
  a3=(a11*a3);
  a3=(a0*a3);
  a4=(a4+a3);
  a15=(a15+a4);
  a4=30.;
  a15=(a15/a4);
  a15=(-a15);
  if (res[0]!=0) res[0][1]=a15;
  a15=cos(a2);
  a3=(a1*a15);
  a9=(a15*a9);
  a2=sin(a2);
  a7=(a2*a7);
  a9=(a9-a7);
  a7=(a1*a9);
  a7=(a3+a7);
  a13=(a9*a13);
  a2=(a2*a12);
  a6=(a15*a6);
  a2=(a2+a6);
  a2=(a2*a14);
  a13=(a13-a2);
  a13=(a11*a13);
  a7=(a7+a13);
  a7=(a0*a7);
  a9=(a11*a9);
  a3=(a3+a9);
  a3=(a0*a3);
  a11=(a11*a15);
  a0=(a0*a11);
  a3=(a3+a0);
  a7=(a7+a3);
  a7=(a7/a4);
  if (res[0]!=0) res[0][2]=a7;
  a1=(a1*a5);
  a10=(a10+a1);
  a10=(-a10);
  if (res[0]!=0) res[0][3]=a10;
  return 0;
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_Task(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_Task_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_Task_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_Task_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_Task_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_Task_release(int mem) {
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_Task_incref(void) {
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_Task_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int g_InverseKinematics_Task_n_in(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_int g_InverseKinematics_Task_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real g_InverseKinematics_Task_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_InverseKinematics_Task_name_in(casadi_int i){
  switch (i) {
    case 0: return "q";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_InverseKinematics_Task_name_out(casadi_int i){
  switch (i) {
    case 0: return "Task";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_InverseKinematics_Task_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_InverseKinematics_Task_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_Task_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 1;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

/* g_InverseKinematics_TaskJacobian:(q[3])->(TaskJacobian[4x3,9nz]) */
static int casadi_f1(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a3, a4, a5, a6, a7, a8, a9;
  a0=3.3333333333333333e-02;
  a1=10.;
  a2=5.0000000000000000e-01;
  a3=arg[0]? arg[0][0] : 0;
  a4=cos(a3);
  a5=(a2*a4);
  a6=arg[0]? arg[0][1] : 0;
  a7=cos(a6);
  a8=(a7*a4);
  a9=sin(a6);
  a10=sin(a3);
  a11=(a9*a10);
  a8=(a8-a11);
  a11=(a2*a8);
  a11=(a5+a11);
  a12=2.5000000000000000e-01;
  a13=arg[0]? arg[0][2] : 0;
  a14=cos(a13);
  a15=(a14*a8);
  a16=sin(a13);
  a17=cos(a6);
  a10=(a17*a10);
  a18=sin(a6);
  a19=(a18*a4);
  a10=(a10+a19);
  a10=(a16*a10);
  a15=(a15-a10);
  a10=(a12*a15);
  a10=(a11+a10);
  a10=(a1*a10);
  a8=(a12*a8);
  a5=(a5+a8);
  a5=(a1*a5);
  a4=(a12*a4);
  a4=(a1*a4);
  a5=(a5+a4);
  a10=(a10+a5);
  a10=(a0*a10);
  a10=(-a10);
  if (res[0]!=0) res[0][0]=a10;
  a10=sin(a3);
  a5=(a2*a10);
  a4=(a7*a10);
  a8=cos(a3);
  a19=(a9*a8);
  a4=(a4+a19);
  a19=(a2*a4);
  a19=(a5+a19);
  a20=(a14*a4);
  a8=(a17*a8);
  a21=(a18*a10);
  a8=(a8-a21);
  a8=(a16*a8);
  a20=(a20+a8);
  a20=(a12*a20);
  a19=(a19+a20);
  a19=(a1*a19);
  a4=(a12*a4);
  a5=(a5+a4);
  a5=(a1*a5);
  a10=(a12*a10);
  a10=(a1*a10);
  a5=(a5+a10);
  a19=(a19+a5);
  a19=(a0*a19);
  a19=(-a19);
  if (res[0]!=0) res[0][1]=a19;
  a15=(a2*a15);
  a11=(a11+a15);
  a11=(-a11);
  if (res[0]!=0) res[0][2]=a11;
  a11=cos(a3);
  a15=cos(a6);
  a19=(a11*a15);
  a5=sin(a3);
  a10=sin(a6);
  a4=(a5*a10);
  a19=(a19-a4);
  a4=(a2*a19);
  a20=(a14*a19);
  a8=sin(a6);
  a21=(a11*a8);
  a6=cos(a6);
  a22=(a5*a6);
  a21=(a21+a22);
  a21=(a16*a21);
  a20=(a20-a21);
  a21=(a12*a20);
  a21=(a4+a21);
  a21=(a1*a21);
  a19=(a12*a19);
  a19=(a1*a19);
  a21=(a21+a19);
  a21=(a0*a21);
  a21=(-a21);
  if (res[0]!=0) res[0][3]=a21;
  a21=cos(a3);
  a10=(a21*a10);
  a3=sin(a3);
  a15=(a3*a15);
  a10=(a10+a15);
  a15=(a2*a10);
  a14=(a14*a10);
  a6=(a21*a6);
  a8=(a3*a8);
  a6=(a6-a8);
  a16=(a16*a6);
  a14=(a14+a16);
  a14=(a12*a14);
  a15=(a15+a14);
  a15=(a1*a15);
  a10=(a12*a10);
  a10=(a1*a10);
  a15=(a15+a10);
  a15=(a0*a15);
  a15=(-a15);
  if (res[0]!=0) res[0][4]=a15;
  a20=(a2*a20);
  a4=(a4+a20);
  a4=(-a4);
  if (res[0]!=0) res[0][5]=a4;
  a4=(a11*a17);
  a20=(a5*a18);
  a4=(a4-a20);
  a20=cos(a13);
  a4=(a4*a20);
  a11=(a11*a9);
  a5=(a5*a7);
  a11=(a11+a5);
  a13=sin(a13);
  a11=(a11*a13);
  a4=(a4-a11);
  a11=(a12*a4);
  a11=(a1*a11);
  a11=(a0*a11);
  a11=(-a11);
  if (res[0]!=0) res[0][6]=a11;
  a7=(a21*a7);
  a9=(a3*a9);
  a7=(a7-a9);
  a7=(a7*a13);
  a3=(a3*a17);
  a21=(a21*a18);
  a3=(a3+a21);
  a3=(a3*a20);
  a7=(a7+a3);
  a12=(a12*a7);
  a1=(a1*a12);
  a0=(a0*a1);
  a0=(-a0);
  if (res[0]!=0) res[0][7]=a0;
  a2=(a2*a4);
  a2=(-a2);
  if (res[0]!=0) res[0][8]=a2;
  return 0;
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_TaskJacobian(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f1(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_TaskJacobian_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_TaskJacobian_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_TaskJacobian_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_TaskJacobian_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_TaskJacobian_release(int mem) {
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_TaskJacobian_incref(void) {
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_TaskJacobian_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int g_InverseKinematics_TaskJacobian_n_in(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_int g_InverseKinematics_TaskJacobian_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real g_InverseKinematics_TaskJacobian_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_InverseKinematics_TaskJacobian_name_in(casadi_int i){
  switch (i) {
    case 0: return "q";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_InverseKinematics_TaskJacobian_name_out(casadi_int i){
  switch (i) {
    case 0: return "TaskJacobian";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_InverseKinematics_TaskJacobian_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_InverseKinematics_TaskJacobian_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_TaskJacobian_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 1;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

/* g_InverseKinematics_TaskJacobian_derivative:(q[3],v[3])->(TaskJacobian_derivative[4x3,9nz]) */
static int casadi_f2(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a4, a40, a41, a42, a43, a44, a5, a6, a7, a8, a9;
  a0=3.3333333333333333e-02;
  a1=10.;
  a2=5.0000000000000000e-01;
  a3=arg[0]? arg[0][0] : 0;
  a4=sin(a3);
  a5=(a2*a4);
  a6=arg[0]? arg[0][1] : 0;
  a7=cos(a6);
  a8=(a7*a4);
  a9=sin(a6);
  a10=cos(a3);
  a11=(a9*a10);
  a8=(a8+a11);
  a11=(a2*a8);
  a11=(a5+a11);
  a12=2.5000000000000000e-01;
  a13=arg[0]? arg[0][2] : 0;
  a14=cos(a13);
  a15=(a14*a8);
  a16=sin(a13);
  a17=cos(a6);
  a10=(a17*a10);
  a18=sin(a6);
  a19=(a18*a4);
  a10=(a10-a19);
  a10=(a16*a10);
  a15=(a15+a10);
  a10=(a12*a15);
  a10=(a11+a10);
  a10=(a1*a10);
  a8=(a12*a8);
  a5=(a5+a8);
  a5=(a1*a5);
  a4=(a12*a4);
  a4=(a1*a4);
  a5=(a5+a4);
  a10=(a10+a5);
  a10=(a0*a10);
  a5=arg[1]? arg[1][0] : 0;
  a10=(a10*a5);
  a4=cos(a3);
  a8=sin(a6);
  a19=(a4*a8);
  a20=sin(a3);
  a21=cos(a6);
  a22=(a20*a21);
  a19=(a19+a22);
  a22=(a2*a19);
  a23=(a14*a19);
  a24=cos(a6);
  a25=(a4*a24);
  a26=sin(a6);
  a27=(a20*a26);
  a25=(a25-a27);
  a25=(a16*a25);
  a23=(a23+a25);
  a25=(a12*a23);
  a25=(a22+a25);
  a25=(a1*a25);
  a19=(a12*a19);
  a19=(a1*a19);
  a25=(a25+a19);
  a25=(a0*a25);
  a19=arg[1]? arg[1][1] : 0;
  a25=(a25*a19);
  a10=(a10+a25);
  a25=(a7*a4);
  a27=(a9*a20);
  a25=(a25-a27);
  a27=sin(a13);
  a25=(a25*a27);
  a20=(a17*a20);
  a4=(a18*a4);
  a20=(a20+a4);
  a4=cos(a13);
  a20=(a20*a4);
  a25=(a25+a20);
  a20=(a12*a25);
  a20=(a1*a20);
  a20=(a0*a20);
  a28=arg[1]? arg[1][2] : 0;
  a20=(a20*a28);
  a10=(a10+a20);
  if (res[0]!=0) res[0][0]=a10;
  a10=cos(a3);
  a20=(a2*a10);
  a29=(a7*a10);
  a30=sin(a3);
  a31=(a9*a30);
  a29=(a29-a31);
  a31=(a2*a29);
  a31=(a20+a31);
  a32=(a14*a29);
  a30=(a17*a30);
  a33=(a18*a10);
  a30=(a30+a33);
  a30=(a16*a30);
  a32=(a32-a30);
  a32=(a12*a32);
  a31=(a31+a32);
  a31=(a1*a31);
  a29=(a12*a29);
  a20=(a20+a29);
  a20=(a1*a20);
  a10=(a12*a10);
  a10=(a1*a10);
  a20=(a20+a10);
  a31=(a31+a20);
  a31=(a0*a31);
  a31=(a31*a5);
  a20=cos(a3);
  a10=(a20*a21);
  a29=sin(a3);
  a32=(a29*a8);
  a10=(a10-a32);
  a32=(a2*a10);
  a30=(a14*a10);
  a33=(a20*a26);
  a34=(a29*a24);
  a33=(a33+a34);
  a33=(a16*a33);
  a30=(a30-a33);
  a30=(a12*a30);
  a32=(a32+a30);
  a32=(a1*a32);
  a10=(a12*a10);
  a10=(a1*a10);
  a32=(a32+a10);
  a32=(a0*a32);
  a32=(a32*a19);
  a31=(a31+a32);
  a32=(a17*a20);
  a10=(a18*a29);
  a32=(a32-a10);
  a32=(a32*a4);
  a29=(a7*a29);
  a20=(a9*a20);
  a29=(a29+a20);
  a29=(a29*a27);
  a32=(a32-a29);
  a32=(a12*a32);
  a32=(a1*a32);
  a32=(a0*a32);
  a32=(a32*a28);
  a31=(a31+a32);
  a31=(-a31);
  if (res[0]!=0) res[0][1]=a31;
  a15=(a2*a15);
  a11=(a11+a15);
  a11=(a11*a5);
  a23=(a2*a23);
  a22=(a22+a23);
  a22=(a22*a19);
  a11=(a11+a22);
  a25=(a2*a25);
  a25=(a25*a28);
  a11=(a11+a25);
  if (res[0]!=0) res[0][2]=a11;
  a11=cos(a6);
  a25=sin(a3);
  a22=(a11*a25);
  a23=sin(a6);
  a15=cos(a3);
  a31=(a23*a15);
  a22=(a22+a31);
  a31=(a2*a22);
  a32=(a14*a22);
  a29=cos(a6);
  a20=(a29*a15);
  a10=sin(a6);
  a30=(a10*a25);
  a20=(a20-a30);
  a20=(a16*a20);
  a32=(a32+a20);
  a20=(a12*a32);
  a20=(a31+a20);
  a20=(a1*a20);
  a22=(a12*a22);
  a22=(a1*a22);
  a20=(a20+a22);
  a20=(a0*a20);
  a20=(a20*a5);
  a22=cos(a3);
  a30=sin(a6);
  a33=(a22*a30);
  a34=sin(a3);
  a35=cos(a6);
  a36=(a34*a35);
  a33=(a33+a36);
  a36=(a2*a33);
  a37=(a14*a33);
  a38=cos(a6);
  a39=(a22*a38);
  a6=sin(a6);
  a40=(a34*a6);
  a39=(a39-a40);
  a39=(a16*a39);
  a37=(a37+a39);
  a39=(a12*a37);
  a39=(a36+a39);
  a39=(a1*a39);
  a33=(a12*a33);
  a33=(a1*a33);
  a39=(a39+a33);
  a39=(a0*a39);
  a39=(a39*a19);
  a20=(a20+a39);
  a39=(a22*a11);
  a33=(a34*a23);
  a39=(a39-a33);
  a39=(a39*a27);
  a33=(a22*a10);
  a40=(a34*a29);
  a33=(a33+a40);
  a33=(a33*a4);
  a39=(a39+a33);
  a33=(a12*a39);
  a33=(a1*a33);
  a33=(a0*a33);
  a33=(a33*a28);
  a20=(a20+a33);
  if (res[0]!=0) res[0][3]=a20;
  a20=cos(a3);
  a33=(a11*a20);
  a40=sin(a3);
  a41=(a23*a40);
  a33=(a33-a41);
  a41=(a2*a33);
  a42=(a14*a33);
  a43=(a29*a40);
  a44=(a10*a20);
  a43=(a43+a44);
  a43=(a16*a43);
  a42=(a42-a43);
  a42=(a12*a42);
  a41=(a41+a42);
  a41=(a1*a41);
  a33=(a12*a33);
  a33=(a1*a33);
  a41=(a41+a33);
  a41=(a0*a41);
  a41=(a41*a5);
  a33=cos(a3);
  a35=(a33*a35);
  a3=sin(a3);
  a30=(a3*a30);
  a35=(a35-a30);
  a30=(a2*a35);
  a14=(a14*a35);
  a6=(a33*a6);
  a38=(a3*a38);
  a6=(a6+a38);
  a16=(a16*a6);
  a14=(a14-a16);
  a14=(a12*a14);
  a30=(a30+a14);
  a30=(a1*a30);
  a35=(a12*a35);
  a35=(a1*a35);
  a30=(a30+a35);
  a30=(a0*a30);
  a30=(a30*a19);
  a41=(a41+a30);
  a29=(a33*a29);
  a10=(a3*a10);
  a29=(a29-a10);
  a29=(a29*a4);
  a23=(a33*a23);
  a11=(a3*a11);
  a23=(a23+a11);
  a23=(a23*a27);
  a29=(a29-a23);
  a29=(a12*a29);
  a29=(a1*a29);
  a29=(a0*a29);
  a29=(a29*a28);
  a41=(a41+a29);
  a41=(-a41);
  if (res[0]!=0) res[0][4]=a41;
  a32=(a2*a32);
  a31=(a31+a32);
  a31=(a31*a5);
  a37=(a2*a37);
  a36=(a36+a37);
  a36=(a36*a19);
  a31=(a31+a36);
  a39=(a2*a39);
  a39=(a39*a28);
  a31=(a31+a39);
  if (res[0]!=0) res[0][5]=a31;
  a31=cos(a13);
  a39=(a17*a25);
  a36=(a18*a15);
  a39=(a39+a36);
  a39=(a31*a39);
  a36=sin(a13);
  a15=(a7*a15);
  a25=(a9*a25);
  a15=(a15-a25);
  a15=(a36*a15);
  a39=(a39+a15);
  a15=(a12*a39);
  a15=(a1*a15);
  a15=(a0*a15);
  a15=(a15*a5);
  a25=(a22*a26);
  a37=(a34*a24);
  a25=(a25+a37);
  a25=(a31*a25);
  a37=(a22*a21);
  a32=(a34*a8);
  a37=(a37-a32);
  a37=(a36*a37);
  a25=(a25+a37);
  a37=(a12*a25);
  a37=(a1*a37);
  a37=(a0*a37);
  a37=(a37*a19);
  a15=(a15+a37);
  a37=(a22*a17);
  a32=(a34*a18);
  a37=(a37-a32);
  a32=sin(a13);
  a37=(a37*a32);
  a22=(a22*a9);
  a34=(a34*a7);
  a22=(a22+a34);
  a13=cos(a13);
  a22=(a22*a13);
  a37=(a37+a22);
  a22=(a12*a37);
  a22=(a1*a22);
  a22=(a0*a22);
  a22=(a22*a28);
  a15=(a15+a22);
  if (res[0]!=0) res[0][6]=a15;
  a15=(a17*a20);
  a22=(a18*a40);
  a15=(a15-a22);
  a15=(a31*a15);
  a40=(a7*a40);
  a20=(a9*a20);
  a40=(a40+a20);
  a40=(a36*a40);
  a15=(a15-a40);
  a15=(a12*a15);
  a15=(a1*a15);
  a15=(a0*a15);
  a15=(a15*a5);
  a24=(a33*a24);
  a26=(a3*a26);
  a24=(a24-a26);
  a31=(a31*a24);
  a8=(a33*a8);
  a21=(a3*a21);
  a8=(a8+a21);
  a36=(a36*a8);
  a31=(a31-a36);
  a31=(a12*a31);
  a31=(a1*a31);
  a31=(a0*a31);
  a31=(a31*a19);
  a15=(a15+a31);
  a7=(a33*a7);
  a9=(a3*a9);
  a7=(a7-a9);
  a7=(a7*a13);
  a3=(a3*a17);
  a33=(a33*a18);
  a3=(a3+a33);
  a3=(a3*a32);
  a7=(a7-a3);
  a12=(a12*a7);
  a1=(a1*a12);
  a0=(a0*a1);
  a0=(a0*a28);
  a15=(a15+a0);
  a15=(-a15);
  if (res[0]!=0) res[0][7]=a15;
  a39=(a2*a39);
  a39=(a39*a5);
  a25=(a2*a25);
  a25=(a25*a19);
  a39=(a39+a25);
  a2=(a2*a37);
  a2=(a2*a28);
  a39=(a39+a2);
  if (res[0]!=0) res[0][8]=a39;
  return 0;
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_TaskJacobian_derivative(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f2(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_TaskJacobian_derivative_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_TaskJacobian_derivative_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_TaskJacobian_derivative_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_TaskJacobian_derivative_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_TaskJacobian_derivative_release(int mem) {
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_TaskJacobian_derivative_incref(void) {
}

CASADI_SYMBOL_EXPORT void g_InverseKinematics_TaskJacobian_derivative_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int g_InverseKinematics_TaskJacobian_derivative_n_in(void) { return 2;}

CASADI_SYMBOL_EXPORT casadi_int g_InverseKinematics_TaskJacobian_derivative_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real g_InverseKinematics_TaskJacobian_derivative_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_InverseKinematics_TaskJacobian_derivative_name_in(casadi_int i){
  switch (i) {
    case 0: return "q";
    case 1: return "v";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_InverseKinematics_TaskJacobian_derivative_name_out(casadi_int i){
  switch (i) {
    case 0: return "TaskJacobian_derivative";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_InverseKinematics_TaskJacobian_derivative_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_InverseKinematics_TaskJacobian_derivative_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int g_InverseKinematics_TaskJacobian_derivative_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 2;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
