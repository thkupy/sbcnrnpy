/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__kht
#define _nrn_initial _nrn_initial__kht
#define nrn_cur _nrn_cur__kht
#define _nrn_current _nrn_current__kht
#define nrn_jacob _nrn_jacob__kht
#define nrn_state _nrn_state__kht
#define _net_receive _net_receive__kht 
#define _f_trates _f_trates__kht 
#define rates rates__kht 
#define states states__kht 
#define trates trates__kht 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gkhtbar _p[0]
#define ik _p[1]
#define gkht _p[2]
#define n _p[3]
#define p _p[4]
#define ek _p[5]
#define Dn _p[6]
#define Dp _p[7]
#define _g _p[8]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_states(void);
 static void _hoc_trates(void);
 static void _hoc_vtrap(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_kht", _hoc_setdata,
 "rates_kht", _hoc_rates,
 "states_kht", _hoc_states,
 "trates_kht", _hoc_trates,
 "vtrap_kht", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_kht
 extern double vtrap( double , double );
 /* declare global and static user variables */
#define nf nf_kht
 double nf = 0.85;
#define ntau ntau_kht
 double ntau = 0;
#define ninf ninf_kht
 double ninf = 0;
#define ptau ptau_kht
 double ptau = 0;
#define pinf pinf_kht
 double pinf = 0;
#define usetable usetable_kht
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gkhtbar_kht", 0, 1e+09,
 "nf_kht", 0, 1,
 "usetable_kht", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ptau_kht", "ms",
 "ntau_kht", "ms",
 "gkhtbar_kht", "mho/cm2",
 "ik_kht", "mA/cm",
 "gkht_kht", "mho/cm2",
 0,0
};
 static double delta_t = 1;
 static double n0 = 0;
 static double p0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "nf_kht", &nf_kht,
 "pinf_kht", &pinf_kht,
 "ninf_kht", &ninf_kht,
 "ptau_kht", &ptau_kht,
 "ntau_kht", &ntau_kht,
 "usetable_kht", &usetable_kht,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"kht",
 "gkhtbar_kht",
 0,
 "ik_kht",
 "gkht_kht",
 0,
 "n_kht",
 "p_kht",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 9, _prop);
 	/*initialize range parameters*/
 	gkhtbar = 0.01592;
 	_prop->param = _p;
 	_prop->param_size = 9;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _kht_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 9, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 kht /home/kuenzel/Dokumente/Python/smallexc/mechanisms/x86_64/kht.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _znexp , _zpexp ;
 static double _zq10 ;
 static double *_t_ninf;
 static double *_t__znexp;
 static double *_t_pinf;
 static double *_t__zpexp;
static int _reset;
static char *modelname = "kht.mod  The high threshold conductance of cochlear nucleus neurons";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_trates(double);
static int rates(double);
static int states();
static int trates(double);
 static void _n_trates(double);
 
static int  states (  ) {
   trates ( _threadargscomma_ v ) ;
   n = n + _znexp * ( ninf - n ) ;
   p = p + _zpexp * ( pinf - p ) ;
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv ) {
   _zq10 = pow( 3.0 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   ninf = pow( ( 1.0 + exp ( - ( _lv + 15.0 ) / 5.0 ) ) , - 0.5 ) ;
   pinf = 1.0 / ( 1.0 + exp ( - ( _lv + 23.0 ) / 6.0 ) ) ;
   ntau = ( 100.0 / ( 11.0 * exp ( ( _lv + 60.0 ) / 24.0 ) + 21.0 * exp ( - ( _lv + 60.0 ) / 23.0 ) ) ) + 0.7 ;
   ptau = ( 100.0 / ( 4.0 * exp ( ( _lv + 60.0 ) / 32.0 ) + 5.0 * exp ( - ( _lv + 60.0 ) / 22.0 ) ) ) + 5.0 ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_trates, _tmin_trates;
 static void _check_trates();
 static void _check_trates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_trates)/300.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 301; _x += _dx, _i++) {
    _f_trates(_x);
    _t_ninf[_i] = ninf;
    _t__znexp[_i] = _znexp;
    _t_pinf[_i] = pinf;
    _t__zpexp[_i] = _zpexp;
   }
   _sav_dt = dt;
   _sav_celsius = celsius;
  }
 }

 static int trates(double _lv){ _check_trates();
 _n_trates(_lv);
 return 0;
 }

 static void _n_trates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 if (isnan(_xi)) {
  ninf = _xi;
  _znexp = _xi;
  pinf = _xi;
  _zpexp = _xi;
  return;
 }
 if (_xi <= 0.) {
 ninf = _t_ninf[0];
 _znexp = _t__znexp[0];
 pinf = _t_pinf[0];
 _zpexp = _t__zpexp[0];
 return; }
 if (_xi >= 300.) {
 ninf = _t_ninf[300];
 _znexp = _t__znexp[300];
 pinf = _t_pinf[300];
 _zpexp = _t__zpexp[300];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 _znexp = _t__znexp[_i] + _theta*(_t__znexp[_i+1] - _t__znexp[_i]);
 pinf = _t_pinf[_i] + _theta*(_t_pinf[_i+1] - _t_pinf[_i]);
 _zpexp = _t__zpexp[_i] + _theta*(_t__zpexp[_i+1] - _t__zpexp[_i]);
 }

 
static int  _f_trates (  double _lv ) {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   _znexp = 1.0 - exp ( _ltinc / ntau ) ;
   _zpexp = 1.0 - exp ( _ltinc / ptau ) ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
    _r = 1.;
 trates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap (  double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   _r =  vtrap (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("kht", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  n = n0;
  p = p0;
 {
   trates ( _threadargscomma_ v ) ;
   p = pinf ;
   n = ninf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gkht = gkhtbar * ( nf * ( pow( n , 2.0 ) ) + ( 1.0 - nf ) * p ) ;
   ik = gkht * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ek = _ion_ek;
 { error =  states();
 if(error){fprintf(stderr,"at line 74 in file kht.mod:\n    \n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_ninf = makevector(301*sizeof(double));
   _t__znexp = makevector(301*sizeof(double));
   _t_pinf = makevector(301*sizeof(double));
   _t__zpexp = makevector(301*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/kuenzel/Dokumente/Python/smallexc/mechanisms/kht.mod";
static const char* nmodl_file_text = 
  "TITLE kht.mod  The high threshold conductance of cochlear nucleus neurons\n"
  "\n"
  "COMMENT\n"
  "\n"
  "NEURON implementation of Jason Rothman's measurements of VCN conductances.\n"
  "\n"
  "This file implements the high threshold potassium current found in several brainstem\n"
  " nuclei of the auditory system, including the spherical and globular bushy cells\n"
  "  (Manis and Marx, 1991; Rothman and Manis, 2003a,b) and multipolar (stellate) \n"
  "  cells of the ventral cochlear nucleus, principal cells of the medial \n"
  "  nucleus of the trapzoid body (Brew and Forsythe, 1995, Wang and Kaczmarek, \n"
  "  1997) and neurons of the medial superior olive. The current is likely mediated by \n"
  "  Kv3.1  potassium channel subunits. The specific \n"
  "  implementation is described in Rothman and Manis, J. Neurophysiol. 2003, in the \n"
  "  appendix. Measurements were made from isolated neurons from adult guinea pig, \n"
  "  under reasonably stringent voltage clamp conditions. The measured current is \n"
  "  sensitive to 4-aminopyridine and TEA, but is spared by mamba snake toxi\n"
  "  dendrotoxin I.\n"
  "\n"
  "\n"
  "Similar conductrances are found in the homologous neurons of the avian auditory \n"
  "system (Reyes and Rubel; Zhang and Trussell; Rathouz and Trussell), and the \n"
  "conductance described here, in the absence of more detailed kinetic measurements\n"
  ", is probably suitable for use in modeling that system.\n"
  "\n"
  "\n"
  "Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.\n"
  "\n"
  "File split implementation, February 28, 2004.\n"
  "\n"
  "Contact: pmanis@med.unc.edu\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "UNITS {\n"
  "        (mA) = (milliamp)\n"
  "        (mV) = (millivolt)\n"
  "        (nA) = (nanoamp)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "        SUFFIX kht\n"
  "        USEION k READ ek WRITE ik\n"
  "        RANGE gkhtbar, gkht, ik\n"
  "        GLOBAL ninf, pinf, ntau, ptau\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "PARAMETER {\n"
  "        v (mV)\n"
  "        celsius = 22 (degC)  : model is defined on measurements made at room temp in Baltimore\n"
  "        dt (ms)\n"
  "        ek = -77 (mV)\n"
  "        gkhtbar = 0.01592 (mho/cm2) <0,1e9>\n"
  "		nf = 0.85 <0,1> :proportion of n vs p kinetics\n"
  "}\n"
  "\n"
  "STATE {\n"
  "        n p\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    ik (mA/cm) \n"
  "    gkht (mho/cm2)\n"
  "    pinf ninf\n"
  "    ptau (ms) ntau (ms)\n"
  "    }\n"
  "\n"
  "LOCAL nexp, pexp\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "    \n"
  "	gkht = gkhtbar*(nf*(n^2) + (1-nf)*p)\n"
  "    ik = gkht*(v - ek)\n"
  "\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL {\n"
  "    trates(v)\n"
  "    p = pinf\n"
  "    n = ninf\n"
  "}\n"
  "\n"
  "PROCEDURE states() {  :Computes state variables m, h, and n\n"
  "	trates(v)      :             at the current v and dt.\n"
  "	n = n + nexp*(ninf-n)\n"
  "	p = p + pexp*(pinf-p)\n"
  "VERBATIM\n"
  "	return 0;\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "LOCAL q10\n"
  "\n"
  "PROCEDURE rates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "\n"
  "	q10 = 3^((celsius - 22)/10) : if you don't like room temp, it can be changed!\n"
  "\n"
  "    ninf =   (1 + exp(-(v + 15) / 5))^-0.5\n"
  "    pinf =  1 / (1 + exp(-(v + 23) / 6))\n"
  "\n"
  "	ntau =  (100 / (11*exp((v+60) / 24) + 21*exp(-(v+60) / 23))) + 0.7\n"
  "    ptau = (100 / (4*exp((v+60) / 32) + 5*exp(-(v+60) / 22))) + 5\n"
  "}\n"
  "\n"
  "PROCEDURE trates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "	LOCAL tinc\n"
  "	TABLE ninf, nexp, pinf, pexp\n"
  "	DEPEND dt, celsius FROM -150 TO 150 WITH 300\n"
  "\n"
  "    rates(v)    : not consistently executed from here if usetable_hh == 1\n"
  "        : so don't expect the tau values to be tracking along with\n"
  "        : the inf values in hoc\n"
  "\n"
  "	tinc = -dt * q10\n"
  "	nexp = 1 - exp(tinc/ntau)\n"
  "	pexp = 1 - exp(tinc/ptau)\n"
  "	}\n"
  "\n"
  "FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.\n"
  "        if (fabs(x/y) < 1e-6) {\n"
  "                vtrap = y*(1 - x/y/2)\n"
  "        }else{\n"
  "                vtrap = x/(exp(x/y) - 1)\n"
  "        }\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
