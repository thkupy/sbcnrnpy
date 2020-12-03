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
 
#define nrn_init _nrn_init__ka
#define _nrn_initial _nrn_initial__ka
#define nrn_cur _nrn_cur__ka
#define _nrn_current _nrn_current__ka
#define nrn_jacob _nrn_jacob__ka
#define nrn_state _nrn_state__ka
#define _net_receive _net_receive__ka 
#define _f_trates _f_trates__ka 
#define rates rates__ka 
#define states states__ka 
#define trates trates__ka 
 
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
#define gkabar _p[0]
#define ik _p[1]
#define gka _p[2]
#define a _p[3]
#define b _p[4]
#define c _p[5]
#define ek _p[6]
#define Da _p[7]
#define Db _p[8]
#define Dc _p[9]
#define _g _p[10]
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
 "setdata_ka", _hoc_setdata,
 "rates_ka", _hoc_rates,
 "states_ka", _hoc_states,
 "trates_ka", _hoc_trates,
 "vtrap_ka", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_ka
 extern double vtrap( double , double );
 /* declare global and static user variables */
#define atau atau_ka
 double atau = 0;
#define ainf ainf_ka
 double ainf = 0;
#define btau btau_ka
 double btau = 0;
#define binf binf_ka
 double binf = 0;
#define ctau ctau_ka
 double ctau = 0;
#define cinf cinf_ka
 double cinf = 0;
#define usetable usetable_ka
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gkabar_ka", 0, 1e+09,
 "usetable_ka", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "atau_ka", "ms",
 "btau_ka", "ms",
 "ctau_ka", "ms",
 "gkabar_ka", "mho/cm2",
 "ik_ka", "mA/cm2",
 "gka_ka", "mho/cm2",
 0,0
};
 static double a0 = 0;
 static double b0 = 0;
 static double c0 = 0;
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ainf_ka", &ainf_ka,
 "binf_ka", &binf_ka,
 "cinf_ka", &cinf_ka,
 "atau_ka", &atau_ka,
 "btau_ka", &btau_ka,
 "ctau_ka", &ctau_ka,
 "usetable_ka", &usetable_ka,
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
"ka",
 "gkabar_ka",
 0,
 "ik_ka",
 "gka_ka",
 0,
 "a_ka",
 "b_ka",
 "c_ka",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 11, _prop);
 	/*initialize range parameters*/
 	gkabar = 0.00477;
 	_prop->param = _p;
 	_prop->param_size = 11;
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

 void _ka_reg() {
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
  hoc_register_prop_size(_mechtype, 11, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ka /home/kuenzel/Dokumente/Python/smallexc/mechanisms/x86_64/ka.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zaexp , _zbexp , _zcexp ;
 static double _zq10 ;
 static double *_t_ainf;
 static double *_t__zaexp;
 static double *_t_binf;
 static double *_t__zbexp;
 static double *_t_cinf;
 static double *_t__zcexp;
static int _reset;
static char *modelname = "klt.mod  The low threshold conductance of cochlear nucleus neurons";

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
   a = a + _zaexp * ( ainf - a ) ;
   b = b + _zbexp * ( binf - b ) ;
   c = c + _zcexp * ( cinf - c ) ;
   
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
   ainf = pow( ( 1.0 / ( 1.0 + exp ( - 1.0 * ( _lv + 31.0 ) / 6.0 ) ) ) , 0.25 ) ;
   binf = 1.0 / pow( ( 1.0 + exp ( ( _lv + 66.0 ) / 7.0 ) ) , 0.5 ) ;
   cinf = 1.0 / pow( ( 1.0 + exp ( ( _lv + 66.0 ) / 7.0 ) ) , 0.5 ) ;
   atau = ( 100.0 / ( 7.0 * exp ( ( _lv + 60.0 ) / 14.0 ) + 29.0 * exp ( - ( _lv + 60.0 ) / 24.0 ) ) ) + 0.1 ;
   btau = ( 1000.0 / ( 14.0 * exp ( ( _lv + 60.0 ) / 27.0 ) + 29.0 * exp ( - ( _lv + 60.0 ) / 24.0 ) ) ) + 1.0 ;
   ctau = ( 90.0 / ( 1.0 + exp ( ( - 66.0 - _lv ) / 17.0 ) ) ) + 10.0 ;
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
    _t_ainf[_i] = ainf;
    _t__zaexp[_i] = _zaexp;
    _t_binf[_i] = binf;
    _t__zbexp[_i] = _zbexp;
    _t_cinf[_i] = cinf;
    _t__zcexp[_i] = _zcexp;
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
  ainf = _xi;
  _zaexp = _xi;
  binf = _xi;
  _zbexp = _xi;
  cinf = _xi;
  _zcexp = _xi;
  return;
 }
 if (_xi <= 0.) {
 ainf = _t_ainf[0];
 _zaexp = _t__zaexp[0];
 binf = _t_binf[0];
 _zbexp = _t__zbexp[0];
 cinf = _t_cinf[0];
 _zcexp = _t__zcexp[0];
 return; }
 if (_xi >= 300.) {
 ainf = _t_ainf[300];
 _zaexp = _t__zaexp[300];
 binf = _t_binf[300];
 _zbexp = _t__zbexp[300];
 cinf = _t_cinf[300];
 _zcexp = _t__zcexp[300];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 ainf = _t_ainf[_i] + _theta*(_t_ainf[_i+1] - _t_ainf[_i]);
 _zaexp = _t__zaexp[_i] + _theta*(_t__zaexp[_i+1] - _t__zaexp[_i]);
 binf = _t_binf[_i] + _theta*(_t_binf[_i+1] - _t_binf[_i]);
 _zbexp = _t__zbexp[_i] + _theta*(_t__zbexp[_i+1] - _t__zbexp[_i]);
 cinf = _t_cinf[_i] + _theta*(_t_cinf[_i+1] - _t_cinf[_i]);
 _zcexp = _t__zcexp[_i] + _theta*(_t__zcexp[_i+1] - _t__zcexp[_i]);
 }

 
static int  _f_trates (  double _lv ) {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   _zaexp = 1.0 - exp ( _ltinc / atau ) ;
   _zbexp = 1.0 - exp ( _ltinc / btau ) ;
   _zcexp = 1.0 - exp ( _ltinc / ctau ) ;
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
 
static int _ode_count(int _type){ hoc_execerror("ka", "cannot be used with CVODE"); return 0;}
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
  a = a0;
  b = b0;
  c = c0;
 {
   trates ( _threadargscomma_ v ) ;
   a = ainf ;
   b = binf ;
   c = cinf ;
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
   gka = gkabar * ( pow( a , 4.0 ) ) * b * c ;
   ik = gka * ( v - ek ) ;
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
 if(error){fprintf(stderr,"at line 61 in file ka.mod:\n    \n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_ainf = makevector(301*sizeof(double));
   _t__zaexp = makevector(301*sizeof(double));
   _t_binf = makevector(301*sizeof(double));
   _t__zbexp = makevector(301*sizeof(double));
   _t_cinf = makevector(301*sizeof(double));
   _t__zcexp = makevector(301*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/kuenzel/Dokumente/Python/smallexc/mechanisms/ka.mod";
static const char* nmodl_file_text = 
  "TITLE klt.mod  The low threshold conductance of cochlear nucleus neurons\n"
  "\n"
  "COMMENT\n"
  "\n"
  "NEURON implementation of Jason Rothman's measurements of VCN conductances.\n"
  "\n"
  "This file implements the transient potassium current found in ventral cochlear\n"
  "nucleus \"Type I\" cells, which are largely \"stellate\" or \"multipolar\" cells  (Manis and\n"
  "Marx, 1991; Rothman and Manis, 2003a,b; Manis et al, 1996). The current is likely\n"
  " mediated by Kv4.2 potassium channel subunits, but this has not been directly\n"
  "demonstrated. The specific implementation is described in Rothman and Manis, J.\n"
  "Neurophysiol. 2003, in the appendix. Measurements were made from isolated \n"
  "neurons from adult guinea pig, under reasonably stringent voltage clamp conditions.\n"
  " The measured current is sensitive to 4-aminopyridine. \n"
  "Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.\n"
  "\n"
  "File split implementaiton, April 1, 2004.\n"
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
  "        SUFFIX ka\n"
  "        USEION k READ ek WRITE ik\n"
  "        RANGE gkabar, gka, ik\n"
  "        GLOBAL ainf, binf, cinf, atau, btau, ctau\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "PARAMETER {\n"
  "        v (mV)\n"
  "        celsius = 22 (degC)  : model is defined on measurements made at room temp in Baltimore\n"
  "        dt (ms)\n"
  "        ek = -77 (mV)\n"
  "        gkabar = 0.00477 (mho/cm2) <0,1e9>\n"
  "}\n"
  "\n"
  "STATE {\n"
  "        a b c\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    ik (mA/cm2) \n"
  "    gka (mho/cm2)\n"
  "    ainf binf cinf\n"
  "    atau (ms) btau (ms) ctau (ms)\n"
  "    }\n"
  "\n"
  "LOCAL aexp, bexp, cexp\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "    \n"
  "	gka = gkabar*(a^4)*b*c\n"
  "    ik = gka*(v - ek)\n"
  "\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL {\n"
  "    trates(v)\n"
  "    a = ainf\n"
  "    b = binf\n"
  "    c = cinf\n"
  "}\n"
  "\n"
  "PROCEDURE states() {  :Computes state variables m, h, and n\n"
  "	trates(v)      :             at the current v and dt.\n"
  "	a = a + aexp*(ainf-a)\n"
  "	b = b + bexp*(binf-b)\n"
  "	c = c + cexp*(cinf-c)\n"
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
  "    ainf = (1 / (1 + exp(-1*(v + 31) / 6)))^0.25\n"
  "    binf = 1 / (1 + exp((v + 66) / 7))^0.5\n"
  "    cinf = 1 / (1 + exp((v + 66) / 7))^0.5\n"
  "\n"
  "    atau =  (100 / (7*exp((v+60) / 14) + 29*exp(-(v+60) / 24))) + 0.1\n"
  "    btau =  (1000 / (14*exp((v+60) / 27) + 29*exp(-(v+60) / 24))) + 1\n"
  "    ctau = (90 / (1 + exp((-66-v) / 17))) + 10\n"
  "}\n"
  "\n"
  "PROCEDURE trates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "	LOCAL tinc\n"
  "	TABLE ainf, aexp, binf, bexp, cinf, cexp\n"
  "	DEPEND dt, celsius FROM -150 TO 150 WITH 300\n"
  "\n"
  "    rates(v)    : not consistently executed from here if usetable_hh == 1\n"
  "        : so don't expect the tau values to be tracking along with\n"
  "        : the inf values in hoc\n"
  "\n"
  "	tinc = -dt * q10\n"
  "	aexp = 1 - exp(tinc/atau)\n"
  "	bexp = 1 - exp(tinc/btau)\n"
  "	cexp = 1 - exp(tinc/ctau)\n"
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
