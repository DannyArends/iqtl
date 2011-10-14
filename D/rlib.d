/**
 * Functions exported by Rlib
 */

import std.stdio;
import std.conv;

alias uint SEXPTYPE;
alias int R_len_t;

struct sxpinfo_struct{
  SEXPTYPE type      =  5;
  uint obj   =  1;
  uint named =  2;
  uint gp    = 16;
  uint mark  =  1;
  uint deb   =  1;
  uint trace =  1;  /* functions and memory tracing */
  uint spare =  1;  /* currently unused */
  uint gcgen =  1;  /* old generation number */
  uint gccls =  3;  /* node class */
};

struct symsxp_struct {
  SEXPREC *pname;
  SEXPREC *value;
  SEXPREC *internal;
};

struct listsxp_struct {
  SEXPREC* carval;
  SEXPREC* cdrval;
  SEXPREC* tagval;
};

struct envsxp_struct {
  SEXPREC* frame;
  SEXPREC* enclos;
  SEXPREC* hashtab;
};

struct closxp_struct {
  SEXPREC* formals;
  SEXPREC* b;
  SEXPREC* env;
};

struct promsxp_struct {
  SEXPREC* value;
  SEXPREC* expr;
  SEXPREC* env;
};

struct vecsxp_struct {
    R_len_t	length;
    R_len_t	truelength;
};

struct primsxp_struct {
    int offset;
};

struct SEXPREC_HEADER{
  union U{
    sxpinfo_struct sxpinfo;
    SEXPREC* attrib;
    SEXPREC* gengc_next_node;
    SEXPREC* gengc_prev_node;
  };
};

struct SEXPREC {
  SEXPREC_HEADER header;
  union U {
	  primsxp_struct primsxp;
	  symsxp_struct symsxp;
	  listsxp_struct listsxp;
	  envsxp_struct envsxp;
	  closxp_struct closxp;
	  promsxp_struct promsxp;
  };
};

alias SEXPREC *SEXP;

const NILSXP = 0;
const SYMSXP = 1;
const LISTSXP = 2;
const CLOSXP = 3;
const ENVSXP = 4;
const PROMSXP = 5;
const LANGSXP = 6;
const SPECIALSXP = 7;
const BUILTINSXP = 8;
const CHARSXP = 9;
const LGLSXP = 10;
const INTSXP = 13;
const REALSXP = 14;
const CPLXSXP = 15;
const STRSXP = 16;
const ANYSXP = 18;
const VECSXP = 19;
const EXPRSXP = 20;
const BCODESXP = 21;
const EXTPTRSXP = 22;
const WEAKREFSXP = 23;
const RAWSXP = 24;
const S4SXP = 25;
const FUNSXP = 99;

version (Windows) {
  import libload;
  import std.loader;
  
  SEXPREC* R_GlobalEnv;
  
  extern(C){
    double function(double, double, double, int) dnorm;
    double function(double, double, double, int, int) qf;

    // Here we have the tricky part, we need to combine calling these 3 functions to produce a new seed
    void   function() R_SeedsSymbol; 
    void   function(int*) seed_in;
    void   function(int*) seed_out;
       
    double function() norm_rand;
    void   function(char *, ...) Rprintf;
    void   function(char *, ...) REprintf;
    double function() unif_rand;
    double function() exp_rand;
    
    //Wrapping the initialization of an embedded R interpretig process
    int    function(int argc, char **argv) Rf_initEmbeddedR;
    void   function(int fatal) Rf_endEmbeddedR;
    
    SEXP     function(SEXP) Rf_protect;
    SEXP     function(SEXPTYPE, R_len_t) Rf_allocVector;
    SEXPREC* function(SEXPREC *x, SEXPREC *y)SETCAR;
    SEXPREC* function(char *)Rf_install;
    SEXP     function(SEXP, SEXP) Rf_eval;
    SEXP     function(SEXP, SEXP) Rf_findFun;
    SEXPREC* function(SEXPREC *, SEXPREC *, int *)R_tryEval;
    SEXPREC* function(SEXPREC *, SEXPREC *, int *)R_tryEvalSilent;
    
    void    function(int)Rf_unprotect;
    void    function(SEXPREC *)Rf_unprotect_ptr;
    
    int*    function(SEXPREC *x)LOGICAL;
    int*    function(SEXPREC *x)INTEGER;
    ubyte*  function(SEXPREC *x)RAW;
    double* function(SEXPREC *x)REAL;
  }
  
  SEXP NEW_LOGICAL(int n){return Rf_allocVector(LGLSXP,n); }
  SEXP NEW_INTEGER(int n){return Rf_allocVector(INTSXP,n);}
  SEXP NEW_NUMERIC(int n){return Rf_allocVector(REALSXP,n);}
  SEXP NEW_CHARACTER(int n){return  Rf_allocVector(STRSXP,n);}
  SEXP NEW_COMPLEX(int n){return Rf_allocVector(CPLXSXP,n);}
  SEXP NEW_LIST(int n){return Rf_allocVector(VECSXP,n);}
  SEXP NEW_STRING(int n){return NEW_CHARACTER(n);}
  SEXP NEW_RAW(int n){return Rf_allocVector(RAWSXP,n);}
  int* INTEGER_DATA(SEXPREC *x){return INTEGER(x);}
  SEXP PROTECT(SEXP e){ return Rf_protect(e); }
  void UNPROTECT(int e){ return Rf_unprotect(e); }
  
  char* function(int, size_t) R_alloc;
  char* function(long, size_t) S_alloc;
  char* function(char*, long, long, int) S_realloc;
  
  int function(SEXPREC *x)Rf_length;
  
  SEXPREC* function(SEXPREC *e)CDR;
  
  //Karl wants to bind:
  //norm_rand; set_seed
  //unif_rand; get_seed
  
  //set_seed ~= seed_in
  //get_seed ~= seed_out
  
  void LoadR(){
    HXModule lib = load_library("R");
    load_function(dnorm)(lib,"Rf_dnorm4");
    load_function(qf)(lib,"Rf_qf");

    load_function(R_SeedsSymbol)(lib,"R_SeedsSymbol");
    load_function(seed_in)(lib,"seed_in");
    load_function(seed_out)(lib,"seed_out");

    load_function(norm_rand)(lib,"norm_rand");
    load_function(Rprintf)(lib,"Rprintf");
    load_function(REprintf)(lib,"REprintf");
    load_function(unif_rand)(lib,"unif_rand");
    load_function(exp_rand)(lib,"exp_rand");
   
    load_function(Rf_initEmbeddedR)(lib,"Rf_initEmbeddedR");
    load_function(Rf_endEmbeddedR)(lib,"Rf_endEmbeddedR");
    
    load_function(Rf_protect)(lib,"Rf_protect");
    load_function(Rf_allocVector)(lib,"Rf_allocVector");
    load_function(SETCAR)(lib,"SETCAR");
    load_function(Rf_install)(lib,"Rf_install");
    load_function(Rf_findFun)(lib,"Rf_findFun");
    load_function(R_tryEval)(lib,"R_tryEval");
    load_function(Rf_eval)(lib,"Rf_eval");
    load_function(R_tryEvalSilent)(lib,"R_tryEvalSilent");
    
    load_function(Rf_unprotect)(lib,"Rf_unprotect");
    load_function(Rf_unprotect_ptr)(lib,"Rf_unprotect_ptr");
    
    load_function(LOGICAL)(lib,"LOGICAL");
    load_function(INTEGER)(lib,"INTEGER");
    load_function(RAW)(lib,"RAW");
    load_function(REAL)(lib,"REAL");
    load_function(CDR)(lib,"CDR");
    load_function(Rf_length)(lib,"Rf_length");
    
    
    load_function(R_alloc)(lib,"R_alloc");
    load_function(S_alloc)(lib,"S_alloc");
    load_function(S_realloc)(lib,"S_realloc");

    writeln("Loaded R functionality");
  }
  
}else{
  pragma(lib, "libR.so");
  
  extern(C){
    double dnorm(double, double, double, int);
    double qf(double, double, double, int, int);
    void R_SeedsSymbol(); // this is the tricky part, we need to combine calling these 3 functions to produce a new seed

    void seed_in(int*);
    void seed_out(int*);

    double norm_rand();
    double unif_rand();
    double exp_rand();

  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - norm_rand: " ~ to!string(norm_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand()));
  writeln("  - unif_rand: " ~ to!string(unif_rand())); 
}
