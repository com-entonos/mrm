///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	3 January 2013
///////////////////////////////////////////////////////////

    #ifndef DECL_FORMULARY_H
    #define DECL_FORMULARY_H

//**************************************************************
//--------------------------------------------------------------
    namespace Physical_Constants {

//      Normalization density 
        const double Bnorm = 1.0;                   // normalized to 1 Tesla

//      Physical constants
        const double mu_0    = 4.0*M_PI*1.0e-7;
        // const double g_lande = 2.5;                   // For FePt     mu0 = 2.5muB
						           //     electron mu0 = 1.0muB
        const double exchJ   = 1.0;
        const double kB      = 1.380648813*1.0e-23;   // J/K
        const double mB      = 9.2740096820*1.0e-24;  // J/T
        const double hbar    = 1.054571726*1.0e-34;   // J*s

        const double c_si    = 2.99792458*1.0e+8;     // m/sec
        const double c       = 2.99792458*1.0e+10;    // cm/sec

        const double e_Cb    = 1.60217657*1.0e-19;    // Cb
        const double e       = 4.80320425*1.0e-10;    // Franklin or statC 

        const double m_kg    = 9.10938215*1.0e-31;    // kgr
        const double m       = 9.10938215*1.0e-28;    // gr

        const double gm_ratio = 1.7608597*1.0e+11;  // gyromagnetic ratio

//      Exchange field coefficient
        const double LambdaExchange = 1450.0;         // HE = LE * M

    }
//--------------------------------------------------------------



//--------------------------------------------------------------
//  Make Gaussian
    template<class T> 
    valarray<T> Gaussian(const size_t N, const T vmin, const T vmax, const T vth){

//      Make axis first
        valarray<T> G(N);
        for (size_t i(0); i < N; ++i) {
            G[i] = static_cast<T>(i);
        }
        G *= (vmax-vmin)/(static_cast<T>(N-1));
        G += vmin;

//      Make Gaussian
        T           C(pow( 1.0/ (sqrt(2.0*M_PI)*vth), 3));
        T           al( (-0.5) / (vth*vth));
        for (size_t i(0); i < N; ++i) {
            G[i]  = exp( al * G[i]*G[i] ); 
        }
        G *= C;

        return G;                              
    }


//--------------------------------------------------------------
    class units {
    public:
//      Contents
        string label;   // e.g. label = sec
        double d;       // e.g. x[label] = c * x[1/wp] => c = 1/wp 

//      Constructors
        units() : label("default"), d(0.0) {}
        units(string _x, double _d) : label(_x), d(_d){}

//      Copy constructor 
        units(const units& other) { 
           label = other.label; 
           d = other.d;
        } 
        ~units(){ }
    };

    class Formulary {
    public:
//    Construct an underlying dictionary for units
      Formulary();

//    Access to unit systems
      units Units(string key) { return D[key]; }
      units Units(string key1, string key2) { return D[key1+"_"+key2]; }
      string Label(string key) { return D[key].label; }
      string Label(string key1, string key2) { return D[key1+"_"+key2].label; }
      double Uconv(string key) { return D[key].d; }
      double Uconv(string key1, string key2) { return D[key1+"_"+key2].d; }

//    Data conversion
      double Sec_to_N(const double _sec);
      double N_to_Sec(const double _Nor);

//      double Kelvin_to_N(const double _kelvin);
//      double N_to_Kelvin(const double _Nor);

      double kOe_to_N(const double _kOe);
      double N_to_kOe(const double _Nor);
       
    private:
        map<string,units> D;

        static const double nmin = 1.0e-8;
    };
//--------------------------------------------------------------

    Formulary& formulary();
//--------------------------------------------------------------


//  Calculate the Magnetization for a given spin distribution
//  Decleration
    MagneticField Eval_Magnetization(const SpinDistribution& SDin);
    MagneticField_U1 Eval_Magnetization(const SpinDistribution_U1& SDin);
//**************************************************************


#endif
