///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	12 March 2014
///////////////////////////////////////////////////////////


    #ifndef DECL_VLASOVMAXWELL_H
    #define DECL_VLASOVMAXWELL_H


//********************************************************
//--------------------------------------------------------
namespace OperatorSpace {
//--------------------------------------------------------
//  Explicit operators used by the Runge-Kutta 
//--------------------------------------------------------
//  Implements df/dt - W*(N x df/dN) = 0
    class Rotation {
    public:
//      Constructors/Destructors
        Rotation(const unsigned _Nl, const unsigned _Nm);
        void operator()(const SpinDistribution& SDin, const MagneticField& MFin, 
	        	      SpinDistribution& SDh);
    private:
        valarray<double> A1;
        Array2D<double>  A2;
        double           A3;
        valarray<double> B1;
    };


//  Implements df/dt + lamda*W*df/dN = 0
    class AdvectionE {
    public:
//      Constructors/Destructors
        AdvectionE(const unsigned _Nl, const unsigned _Nm);
        void operator()(const SpinDistribution& SDin, const MagneticField& MFin, 
			      SpinDistribution& SDh,  const double _lambda);
    private:
        Array2D<double>   A1, A2;
        valarray<double>  B1, B2;
        valarray<double>  C1, C3;
        Array2D<double>   C2, C4;
    };

//  Implements df/dt - 2*lamda(N*W)f = 0
    class AdvectionA {
    public:
//      Constructors/Destructors
        AdvectionA(const unsigned _Nl, const unsigned _Nm);
        void operator()(const SpinDistribution& SDin, const MagneticField& MFin, 
		              SpinDistribution& SDh, const double _lambda);
    private:
        Array2D<double>   A1, A2;
        valarray<double>  B1, B2;
        valarray<double>  C1, C3;
        Array2D<double>   C2, C4;
    };

//--------------------------------------------------------
//  Wrapper Functors for the Explicit Operators 
//  Incorporate exchange field
//--------------------------------------------------------
//  Collects the Explicit operators
    class RKFunctor : public Algorithms::AbstFunctor<LocalState> {
    public:
//      Constructor
        RKFunctor(const unsigned _Nl, const unsigned _Nm); 
        ~RKFunctor(){ }

//      Collect all the operators and apply on Yin
        void operator()(const LocalState& Yin, LocalState& Yslope);

    private:
        // Rotation Hfld;
        AdvectionE Efld;
        AdvectionA Afld;

        // Exchange field
        MagneticField Htotal;
    };

//  Wraps only the Rotation operator
    class RK_Rotation: public Algorithms::AbstFunctor<LocalState> {
    public:
//      Constructor
        RK_Rotation(const unsigned _Nl, const unsigned _Nm);
        ~RK_Rotation(){ }

        void operator()(const LocalState& Yin, LocalState& Yslope);

    private:
        Rotation Hfld;

        // Exchange field
        MagneticField Htotal;
    };



//--------------------------------------------------------
//  IMPLICIT (Crank-Nickolson) OPERATORS:

//  Implements df/dt + (gamma*T*lamda/mu0)(l+1)*l*f = 0
    class ImplicitDiffusion {
    public:
//      Constructors/Destructors
	ImplicitDiffusion(const unsigned _Nl); 
        void operator()(const double _T, const double _dt, 
                        SpinDistribution& SDin, const double _lambda);

    private:
        valarray<double>  A1;

//      Choose theta for the implicit scheme
        static const double theta = 0.5;
    };



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
// IMPLICIT ROTATION PERFORMS VERY POORLY AND IS NOT USED
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

    class ImplicitRotation {
        public:
//          Constructors/Destructors
	    ImplicitRotation(const unsigned _Nl, const unsigned _Nm); 
            void operator()(const double _T, const double _dt, 
    			    const MagneticField& MFin,
                            SpinDistribution& SDin,
         		    const double _lambda);

        private:
            double    	      Am;
            valarray<double>  C0, Bl;
            valarray<double>  Bm;
            Array2D<double>   Cm;

//          For the suitable implicit scheme
            static const double theta = 0.5;
    };

//--------------------------------------------------------
}
//********************************************************

//********************************************************
//--------------------------------------------------------
namespace OperatorSpaceU1 {
//--------------------------------------------------------
//  Explicit operators used by the Runge-Kutta 
//--------------------------------------------------------

//  Implements df/dt + lamda*W*df/dN = 0
    class Damping {
    public:
//      Constructors/Destructors
        Damping(){}
        void operator()(const SpinDistribution_U1& SDin, const MagneticField_U1& MFin, 
                  SpinDistribution_U1& SDh,  const double _lambda);
        void operator()(const LocalState_U1& Yin, SpinDistribution_U1& SDh);
    private:
    };


//  Implements df/dt + lamda*W*df/dN = 0
//  for the second order m
    class Precession2 {
    public:
//      Constructors/Destructors
        Precession2(){}
        void operator()(const valarray<double>& SDin, const valarray<double>& MFin, 
                      valarray<double>& SDh);
    private:
    };

//--------------------------------------------------------
//  Wrapper Functors for the Explicit Operators 
//  Incorporate exchange field
//--------------------------------------------------------
//  Collects the Explicit operators
    class RKFunctor : public Algorithms::AbstFunctor<LocalState_U1> {
    public:
//      Constructor
        RKFunctor(const unsigned _Nl, const unsigned _Nm); 
        ~RKFunctor(){ }

//      Collect all the operators and apply on Yin
        void operator()(const LocalState_U1& Yin, LocalState_U1& Yslope);

    private:
        Damping Dfld;

        // Exchange field
        MagneticField_U1 Htotal;
    };

//  Wraps only the Rotation operator
    class RK_Precession: public Algorithms::AbstFunctor< valarray<double> > {
    private:
        Precession2 Pfld;

        // Exchange field
        valarray<double> Htotal;
    public:
//      Constructor
        RK_Precession():Htotal(3){}
        ~RK_Precession(){ }

        void ResetFields(const LocalState_U1& Yin);
        void operator()(const valarray<double>& vin, valarray<double>& vslope);

    };

//--------------------------------------------------------
//  Implements df/dt + (gamma*T*lamda/mu0)(l+1)*l*f = 0
    class ExactDiffusion {
    public:
//      Constructors/Destructors
    ExactDiffusion(const unsigned _Nl); 
    void operator()(const double _T, const double _dt, 
                    SpinDistribution_U1& SDin, const double _lambda);
    void operator()(const double _T, const double _dt, 
                    SpinDistribution_U1& SDin, const double _lambda, const unsigned fixed_l);

    private:
        valarray<double>  A1;
    };

//--------------------------------------------------------
}
//********************************************************


    #endif
