///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//  Last Modified: 12 March 2014 
///////////////////////////////////////////////////////////


    #ifndef DECL_STATE_H
    #define DECL_STATE_H

//--------------------------------------------------------------
//  Distribution function of a single spin in spherical harmonics. 
    class SpinDistribution {

    private:
        // Even memory locations for Real, even for Imaginary
        valarray<double>  *sd;    

        unsigned sz, Nl, Nm; 
        unsigned twoNl_1; 

    public:
//      Constructors/Destructors
        SpinDistribution(unsigned _Nl, unsigned _Nm);
        SpinDistribution(const SpinDistribution& other);
        ~SpinDistribution();

        valarray<double>& varray() const {return (*sd);}

//      Shape information
        unsigned m0()   const {return Nm-1;  }
        unsigned l0()   const {return Nl-1;  }
        unsigned dim()  const {return sz;}      // Total number of harmonics

//      Serial access -> can't reference complex numbers directly
        inline double& real(size_t i) 		{return (*sd)[2*i];}      
        inline double  real(size_t i) const 	{return (*sd)[2*i];}            
        inline double& imag(size_t i) 		{return (*sd)[2*i+1];}      
        inline double  imag(size_t i) const 	{return (*sd)[2*i+1];}            
//                    -> but can extract them
        inline complex<double> operator()(unsigned i) const { 
            return complex<double>((*sd)[2*i],(*sd)[2*i+1]);
        }

//      Random access -> can't reference complex numbers directly
        inline double& real(unsigned l, unsigned m) 
                        {return (*sd)[m*(twoNl_1-m)+2*l];} 
        inline double  real(unsigned l, unsigned m) const
                        {return (*sd)[m*(twoNl_1-m)+2*l];} 
        inline double& imag(unsigned l, unsigned m) 
                        {return (*sd)[m*(twoNl_1-m)+2*l+1];} 
        inline double  imag(unsigned l, unsigned m) const
                        {return (*sd)[m*(twoNl_1-m)+2*l+1];} 
//                    -> but can extract them
        inline complex<double> operator()(unsigned l, unsigned m) const { 
            unsigned i(m*(twoNl_1-m)+2*l);
            return complex<double>((*sd)[i],(*sd)[i+1]);
        }

//      Operators
//      ---------
//      Copy assignment
        SpinDistribution& assign(unsigned l, unsigned m, const complex<double>& c);
        SpinDistribution& operator=(const double& d);
        SpinDistribution& operator=(const complex<double>& c);
        SpinDistribution& operator=(const valarray<double>& vd);
        SpinDistribution& operator=(const valarray< complex<double> >& vc);
        SpinDistribution& operator=(const SpinDistribution& other);

//      Multiplication
        SpinDistribution& mul_real(unsigned l, unsigned m, const double& d);
        SpinDistribution& mul_imag(unsigned l, unsigned m, const double& d);
        SpinDistribution& mul_cplx(unsigned l, unsigned m, const complex<double>& c);
        SpinDistribution& operator*=(const double& d);
        SpinDistribution& operator*=(const complex<double>& c);
        SpinDistribution& operator*=(const valarray<double>& vd);
        SpinDistribution& operator*=(const valarray< complex<double> >& vc);
        SpinDistribution& operator*=(const SpinDistribution& sdmulti);

//      Addition
        SpinDistribution& add_cplx(unsigned l, unsigned m, const complex<double>& c);
        SpinDistribution& operator+=(const double& d);
        SpinDistribution& operator+=(const complex<double>& c);
        SpinDistribution& operator+=(const valarray<double>& vd);
        SpinDistribution& operator+=(const valarray< complex<double> >& vc);
        SpinDistribution& operator+=(const SpinDistribution& sdadd);

//      Subtraction
        SpinDistribution& sub_cplx(unsigned l, unsigned m, const complex<double>& c);
        SpinDistribution& operator-=(const double& d);
        SpinDistribution& operator-=(const complex<double>& c);
        SpinDistribution& operator-=(const valarray<double>& vd);
        SpinDistribution& operator-=(const valarray< complex<double> >& vc);
        SpinDistribution& operator-=(const SpinDistribution& sdadd);

    };
//--------------------------------------------------------------


//  Magnetic Field H decleration
    class MagneticField {

    private:
//      Just a real valarray with 3 variables (x,y,z) and interfaces. 
//      Convention: (0,1,2) --> (x,y,z)
        valarray<double>  *heta;    

    public:
//      Constructors/Destructors
        MagneticField(unsigned _dims = 3);
        MagneticField(const MagneticField& other);
        ~MagneticField();

        valarray<double>& varray() const {return (*heta);}

//      Shape information
        unsigned dim()  const {return (*heta).size();}     // valarray size 

//      Access
        inline double  operator()(unsigned i) const { return (*heta)[i]; }
        inline double& operator()(unsigned i)       { return (*heta)[i]; }

        inline double  Hx() const { return (*heta)[0]; } 
        inline double& Hx()       { return (*heta)[0]; } 
        inline double  Hy() const { return (*heta)[1]; } 
        inline double& Hy()       { return (*heta)[1]; } 
        inline double  Hz() const { return (*heta)[2]; } 
        inline double& Hz()       { return (*heta)[2]; } 

        complex<double> Hp(unsigned i, unsigned j) const;  // e.g. Hp(0,1) = Hx + i*Hy
        complex<double> Hm(unsigned i, unsigned j) const;  // e.g. Hm(0,1) = [Hp(0,1)]* = Hx - i*Hy

//      Operators
        MagneticField& operator=(const double& d);
        MagneticField& operator=(const valarray<double>& vd);
        MagneticField& operator=(const MagneticField& heq);

        MagneticField& operator*=(const double& d);
        MagneticField& operator*=(const valarray<double>& vd);
        MagneticField& operator*=(const MagneticField& hmul);

        MagneticField& operator+=(const double& d);
        MagneticField& operator+=(const valarray<double>& vd);
        MagneticField& operator+=(const MagneticField& hadd);

        MagneticField& operator-=(const double& d);
        MagneticField& operator-=(const valarray<double>& vd);
        MagneticField& operator-=(const MagneticField& hsub);

    };
//--------------------------------------------------------------

//  Collection of the local magnetic field and the spin distribution
    class LocalState {
    private:
        SpinDistribution *sp;
        MagneticField    *h;

    public:
//      Constructors/Destructors
        LocalState(unsigned _Nl, unsigned _Nm, unsigned _dims = 3);
        LocalState(const LocalState& other);            // Includes parameters
        ~LocalState();

        SpinDistribution& SD() const { return *sp; } 
        MagneticField& MF()    const { return *h; } 

//      Access
        SpinDistribution& SD()       { return *sp; } 
        MagneticField& MF()          { return *h; } 
    
//      Operators
        LocalState& operator=(const double& d);
        LocalState& operator=(const LocalState& other); // Includes parameters
        LocalState& operator*=(const double& d);
        LocalState& operator*=(const LocalState& smul);
        LocalState& operator+=(const double& d);
        LocalState& operator+=(const LocalState& sadd);
        LocalState& operator-=(const double& d);
        LocalState& operator-=(const LocalState& ssub);

//      Parameters
        double kB_Tc_over_mu;  // Included in copy construction and copy assignment
        double lambda;  
    };
//--------------------------------------------------------------
//**************************************************************

//--------------------------------------------------------------
//  Distribution function of a single spin in spherical harmonics. 
    class SpinDistribution_U1 {

    public:
        // First all the real data, then all imaginary data
        valarray<double>  *sd;    

        unsigned sz, Nl, Nm; 
        unsigned twoNl_1; 

//      Constructors/Destructors
        SpinDistribution_U1(unsigned _Nl, unsigned _Nm);
        SpinDistribution_U1(const SpinDistribution_U1& other);
        ~SpinDistribution_U1();

        valarray<double>& varray() const {return (*sd);}

//      Shape information
        unsigned m0()   const {return Nm-1;  }
        unsigned l0()   const {return Nl-1;  }
        unsigned dim()  const {return sz;}         // Total number of harmonics "0th" excluded
        unsigned size() const {return 2*sz-Nl+1;}  // Total number of elements

//      Serial access -> can't reference complex numbers directly
//      The usufulness of this is unclear
        inline double& real(size_t i)       {return (*sd)[i];}         //    -1 < i < sz    
        inline double  real(size_t i) const     {return (*sd)[i];}         //    -1 < i < sz  
        inline double& imag(size_t i)       {return (*sd)[i+sz-Nl+1];} //  Nl-2 < i < sz       
        inline double  imag(size_t i) const     {return (*sd)[i+sz-Nl+1];} //  Nl-2 < i < sz          
//                    -> but can extract them
        inline complex<double> operator()(unsigned i) const { 
            return complex<double>((*sd)[i],(*sd)[i+sz-Nl+1]); // Nl-2 < i < sz
        }

//      Random access -> can't reference complex numbers directly
//      This access may be expensive, use for probing 
        inline double& real(unsigned l, unsigned m) 
                        {return (*sd)[l-1+(2*Nl-m-1)*m/2];} 
        inline double  real(unsigned l, unsigned m) const
                        {return (*sd)[l-1+(2*Nl-m-1)*m/2];} 
        inline double& imag(unsigned l, unsigned m)             // 0 < m
                        {return (*sd)[sz+l-Nl+(2*Nl-m-1)*m/2];} 
        inline double  imag(unsigned l, unsigned m) const       // 0 < m
                        {return (*sd)[sz+l-Nl+(2*Nl-m-1)*m/2];}  
//                    -> but can extract them
        inline complex<double> operator()(unsigned l, unsigned m) const { // 0 < m
            unsigned i(l-1+(2*Nl-m-1)*m/2);       // e.g.: for (unsigned l(1); l < l0()+1; ++l)
            return complex<double>((*sd)[i],(*sd)[i+sz-Nl+1]); //      for (unsigned m(1); m < l; ++m)
        }

//      Operators
//      ---------
//      Copy assignment
/*void*/SpinDistribution_U1& assign(unsigned l, unsigned m, const complex<double>& c);
        SpinDistribution_U1& operator=(const double& d);
/*void*/SpinDistribution_U1& operator=(const complex<double>& c);
/*void*/SpinDistribution_U1& operator=(const valarray<double>& vd);
/*void*/SpinDistribution_U1& operator=(const valarray< complex<double> >& vc);
        SpinDistribution_U1& operator=(const SpinDistribution_U1& other);

//      Multiplication
/*void*/SpinDistribution_U1& mul_real(unsigned l, unsigned m, const double& d);
/*void*/SpinDistribution_U1& mul_imag(unsigned l, unsigned m, const double& d);
/*void*/SpinDistribution_U1& mul_cplx(unsigned l, unsigned m, const complex<double>& c);
        SpinDistribution_U1& operator*=(const double& d);
/*void*/SpinDistribution_U1& operator*=(const complex<double>& c);
/*void*/SpinDistribution_U1& operator*=(const valarray<double>& vd);
/*void*/SpinDistribution_U1& operator*=(const valarray< complex<double> >& vc);
        SpinDistribution_U1& operator*=(const SpinDistribution_U1& sdmulti);

//      Addition
/*void*/SpinDistribution_U1& add_cplx(unsigned l, unsigned m, const complex<double>& c);
        SpinDistribution_U1& operator+=(const double& d);
/*void*/SpinDistribution_U1& operator+=(const complex<double>& c);
/*void*/SpinDistribution_U1& operator+=(const valarray<double>& vd);
/*void*/SpinDistribution_U1& operator+=(const valarray< complex<double> >& vc);
        SpinDistribution_U1& operator+=(const SpinDistribution_U1& sdadd);

//      Subtraction
/*void*/SpinDistribution_U1& sub_cplx(unsigned l, unsigned m, const complex<double>& c);
        SpinDistribution_U1& operator-=(const double& d);
/*void*/SpinDistribution_U1& operator-=(const complex<double>& c);
/*void*/SpinDistribution_U1& operator-=(const valarray<double>& vd);
/*void*/SpinDistribution_U1& operator-=(const valarray< complex<double> >& vc);
        SpinDistribution_U1& operator-=(const SpinDistribution_U1& sdadd);

    };
//--------------------------------------------------------------

//  Magnetic Field H decleration
    class MagneticField_U1 {

    public:
//      Just a real valarray with 3 variables (x,y,z) and interfaces. 
//      Convention: (0,1,2) --> (x,y,z)
        valarray<double>  *heta;    

//      Constructors/Destructors
        MagneticField_U1(unsigned _dims = 3);
        MagneticField_U1(const MagneticField_U1& other);
        ~MagneticField_U1();

        valarray<double>& varray() const {return (*heta);}

//      Shape information
        unsigned dim()  const {return (*heta).size();}     // valarray size 

//      Access
        inline double  operator()(unsigned i) const { return (*heta)[i]; }
        inline double& operator()(unsigned i)       { return (*heta)[i]; }

        inline double  Hx() const { return (*heta)[0]; } 
        inline double& Hx()       { return (*heta)[0]; } 
        inline double  Hy() const { return (*heta)[1]; } 
        inline double& Hy()       { return (*heta)[1]; } 
        inline double  Hz() const { return (*heta)[2]; } 
        inline double& Hz()       { return (*heta)[2]; } 

        complex<double> Hp(unsigned i, unsigned j) const;  // e.g. Hp(0,1) = Hx + i*Hy
        complex<double> Hm(unsigned i, unsigned j) const;  // e.g. Hm(0,1) = [Hp(0,1)]* = Hx - i*Hy

//      Operators
        MagneticField_U1& operator=(const double& d);
        MagneticField_U1& operator=(const valarray<double>& vd);
        MagneticField_U1& operator=(const MagneticField_U1& heq);

        MagneticField_U1& operator*=(const double& d);
        MagneticField_U1& operator*=(const valarray<double>& vd);
        MagneticField_U1& operator*=(const MagneticField_U1& hmul);

        MagneticField_U1& operator+=(const double& d);
        MagneticField_U1& operator+=(const valarray<double>& vd);
        MagneticField_U1& operator+=(const MagneticField_U1& hadd);

        MagneticField_U1& operator-=(const double& d);
        MagneticField_U1& operator-=(const valarray<double>& vd);
        MagneticField_U1& operator-=(const MagneticField_U1& hsub);

    };
//--------------------------------------------------------------

//  Collection of the local magnetic field and the spin distribution
    class LocalState_U1 {
    private:
        SpinDistribution_U1 *sp;
        MagneticField_U1    *h;

    public:
//      Constructors/Destructors
        LocalState_U1(unsigned _Nl, unsigned _Nm, unsigned _dims = 3);
        LocalState_U1(const LocalState_U1& other);            // Includes parameters
        ~LocalState_U1();

        SpinDistribution_U1& SD() const { return *sp; } 
        MagneticField_U1& MF()    const { return *h; } 

//      Access
        SpinDistribution_U1& SD()       { return *sp; } 
        MagneticField_U1& MF()          { return *h; } 
    
//      Operators
        LocalState_U1& operator=(const double& d);
        LocalState_U1& operator=(const LocalState_U1& other); // Includes parameters
        LocalState_U1& operator*=(const double& d);
        LocalState_U1& operator*=(const LocalState_U1& smul);
        LocalState_U1& operator+=(const double& d);
        LocalState_U1& operator+=(const LocalState_U1& sadd);
        LocalState_U1& operator-=(const double& d);
        LocalState_U1& operator-=(const LocalState_U1& ssub);

//      Parameters
        double kB_Tc_over_mu;  // Included in copy construction and copy assignment
        double lambda;  
    };
//--------------------------------------------------------------
//**************************************************************
    #endif
