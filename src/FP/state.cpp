///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//  Last Modified: 12 March 2014
///////////////////////////////////////////////////////////

//  Contents:
//  ---------
//  1) LocalState member class: SpinDistribution
//  2) LocalState member class: MagneticField
//  3) LocalState
//  ---------
//  1) LocalState_U1 member class: SpinDistribution_U1
//  2) LocalState_U1 member class: MagneticField_U1
//  3) LocalState_U1


// Standard Libraries
   #include <iostream>
   #include <vector>
   #include <valarray>
   #include <complex>

// My Libraries
   #include "lib-array.h"

// Declerations
   #include "state.h"

//--------------------------------------------------------------
//**************************************************************

//  Constructor
    SpinDistribution::SpinDistribution(unsigned _Nl, unsigned _Nm)
        : Nl(_Nl),
          Nm(_Nm),
          twoNl_1(2*_Nl-1){

        sz = (Nm*(2*Nl-Nm+1))/2;
        sd = new valarray<double>(2*sz);
    }

//  Copy constructor
    SpinDistribution::SpinDistribution(const SpinDistribution& other){
        Nl = other.l0()+1;
        Nm = other.m0()+1;
        twoNl_1 = 2*Nl-1;
        sz = other.dim();
        sd = new valarray<double>(2*sz);
        *sd = other.varray();
    }

//  Destructor
    SpinDistribution:: ~SpinDistribution(){
//        int pauseandthink;
//        cout << "Killed distribution \n";
//        cin >> pauseandthink;
        delete sd; 
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
    SpinDistribution& 
    SpinDistribution::assign(unsigned l, unsigned m, const complex<double>& c){
        unsigned i(m*(twoNl_1-m)+2*l);
        (*sd)[i]   =c.real();
        (*sd)[i+1] =c.imag();
        return *this;
    } 

    SpinDistribution& SpinDistribution::operator=(const double& d){
        for (unsigned i(0); i< 2*dim(); i+=2) { 
            (*sd)[i]   = d;
            (*sd)[i+1] = 0.0;
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator=(const complex<double>& c){
        for (unsigned i(0); i< 2*dim(); i+=2) { 
            (*sd)[i]   =c.real();
            (*sd)[i+1] =c.imag();
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator=(const valarray<double>& vd){
        if ( dim() == vd.size()) {   
            for (unsigned i(0); i< dim(); ++i) { 
                (*sd)[2*i]   = vd[i];
                (*sd)[2*i+1] = 0.0;
            }
        }
        else  {
            cout << "In SpinDistribution '=' " << dim()<< "!= " << vd.size()<< ", terminating... \n";
            exit(1);
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator=(const valarray< complex<double> >& vc){
        if ( dim() == vc.size()) {  
            for (unsigned i(0); i< dim(); ++i) { 
                (*sd)[2*i]   = (vc[i]).real();
                (*sd)[2*i+1] = (vc[i]).imag();
            }
        }
        else  {
            cout << "In SpinDistribution '=' " << dim()<< "!= " << vc.size()<< ", terminating... \n";
            exit(1);
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator=(const SpinDistribution& other){
        if (this != &other) {   //protect from self-assignment
           *sd = other.varray(); 
        }
        else  {
            cout << "WARNING: self-assignment prevented\n"; 
        }
        return *this;
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  multiply with a real number
    SpinDistribution& SpinDistribution::mul_real(unsigned l, unsigned m, const double& d){
        unsigned i(m*(twoNl_1-m)+2*l);
        (*sd)[i]   *= d;
        (*sd)[i+1] *= d;
        return *this;
    } 

//  multiply with an imaginery number
    SpinDistribution& SpinDistribution::mul_imag(unsigned l, unsigned m, const double& d){
        unsigned i(m*(twoNl_1-m)+2*l);    // sd[i] + i*sd[i+1] = A+iB, 
        double tmp(d*(*sd)[i]);               // tmp = d*A 
        (*sd)[i]   = (-1.0)*d*(*sd)[i+1];     // sd[i] = -d*B
        (*sd)[i+1] = tmp;              // sd[i+1] = d*A
        return *this;
    }
    
//  multiply with a complex number
    SpinDistribution& 
    SpinDistribution::mul_cplx(unsigned l, unsigned m, const complex<double>& c){
        unsigned i(m*(twoNl_1-m)+2*l);  // sd[i] + i*sd[i+1] = A+iB, 
        double sdi((*sd)[i]*c.imag());            
        (*sd)[i]    *= c.real();                     // sd[i]   = Aa 
        (*sd)[i]    -= (*sd)[i+1]*c.imag();          // sd[i]   = Aa - Bb
        (*sd)[i+1]  *= c.real();                     // sd[i+1] = Ba
        (*sd)[i+1]  += sdi;                          // sd[i+1] = Ba + Ab
        return *this;
    }


//  *= 
    SpinDistribution& SpinDistribution::operator*=(const double& d){
        (*sd) *=d;
        return *this;
    }
    SpinDistribution& SpinDistribution::operator*=(const complex<double>& c){
        for (unsigned i(0); i< 2*dim(); i+=2) {         // suppose sd[i]+i*sd[i+1]=A+iB,   c=a+ib 
            double sdi((*sd)[i]*c.imag());            
            (*sd)[i]    *= c.real();                     // sd[i]   = Aa 
            (*sd)[i]    -= (*sd)[i+1]*c.imag();          // sd[i]   = Aa - Bb
            (*sd)[i+1]  *= c.real();                     // sd[i+1] = Ba
            (*sd)[i+1]  += sdi;                          // sd[i+1] = Ba + Ab
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator*=(const valarray<double>& vd){
        if ( dim() == vd.size()) {   
            for (unsigned i(0); i< dim(); ++i) {    
	        (*sd)[2*i]   *= vd[i];      
	        (*sd)[2*i+1] *= vd[i];     
	    }
        }
        else {
            cout << "In SpinDistribution '*=' " << dim()<< "!= " << vd.size()<< ", terminating... \n"; 
            exit(1);
        }
	return *this;
    }
    SpinDistribution& SpinDistribution::operator*=(const valarray< complex<double> >& vc){
        if ( dim() == vc.size()) {   
            for (unsigned i(0); i< dim(); ++i) {           // sd[i]+i*sd[i+1]=A+iB, sdmulti[i]+i*sdmulti[i+1]=a+ib 
                double sdi((*sd)[2*i] * vc[i].imag());            
	        (*sd)[2*i]    *=  vc[i].real();             // sd[i]   = Aa 
	        (*sd)[2*i]    -= (*sd)[2*i+1]*vc[i].imag(); // sd[i]   = Aa - Bb
	        (*sd)[2*i+1]  *= vc[i].real();              // sd[i+1] = Ba
	        (*sd)[2*i+1]  += sdi;                       // sd[i+1] = Ba + Ab
	    }
        }
        else {
            cout << "In SpinDistribution '*=' " << dim()<< "!= " << vc.size()<< ", terminating... \n"; 
            exit(1);
        }
	return *this;
    }
    SpinDistribution& SpinDistribution::operator*=(const SpinDistribution& sdmulti){
        for (unsigned i(0); i< 2*dim(); i+=2) {             // sd[i]+i*sd[i+1]=A+iB, sdmulti[i]+i*sdmulti[i+1]=a+ib 
            double sdi((*sd)[i] * sdmulti.varray()[i+1]);            
	    (*sd)[i]    *=  sdmulti.varray()[i];             // sd[i]   = Aa 
	    (*sd)[i]    -= (*sd)[i+1]*sdmulti.varray()[i+1]; // sd[i]   = Aa - Bb
	    (*sd)[i+1]  *= sdmulti.varray()[i];              // sd[i+1] = Ba
	    (*sd)[i+1]  += sdi;                              // sd[i+1] = Ba + Ab
	}
	return *this;
    }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  Add complex number at a sd(l,m)
    SpinDistribution& 
    SpinDistribution::add_cplx(unsigned l, unsigned m, const complex<double>& c){
        unsigned i(m*(twoNl_1-m)+2*l);
        (*sd)[i]   +=c.real();
        (*sd)[i+1] +=c.imag();
        return *this;
    } 

//  +=
    SpinDistribution& SpinDistribution::operator+=(const double& d){
        for (unsigned i(0); i< 2*dim(); i+=2) { 
            (*sd)[i]   += d;
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator+=(const complex<double>& c){
        for (unsigned i(0); i< 2*dim(); i+=2) { 
            (*sd)[i]   +=c.real();
            (*sd)[i+1] +=c.imag();
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator+=(const valarray<double>& vd){
        if ( dim() == vd.size()) {   
            for (unsigned i(0); i< dim(); ++i) { 
                (*sd)[2*i]   += vd[i];
            }
        }
        else {
            cout << "In SpinDistribution '+=' " << dim()<< "!= " << vd.size()<< ", terminating... \n"; 
            exit(1);
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator+=(const valarray< complex<double> >& vc){
        if ( dim() == vc.size()) {   
            for (unsigned i(0); i< dim(); ++i) { 
                (*sd)[2*i]   += (vc[i]).real();
                (*sd)[2*i+1] += (vc[i]).imag();
            }
        }
        else {
            cout << "In SpinDistribution '+=' " << dim()<< "!= " << vc.size()<< ", terminating... \n"; 
            exit(1);
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator+=(const SpinDistribution& other){
        *sd += other.varray(); 
        return *this;
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  Subtract complex number at a sd(l,m)
    SpinDistribution& 
    SpinDistribution::sub_cplx(unsigned l, unsigned m, const complex<double>& c){
        unsigned i(m*(twoNl_1-m)+2*l);
        (*sd)[i]   -=c.real();
        (*sd)[i+1] -=c.imag();
        return *this;
    } 

//  +=
    SpinDistribution& SpinDistribution::operator-=(const double& d){
        for (unsigned i(0); i< 2*dim(); i+=2) { 
            (*sd)[i]   -= d;
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator-=(const complex<double>& c){
        for (unsigned i(0); i< 2*dim(); i+=2) { 
            (*sd)[i]   -=c.real();
            (*sd)[i+1] -=c.imag();
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator-=(const valarray<double>& vd){
        if ( dim() == vd.size()) {   
            for (unsigned i(0); i< dim(); ++i) { 
                (*sd)[2*i]   -= vd[i];
            }
        }
        else {
            cout << "In SpinDistribution '-=' " << dim()<< "!= " << vd.size()<< ", terminating... \n"; 
            exit(1);
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator-=(const valarray< complex<double> >& vc){
        if ( dim() == vc.size()) {   
            for (unsigned i(0); i< dim(); ++i) { 
                (*sd)[2*i]   -= (vc[i]).real();
                (*sd)[2*i+1] -= (vc[i]).imag();
            }
        }
        else {
            cout << "In SpinDistribution '-=' " << dim()<< "!= " << vc.size()<< ", terminating... \n"; 
            exit(1);
        }
        return *this;
    }
    SpinDistribution& SpinDistribution::operator-=(const SpinDistribution& other){
        *sd -= other.varray(); 
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************

//  Constructor
    MagneticField::MagneticField(unsigned _dims) {
        heta = new valarray<double>(_dims);
    }

//  Copy constructor
    MagneticField::MagneticField(const MagneticField& other){
        heta = new valarray<double>(other.dim());
        *heta = other.varray();
    }

//  Destructor
    MagneticField:: ~MagneticField(){
        delete heta; 
    }


//--------------------------------------------------------------
//   Access
//--------------------------------------------------------------
     complex<double> MagneticField::Hp(unsigned i, unsigned j) const { 
        return complex<double>( (*heta)[i], (*heta)[j] );
     }
     complex<double> MagneticField::Hm(unsigned i, unsigned j) const { 
        return complex<double>( (*heta)[i],(-1.0)*(*heta)[j] );
     }


//--------------------------------------------------------------
//   Operators
//--------------------------------------------------------------
//   =
     MagneticField& MagneticField::operator=(const double& d){
        *heta = d;
        return *this;
     }
     MagneticField& MagneticField::operator=(const valarray<double>& vd){
        *heta = vd;
        return *this;
     }
    MagneticField& MagneticField::operator=(const MagneticField& heq){
        if (this != &heq) {   //self-assignment
            *heta = heq.varray();
        }
        return *this;
    }

//   +=
     MagneticField& MagneticField::operator*=(const double& d){
        *heta *= d;
        return *this;
     }
     MagneticField& MagneticField::operator*=(const valarray<double>& vd){
        *heta *= vd;
        return *this;
     }
     MagneticField& MagneticField::operator*=(const MagneticField& hmul){
        *heta *= hmul.varray();
        return *this;
     }

//   +=
     MagneticField& MagneticField::operator+=(const double& d){
        *heta += d;
        return *this;
     }
     MagneticField& MagneticField::operator+=(const valarray<double>& vd){
        *heta += vd;
        return *this;
     }
     MagneticField& MagneticField::operator+=(const MagneticField& hadd){
        *heta += hadd.varray();
        return *this;
     }

//   -=
     MagneticField& MagneticField::operator-=(const double& d){
        *heta -= d;
        return *this;
     }
     MagneticField& MagneticField::operator-=(const valarray<double>& vd){
        *heta -= vd;
        return *this;
     }
     MagneticField& MagneticField::operator-=(const MagneticField& hsub){
        *heta -= hsub.varray();
        return *this;
     }
//--------------------------------------------------------------
//**************************************************************


//  Constructor: Call individual constructors
    LocalState:: LocalState( unsigned _Nl, unsigned _Nm, unsigned _dims){
        if ( _Nm > _Nl ) { 
            cout << "m-index > l-index is not acceptable\n";
            exit(1);
        }
        sp = new SpinDistribution(_Nl,_Nm);
        h =  new MagneticField(_dims);

    }

//  Copy constructor: Call individual copy constructors
    LocalState:: LocalState(const LocalState& other) {
        sp = new SpinDistribution(other.SD());
        h =  new MagneticField(other.MF());
        
        kB_Tc_over_mu   = other.kB_Tc_over_mu;
        lambda         = other.lambda;
    }
//  Destructor
    LocalState:: ~LocalState(){
        "Destroying local state\n";
        delete sp;
        delete h;
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  =
    LocalState& LocalState::operator=(const double& d){
        *sp = d;
        *h  = d;
        return *this;
    }
//  =
    LocalState& LocalState::operator=(const LocalState& other){
        if (this != &other) {   //self-assignment
            *sp = other.SD();
            *h  = other.MF();
            kB_Tc_over_mu   = other.kB_Tc_over_mu;
            lambda         = other.lambda;
        }
        return *this;
    }

//  *=
    LocalState& LocalState::operator*=(const double& d){
        *sp *= d;
        *h  *= d;
        return *this;
    }
//  *=
    LocalState& LocalState::operator*=(const LocalState& smul){
        *sp *= smul.SD();
        *h  *= smul.MF();
        return *this;
    }

//  +=
    LocalState& LocalState::operator+=(const double& d){
        *sp += d;
        *h  += d;
        return *this;
    }
//  +=
    LocalState& LocalState::operator+=(const LocalState& sadd){
        *sp += sadd.SD();
        *h  += sadd.MF();
        return *this;
    }
//  +=
    LocalState& LocalState::operator-=(const double& d){
        *sp -= d;
        *h  -= d;
        return *this;
    }
//  +=
    LocalState& LocalState::operator-=(const LocalState& ssub){
        *sp -= ssub.SD();
        *h  -= ssub.MF();
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************
//**************************************************************

//  Constructor
    SpinDistribution_U1::SpinDistribution_U1(unsigned _Nl, unsigned _Nm)
        : Nl(_Nl),
          Nm(_Nm),
          twoNl_1(2*_Nl-1){

        sz = Nl-1+(2*Nl-Nm)*(Nm-1)/2;
        sd = new valarray<double>(2*sz-Nl+1);
    }

//  Copy constructor
    SpinDistribution_U1::SpinDistribution_U1(const SpinDistribution_U1& other){
        Nl = other.l0()+1;
        Nm = other.m0()+1;
        twoNl_1 = 2*Nl-1;
        sz = other.dim();
        sd = new valarray<double>(Nl-1+(2*Nl-Nm)*(Nm-1));
        *sd = other.varray();
    }

//  Destructor
    SpinDistribution_U1:: ~SpinDistribution_U1(){
        delete sd; 
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
    SpinDistribution_U1& 
    SpinDistribution_U1::assign(unsigned l, unsigned m, const complex<double>& c){
      //  unsigned i(l-1+(2*Nl-m-1)*m/2);
      //  (*sd)[i]   =c.real();
      //  (*sd)[i+sz-Nl+1] =c.imag();
        return *this;
    }

    SpinDistribution_U1& SpinDistribution_U1::operator=(const double& d){
        for (unsigned i(0);       i< dim();  ++i) { (*sd)[i] = d;  }
        for (unsigned i( dim() ); i< size(); ++i) { (*sd)[i] = 0.0;}

        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator=(const complex<double>& c){
      //  for (unsigned i(0);       i< dim();  ++i) { (*sd)[i] = c.real();  }
      //  for (unsigned i( dim() ); i< size(); ++i) { (*sd)[i] = c.imag();}

        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator=(const valarray<double>& vd){
      //  if ( dim() == vd.size()) {   
      //      for (unsigned i(0);       i< dim();  ++i) { (*sd)[i] = vd[i];  }
      //      for (unsigned i( dim() ); i< size(); ++i) { (*sd)[i] = 0.0;}
      //  }
      //  else  {
      //      cout << "In SpinDistribution_U1 '=' " << dim()<< "!= " << vd.size()<< ", terminating... \n";
      //      exit(1);
      //  }
        return *this;
    }
//  Complex is not interesting ---------------------------------------------------------
    SpinDistribution_U1& SpinDistribution_U1::operator=(const valarray< complex<double> >& vc){
      //  if ( dim() == vc.size()) {  
      //      for (unsigned i(0);       i< dim();  ++i) { (*sd)[i] = (vc[i]).real();  }
      //      for (unsigned i( dim() ); i< size(); ++i) { (*sd)[i] = (vc[i+Nl-1]).imag();  }
      //  }
      //  else  {
      //      cout << "In SpinDistribution_U1 '=' " << dim()<< "!= " << vc.size()<< ", terminating... \n";
      //      exit(1);
      //  }
        return *this;
    }//----------------------------------------------------------------------------------

    SpinDistribution_U1& SpinDistribution_U1::operator=(const SpinDistribution_U1& other){
        if (this != &other) {   //protect from self-assignment
           *sd = other.varray(); 
        }
        else  {
            cout << "WARNING: self-assignment prevented\n"; 
        }
        return *this;
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  multiply with a real number
    SpinDistribution_U1& SpinDistribution_U1::mul_real(unsigned l, unsigned m, const double& d){
      //  unsigned i(l-1+(2*Nl-m-1)*m/2);
      //  (*sd)[i]   *= d;
      //  (*sd)[i+sz-Nl+1] *= d;
        return *this;
    } 

//  multiply with an imaginery number
    SpinDistribution_U1& SpinDistribution_U1::mul_imag(unsigned l, unsigned m, const double& d){
      //  unsigned i(l-1+(2*Nl-m-1)*m/2);         // sd[i] + i*sd[i+sz-Nl+1] = A+iB, 
      //  double tmp(d*(*sd)[i]);                 // tmp = d*A 
      //  (*sd)[i]   = (-1.0)*d*(*sd)[i+sz-Nl+1]; // sd[i] = -d*B
      //  (*sd)[i+sz-Nl+1] = tmp;                 // sd[i+1] = d*A
        return *this;
    }
    
//  multiply with a complex number
    SpinDistribution_U1& 
    SpinDistribution_U1::mul_cplx(unsigned l, unsigned m, const complex<double>& c){
      //  unsigned i(l-1+(2*Nl-m-1)*m/2);              // sd[i] + i*sd[i+sz-Nl+1] = A+iB, 
      //  double sdi((*sd)[i]*c.imag());            
      //  (*sd)[i]    *= c.real();                     // sd[i]   = Aa 
      //  (*sd)[i]    -= (*sd)[i+sz-Nl+1]*c.imag();    // sd[i]   = Aa - Bb
      //  (*sd)[i+sz-Nl+1] *= c.real();                     // sd[i+1] = Ba
      //  (*sd)[i+sz-Nl+1]  += sdi;                          // sd[i+1] = Ba + Ab
        return *this;
    }


//  *= 
    SpinDistribution_U1& SpinDistribution_U1::operator*=(const double& d){
        (*sd) *=d;
        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator*=(const complex<double>& c){
      //  for (unsigned i(0); i< 2*dim(); i+=2) {         // suppose sd[i]+i*sd[i+1]=A+iB,   c=a+ib 
      //      double sdi((*sd)[i]*c.imag());            
      //      (*sd)[i]    *= c.real();                     // sd[i]   = Aa 
      //      (*sd)[i]    -= (*sd)[i+1]*c.imag();          // sd[i]   = Aa - Bb
      //      (*sd)[i+1]  *= c.real();                     // sd[i+1] = Ba
      //      (*sd)[i+1]  += sdi;                          // sd[i+1] = Ba + Ab
      //  }
        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator*=(const valarray<double>& vd){
      //  if ( dim() == vd.size()) {   
      //      for (unsigned i(0); i< dim(); ++i) {    
      //         (*sd)[2*i]   *= vd[i];      
      //        (*sd)[2*i+1] *= vd[i];     
      //    }
      //  }
      //  else {
      //      cout << "In SpinDistribution_U1 '*=' " << dim()<< "!= " << vd.size()<< ", terminating... \n"; 
      //      exit(1);
      //  }
    return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator*=(const valarray< complex<double> >& vc){
      //  if ( dim() == vc.size()) {   
      //      for (unsigned i(0); i< dim(); ++i) {           // sd[i]+i*sd[i+1]=A+iB, sdmulti[i]+i*sdmulti[i+1]=a+ib 
      //          double sdi((*sd)[2*i] * vc[i].imag());            
      //        (*sd)[2*i]    *=  vc[i].real();             // sd[i]   = Aa 
      //        (*sd)[2*i]    -= (*sd)[2*i+1]*vc[i].imag(); // sd[i]   = Aa - Bb
      //        (*sd)[2*i+1]  *= vc[i].real();              // sd[i+1] = Ba
      //        (*sd)[2*i+1]  += sdi;                       // sd[i+1] = Ba + Ab
      //    }
      // }
      //  else {
      //      cout << "In SpinDistribution_U1 '*=' " << dim()<< "!= " << vc.size()<< ", terminating... \n"; 
      //      exit(1);
      //  }
    return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator*=(const SpinDistribution_U1& sdmulti){
        for (unsigned i(0); i< Nl-1; ++i) {             // sd[i]+i*sd[i+1]=A+iB, sdmulti[i]+i*sdmulti[i+1]=a+ib 
           (*sd)[i]    *=  sdmulti.varray()[i];             // sd[i]   = Aa 
        }

        for (unsigned i(Nl-1); i < dim(); ++i) {               // sd[i]+i*sd[i+1]=A+iB, sdmulti[i]+i*sdmulti[i+1]=a+ib 
            double sdi((*sd)[i] * sdmulti.varray()[i+sz-Nl+1]);            
           (*sd)[i]          *=  sdmulti.varray()[i];                         // sd[i]         = Aa 
           (*sd)[i]          -= (*sd)[i+sz-Nl+1]*sdmulti.varray()[i+sz-Nl+1]; // sd[i]         = Aa - Bb
           (*sd)[i+sz-Nl+1]  *= sdmulti.varray()[i];                          // sd[i+sz-Nl+1] = Ba
           (*sd)[i+sz-Nl+1]  += sdi;                                          // sd[i+sz-Nl+1] = Ba + Ab
      }
    return *this;
    }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  Add complex number at a sd(l,m)
    SpinDistribution_U1& 
    SpinDistribution_U1::add_cplx(unsigned l, unsigned m, const complex<double>& c){
    //    unsigned i(m*(twoNl_1-m)+2*l);
    //    (*sd)[i]   +=c.real();
    //    (*sd)[i+1] +=c.imag();
        return *this;
    } 

//  +=
    SpinDistribution_U1& SpinDistribution_U1::operator+=(const double& d){
        for (unsigned i(0); i< dim(); ++i) { 
            (*sd)[i]   += d;
        }
        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator+=(const complex<double>& c){
    //    for (unsigned i(0); i< 2*dim(); i+=2) { 
    //        (*sd)[i]   +=c.real();
    //        (*sd)[i+1] +=c.imag();
    //    }
        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator+=(const valarray<double>& vd){
    //    if ( dim() == vd.size()) {   
    //        for (unsigned i(0); i< dim(); ++i) { 
    //            (*sd)[2*i]   += vd[i];
    //        }
    //    }
    //    else {
    //        cout << "In SpinDistribution_U1 '+=' " << dim()<< "!= " << vd.size()<< ", terminating... \n"; 
    //        exit(1);
    //    }
        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator+=(const valarray< complex<double> >& vc){
    //    if ( dim() == vc.size()) {   
    //        for (unsigned i(0); i< dim(); ++i) { 
    //            (*sd)[2*i]   += (vc[i]).real();
    //            (*sd)[2*i+1] += (vc[i]).imag();
    //        }
    //    }
    //    else {
    //        cout << "In SpinDistribution_U1 '+=' " << dim()<< "!= " << vc.size()<< ", terminating... \n"; 
    //        exit(1);
    //    }
        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator+=(const SpinDistribution_U1& other){
        *sd += other.varray(); 
        return *this;
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  Subtract complex number at a sd(l,m)
    SpinDistribution_U1& 
    SpinDistribution_U1::sub_cplx(unsigned l, unsigned m, const complex<double>& c){
    //    unsigned i(m*(twoNl_1-m)+2*l);
    //   (*sd)[i]   -=c.real();
    //    (*sd)[i+1] -=c.imag();
        return *this;
    } 

//  +=
    SpinDistribution_U1& SpinDistribution_U1::operator-=(const double& d){
        for (unsigned i(0); i< dim(); ++i) { 
            (*sd)[i]   -= d;
        }
        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator-=(const complex<double>& c){
    //    for (unsigned i(0); i< 2*dim(); i+=2) { 
    //        (*sd)[i]   -=c.real();
    //        (*sd)[i+1] -=c.imag();
    //    }
        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator-=(const valarray<double>& vd){
    //    if ( dim() == vd.size()) {   
    //        for (unsigned i(0); i< dim(); ++i) { 
    //            (*sd)[2*i]   -= vd[i];
    //        }
    //    }
    //    else {
    //        cout << "In SpinDistribution_U1 '-=' " << dim()<< "!= " << vd.size()<< ", terminating... \n"; 
    //        exit(1);
    //    }
        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator-=(const valarray< complex<double> >& vc){
    //    if ( dim() == vc.size()) {   
    //        for (unsigned i(0); i< dim(); ++i) { 
    //            (*sd)[2*i]   -= (vc[i]).real();
    //            (*sd)[2*i+1] -= (vc[i]).imag();
    //        }
    //    }
    //    else {
    //        cout << "In SpinDistribution_U1 '-=' " << dim()<< "!= " << vc.size()<< ", terminating... \n"; 
    //        exit(1);
    //    }
        return *this;
    }
    SpinDistribution_U1& SpinDistribution_U1::operator-=(const SpinDistribution_U1& other){
        *sd -= other.varray(); 
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************

//  Constructor
    MagneticField_U1::MagneticField_U1(unsigned _dims) {
        heta = new valarray<double>(_dims);
    }

//  Copy constructor
    MagneticField_U1::MagneticField_U1(const MagneticField_U1& other){
        heta = new valarray<double>(other.dim());
        *heta = other.varray();
    }

//  Destructor
    MagneticField_U1:: ~MagneticField_U1(){
        delete heta; 
    }


//--------------------------------------------------------------
//   Access
//--------------------------------------------------------------
     complex<double> MagneticField_U1::Hp(unsigned i, unsigned j) const { 
        return complex<double>( (*heta)[i], (*heta)[j] );
     }
     complex<double> MagneticField_U1::Hm(unsigned i, unsigned j) const { 
        return complex<double>( (*heta)[i],(-1.0)*(*heta)[j] );
     }


//--------------------------------------------------------------
//   Operators
//--------------------------------------------------------------
//   =
     MagneticField_U1& MagneticField_U1::operator=(const double& d){
        *heta = d;
        return *this;
     }
     MagneticField_U1& MagneticField_U1::operator=(const valarray<double>& vd){
        *heta = vd;
        return *this;
     }
    MagneticField_U1& MagneticField_U1::operator=(const MagneticField_U1& heq){
        if (this != &heq) {   //self-assignment
            *heta = heq.varray();
        }
        return *this;
    }

//   +=
     MagneticField_U1& MagneticField_U1::operator*=(const double& d){
        *heta *= d;
        return *this;
     }
     MagneticField_U1& MagneticField_U1::operator*=(const valarray<double>& vd){
        *heta *= vd;
        return *this;
     }
     MagneticField_U1& MagneticField_U1::operator*=(const MagneticField_U1& hmul){
        *heta *= hmul.varray();
        return *this;
     }

//   +=
     MagneticField_U1& MagneticField_U1::operator+=(const double& d){
        *heta += d;
        return *this;
     }
     MagneticField_U1& MagneticField_U1::operator+=(const valarray<double>& vd){
        *heta += vd;
        return *this;
     }
     MagneticField_U1& MagneticField_U1::operator+=(const MagneticField_U1& hadd){
        *heta += hadd.varray();
        return *this;
     }

//   -=
     MagneticField_U1& MagneticField_U1::operator-=(const double& d){
        *heta -= d;
        return *this;
     }
     MagneticField_U1& MagneticField_U1::operator-=(const valarray<double>& vd){
        *heta -= vd;
        return *this;
     }
     MagneticField_U1& MagneticField_U1::operator-=(const MagneticField_U1& hsub){
        *heta -= hsub.varray();
        return *this;
     }
//--------------------------------------------------------------
//**************************************************************


//  Constructor: Call individual constructors
    LocalState_U1:: LocalState_U1( unsigned _Nl, unsigned _Nm, unsigned _dims){
        if ( _Nm > _Nl ) { 
            cout << "m-index > l-index is not acceptable\n";
            exit(1);
        }
        sp = new SpinDistribution_U1(_Nl,_Nm);
        h =  new MagneticField_U1(_dims);

    }

//  Copy constructor: Call individual copy constructors
    LocalState_U1:: LocalState_U1(const LocalState_U1& other) {
        sp = new SpinDistribution_U1(other.SD());
        h =  new MagneticField_U1(other.MF());
        
        kB_Tc_over_mu   = other.kB_Tc_over_mu;
        lambda         = other.lambda;
    }
//  Destructor
    LocalState_U1:: ~LocalState_U1(){
        "Destroying local state\n";
        delete sp;
        delete h;
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  =
    LocalState_U1& LocalState_U1::operator=(const double& d){
        *sp = d;
        *h  = d;
        return *this;
    }
//  =
    LocalState_U1& LocalState_U1::operator=(const LocalState_U1& other){
        if (this != &other) {   //self-assignment
            *sp = other.SD();
            *h  = other.MF();
            kB_Tc_over_mu   = other.kB_Tc_over_mu;
            lambda         = other.lambda;
        }
        return *this;
    }

//  *=
    LocalState_U1& LocalState_U1::operator*=(const double& d){
        *sp *= d;
        *h  *= d;
        return *this;
    }
//  *=
    LocalState_U1& LocalState_U1::operator*=(const LocalState_U1& smul){
        *sp *= smul.SD();
        *h  *= smul.MF();
        return *this;
    }

//  +=
    LocalState_U1& LocalState_U1::operator+=(const double& d){
        *sp += d;
        *h  += d;
        return *this;
    }
//  +=
    LocalState_U1& LocalState_U1::operator+=(const LocalState_U1& sadd){
        *sp += sadd.SD();
        *h  += sadd.MF();
        return *this;
    }
//  +=
    LocalState_U1& LocalState_U1::operator-=(const double& d){
        *sp -= d;
        *h  -= d;
        return *this;
    }
//  +=
    LocalState_U1& LocalState_U1::operator-=(const LocalState_U1& ssub){
        *sp -= ssub.SD();
        *h  -= ssub.MF();
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************
