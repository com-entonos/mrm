///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	12 March 2014
///////////////////////////////////////////////////////////


//  Contents:
//  ---------
//  1) RKFunctor/RK_Rotation member class: Rotation
//  2) RKFunctor member class: AdvectionE 
//  3) RKFunctor member class: AdvectionA
//  4) RKFunctor 
//  5) RK_Rotation
//  6) ImplicitDiffusion 
//  X) ImplicitRotation (only in z-direction, too inaccurate)


//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <algorithm>
    #include <cstdlib>

    #include <math.h>
    #include <map>

//  My libraries
    #include "lib-array.h"
    #include "lib-algorithms.h"

//  Declerations
    #include "state.h"
    #include "formulary.h"
    #include "operators.h"


//**************************************************************
//--------------------------------------------------------------


//  Calculate each of the coefficients required for applying the
//  operator for the rotation. Reference: notes
    OperatorSpace::Rotation::
    Rotation(const unsigned _Nl, const unsigned _Nm) 
        : A1(_Nm),
          A2(_Nl,_Nm),
          A3(-0.5),
          B1(_Nl) { 

//      Remember that A1 is meant to be imaginary
        for (unsigned m(0); m < _Nm; ++m) {
            A1[m] = (-1.0)*double(m); 
        }

//      A2 
        for (unsigned m(0); m < _Nm; ++m) {
            for (unsigned l(m); l < _Nl; ++l) {
                double el(l);
                double em(m);
                A2(l,m) = (0.5)*(em + el)*(el-em+1.0); 
            }
        }

//      B1
        for (unsigned l(0); l < _Nl; ++l) {
            double el(l);
            B1[l] = el*(el+1.0); 
        }
    }
//--------------------------------------------------------------

//  Apply each of the operators to every spherical harmonic
//  to find the corresponding slope. Follow the notes for the
//  effect of the rotation term
    void OperatorSpace::Rotation::
    operator()(const SpinDistribution& SDin, const MagneticField& MFin, //--> Input
                     SpinDistribution& SDh) {				//--> Output

        complex<double> Wp(MFin.Hy(),MFin.Hx());
        complex<double> Wm(MFin.Hy(),(-1.0)*MFin.Hx());
        complex<double> Wz(0.0,MFin.Hz());

//      A1: B(l,m), m > 0 
        for (unsigned m(1); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0()+1; ++l) {
                SDh.add_cplx( l,m,A1[m]*Wz*SDin(l,m) );
            }
        }    


//      B1: B(l,m-1), m = 1
        for (unsigned l(1); l < SDin.l0()+1; ++l) {
            SDh.real(l,0) += B1[l] * (   Wm.real()*SDin.real(l,1) 
                                       - Wm.imag()*SDin.imag(l,1) );
        }

//      A2: B(l,m-1), m > 1
        for (unsigned m(2); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0()+1; ++l) {
                SDh.add_cplx( l,m-1,A2(l,m)*Wm*SDin(l,m) );
            }
        }

//      A3: B(l,m+1), m < m0, l > m
        for (unsigned m(0); m < SDin.m0(); ++m) {
            for (unsigned l(m+1); l < SDin.l0()+1; ++l) {
                SDh.add_cplx( l,m+1,A3*Wp*SDin(l,m) );
            }
        }
        
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------

//  Calculate each of the coefficients required for applying the
//  operator that looks like an electric field. Reference: notes
    OperatorSpace::AdvectionE::
    AdvectionE(const unsigned _Nl, const unsigned _Nm) 
        : A1(_Nl, _Nm), A2(_Nl,_Nm),
          B1(_Nl),      B2(_Nl),
          C1(_Nl),      C3(_Nl),
          C2(_Nl, _Nm), C4(_Nl,_Nm) {

//      A1, A2 
        for (unsigned m(0); m < _Nm; ++m) {
            for (unsigned l(m); l < _Nl; ++l) {
                double el(l);
                double em(m);
                A1(l,m) = el*(el-em+1.0)/(2.0*el+1.0); 
                A2(l,m) = (-1.0)*(el+1.0)*(em+el)/(2.0*el+1.0); 

                C2(l,m) = 0.5* (el-em+2.0)*(el-em+1.0)* el   /(2.0*el+1.0);
                C4(l,m) = 0.5* (el+em-1.0)*(el+em)    *(el+1)/(2.0*el+1.0);
            }
        }

//      B1, B2, C1, C3
        for (unsigned l(0); l < _Nl; ++l) {
            double el(l);
            B1[l] = el*(el+1)* el   /(2.0*el+1.0); 
            B2[l] = el*(el+1)*(el+1)/(2.0*el+1.0);
            C1[l] = (-1.0)*0.5* el     /(2.0*el+1.0); 
            C3[l] = (-1.0)*0.5*(el+1.0)/(2.0*el+1.0);
        }

    }
//--------------------------------------------------------------

//  Apply each of the operators to every spherical harmonic
//  to find the corresponding slope. Follow the notes for the
//  effect of the term that looks like an effective electric
//  field. 
    void OperatorSpace::AdvectionE::
    operator()(const SpinDistribution& SDin, const MagneticField& MFin,  //--> Input
                     SpinDistribution& SDh,  const double _lambda) { 				 //--> Output                            

        complex<double> Wp(MFin.Hx(),MFin.Hy());         Wp *= _lambda;
        complex<double> Wm(MFin.Hx(),(-1.0)*MFin.Hy());  Wm *= _lambda;
        complex<double> Wz(MFin.Hz(),0);                 Wz *= _lambda;


//      A1: E(l+1,m),  l < l0
        for (unsigned m(0); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0(); ++l) {
                SDh.add_cplx( l+1,m,A1(l,m)*Wz*SDin(l,m) );
            }
        }    

//      A2:  E(l-1,m),  l > m
        for (unsigned m(0); m < SDin.m0()+1; ++m) {
            for (unsigned l(m+1); l < SDin.l0()+1; ++l) {
                SDh.add_cplx( l-1,m,A2(l,m)*Wz*SDin(l,m) );
            }
        }         

//      B1: E(l+1,0),  l < l0, m = 1
        for (unsigned l(1); l < SDin.l0(); ++l) {
            SDh.real(l+1,0) += B1[l]* (   Wp.real()*SDin.real(l,1) 
                                        - Wp.imag()*SDin.imag(l,1) ) ;
        }
//      B2: E(l-1,0),   m = 1
        for (unsigned l(1); l < SDin.l0()+1; ++l) {
            SDh.real(l-1,0) += B2[l]* (   Wp.real()*SDin.real(l,1) 
                                        - Wp.imag()*SDin.imag(l,1) ) ;
        }

//      C1: E(l+1,m+1),  m < m0, l < l0   
        for (unsigned m(0); m < SDin.m0(); ++m) {
            for (unsigned l(m); l < SDin.l0(); ++l) {
                SDh.add_cplx( l+1,m+1,C1[l]*Wm*SDin(l,m) );
            }
        }    

//      C2: E(l+1,m-1),  m > 1, l < l0
        for (unsigned m(2); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0(); ++l) {
                SDh.add_cplx( l+1,m-1,C2(l,m)*Wp*SDin(l,m) );
            }
        }    

//      C3: E(l-1,m+1),  m < m0, l > m+1
        for (unsigned m(0); m < SDin.m0(); ++m) {
            for (unsigned l(m+2); l < SDin.l0()+1; ++l) {
                SDh.add_cplx( l-1,m+1,C3[l]*Wm*SDin(l,m) );
            }
        }    

//      C4: E(l-1,m-1),  m > 1 
        for (unsigned m(2); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0()+1; ++l) {
                SDh.add_cplx( l-1,m-1,C4(l,m)*Wp*SDin(l,m) );
            }
        }

    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------

//  Calculate each of the coefficients required for applying the
//  operator that looks like spatial advection. Reference: notes
    OperatorSpace::AdvectionA::
    AdvectionA(const unsigned _Nl, const unsigned _Nm) 
        : A1(_Nl, _Nm), A2(_Nl,_Nm),
          B1(_Nl),      B2(_Nl),
          C1(_Nl),      C3(_Nl),
          C2(_Nl, _Nm), C4(_Nl,_Nm){

//      A1, A2,C2, C4 
        for (unsigned m(0); m < _Nm; ++m) {
            double em(m);
            for (unsigned l(m); l < _Nl; ++l) {
                double el(l);
                A1(l,m) = 2.0*(el-em+1.0)/(2.0*el+1.0); 
                A2(l,m) = 2.0*(em+el)/(2.0*el+1.0); 

                C2(l,m) =        (el-em+2.0)*(el-em+1.0)/(2.0*el+1.0); 
                C4(l,m) = (-1.0)*(em+el)    *(el+em-1.0)/(2.0*el+1.0); 
            }
        }

//      B1, B2, C1, C3
        for (unsigned l(0); l < _Nl; ++l) {
            double el(l);
            B1[l] =        2.0*el*(el+1)/(2.0*el+1.0); 
            B2[l] = (-1.0)*2.0*el*(el+1)/(2.0*el+1.0);
            C1[l] = (-1.0)/(2.0*el+1.0); 
            C3[l] =    1.0/(2.0*el+1.0); 
        }

    }
//--------------------------------------------------------------

//  Apply each of the operators to every spherical harmonic
//  to find the corresponding slope. Follow the notes for the
//  effect of the term that looks like spatial advection.
    void OperatorSpace::AdvectionA::
    operator()(const SpinDistribution& SDin, const MagneticField& MFin, //--> Input
                     SpinDistribution& SDh, const double _lambda) {	//--> Output

        complex<double> Wp(MFin.Hx(),MFin.Hy());           Wp *= _lambda;
        complex<double> Wm(MFin.Hx(),(-1.0)*MFin.Hy());    Wm *= _lambda;
        complex<double> Wz(MFin.Hz(),0.0);                 Wz *= _lambda;

//      A1: A(l+1,m),  l < l0
        for (unsigned m(0); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0(); ++l) {
                SDh.add_cplx( l+1,m,A1(l,m)*Wz*SDin(l,m) );
            }
        }    
//      A2: A(l-1,m),  l > m
        for (unsigned m(0); m < SDin.m0()+1; ++m) {
            for (unsigned l(m+1); l < SDin.l0()+1; ++l) {
                SDh.add_cplx( l-1,m,A2(l,m)*Wz*SDin(l,m) );
            }
        }         

//      B1: A(l+1,0),  l < l0, m = 1
        for (unsigned l(1); l < SDin.l0(); ++l) {
            SDh.real(l+1,0) +=  B1[l] * (  Wp.real()*SDin.real(l,1)
                                         - Wp.imag()*SDin.imag(l,1) );
        }
//      B2: A(l-1,0),  m = 1
        for (unsigned l(1); l < SDin.l0()+1; ++l) {
            SDh.real(l-1,0) +=  B2[l] * (  Wp.real()*SDin.real(l,1)
                                         - Wp.imag()*SDin.imag(l,1) );
        }

//      C1: A(l+1,m+1),  m < m0, l < l0
        for (unsigned m(0); m < SDin.m0(); ++m) {
            for (unsigned l(m); l < SDin.l0(); ++l) {
                SDh.add_cplx( l+1,m+1, C1[l]*Wm*SDin(l,m) );
            }
        }    
//      C2: A(l+1,m-1),  m > 1, l < l0
        for (unsigned m(2); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0(); ++l) {
                SDh.add_cplx( l+1,m-1, C2(l,m)*Wp*SDin(l,m) );
            }
        }    
//      C3: A(l-1,m+1),  m < m0, l > m+1
        for (unsigned m(0); m < SDin.m0(); ++m) {
            for (unsigned l(m+2); l < SDin.l0()+1; ++l) {
                SDh.add_cplx( l-1,m+1, C3[l]*Wm*SDin(l,m) );
            }
        }    
//      C4: A(l-1,m-1),  m > 1
        for (unsigned m(2); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0()+1; ++l) {
                SDh.add_cplx( l-1,m-1, C4(l,m)*Wp*SDin(l,m) );
            }
        }    


    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------

//  A generic Runge Kutta Functor initialized with the collection
//  of operators that have been implemented. In a Runge-Kutta method
//  solving df/dx = F(f,x) = F1(f,x)+F2(f,x)+... the RKfunctor is the
//  operator F, where the implementations of the specific F1, F2 etc
//  is described elsewhere.
    OperatorSpace::RKFunctor::
    RKFunctor(const unsigned _Nl, const unsigned _Nm) 
        :  //Hfld(_Nl,_Nm) , ---> subcycle for rotation
          Efld(_Nl,_Nm),
          Afld(_Nl,_Nm){
    }
//--------------------------------------------------------------

//  Collect all of the suboperators as in F = F1 + F2 ...
    void OperatorSpace::RKFunctor::
    operator()(const LocalState& Yin, LocalState& Yslope){

      //  MagneticField Htotal; //( Eval_Magnetization(Yin.SD()) ); //    Magnetization
        Htotal.Hx() = (-2.0) * Yin.SD().real(1,1); 
        Htotal.Hy() =   2.0  * Yin.SD().imag(1,1); 
        Htotal.Hz() =          Yin.SD().real(1,0); 
        Htotal *= (4.0*M_PI) * Yin.kB_Tc_over_mu;                    // => Exchange field
        Htotal += Yin.MF();			              // => Total field

        Yslope = 0.0;
        //Hfld( Yin.SD(), Htotal, Yslope.SD() );  ---> subcycle for rotation
        Efld( Yin.SD(), Htotal, Yslope.SD(), Yin.lambda );
        Afld( Yin.SD(), Htotal, Yslope.SD(), Yin.lambda );
    }
//--------------------------------------------------------------
//**************************************************************

//--------------------------------------------------------------

//  A Runge-Kutta Functor as above, but it only implements the 
//  rotation. It is used in case we need to subcycle for the
//  magnetic field`
    OperatorSpace::RK_Rotation::
    RK_Rotation(const unsigned _Nl, const unsigned _Nm) 
        : Hfld(_Nl,_Nm) {
    }
//--------------------------------------------------------------

//  Collect all of the suboperators as in F = F1 + F2 ...
    void OperatorSpace::RK_Rotation::
    operator()(const LocalState& Yin, LocalState& Yslope){

        //MagneticField Htotal; //( Eval_Magnetization(Yin.SD()) ); //    Magnetization
        Htotal.Hx() = (-2.0) * Yin.SD().real(1,1); 
        Htotal.Hy() =   2.0  * Yin.SD().imag(1,1); 
        Htotal.Hz() =          Yin.SD().real(1,0); 
        Htotal *= (4.0*M_PI)* Yin.kB_Tc_over_mu; 	            // => Exchange field
        Htotal += Yin.MF();			                    // => Total field

        Yslope = 0.0;
        Hfld( Yin.SD(), Htotal, Yslope.SD() );
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------

//  The diffusion term is initialized with the damping frequency
//  for each spherical harmonic.
    OperatorSpace::ImplicitDiffusion::
    ImplicitDiffusion(const unsigned _Nl) 
        : A1(_Nl){

//      A1 
        for (unsigned l(0); l < _Nl; ++l) {
            double el(l);
            A1[l] = (el+1.0)*el; 
        }

    }
//--------------------------------------------------------------

//  A Crank-Nicholson implementation of the diffusion 
    void OperatorSpace::ImplicitDiffusion::
    operator()(const double _T, const double _dt, SpinDistribution& SDin, const double _lambda){

        valarray<double> adv(A1); 
        adv *= _lambda * _T *_dt;
        for (unsigned l(0); l < adv.size(); ++l) {
            adv[l] = (1.0 - (1.0-theta)*adv[l])/ (1.0 + theta*adv[l]);
        }

        for (unsigned m(0); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0()+1; ++l) {
                SDin.mul_real( l,m,adv[l] );        
            }
        }    

    }
//--------------------------------------------------------------



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//  The diffusion term is initialized with the damping frequency
//  for each spherical harmonic.
    OperatorSpace::ImplicitRotation::
    ImplicitRotation(const unsigned _Nl, const unsigned _Nm)
        : Bl(_Nl),
          C0(_Nl),
          Bm(_Nm),
          Cm(_Nl, _Nm) {

//      Am 
        Am = 0.5;

//      Cm
        for (unsigned m(0); m < _Nm; ++m) {
            double em(m);
            Bm[m] = em;
            for (unsigned l(m); l < _Nl; ++l) {
                double el(l);
                Bl[l] = (el+1.0)*el; 
                Cm(l,m) = (-0.5)*(el+em+1)*(el-em); 
            }
        }

        for (unsigned l(0); l < _Nl; ++l) {
            double el(l);
            C0[l] = (-1.0)*el*(el+1.0);
        }

    }
//--------------------------------------------------------------

    void OperatorSpace::ImplicitRotation::
    operator()(const double _T, const double _dt, 
    	       const MagneticField& MFin, SpinDistribution& SDin,
 	       const double _lambda) {

        complex<double> Wp( MFin.Hy()*_dt, MFin.Hx()*_dt );
        complex<double> Wm( MFin.Hy()*_dt, (-1.0)*MFin.Hx()*_dt );
        complex<double> Wz( 0.0, MFin.Hz()*_dt );
        double Wt(_lambda*_T*_dt);

        Array2D< complex<double> > bm(Bl.size(), Bm.size());     
        for (unsigned m(0); m < SDin.m0()+1; ++m ){
            for (unsigned l(m); l < SDin.l0()+1; ++l) {
               bm(l,m) = 0.5*(/*Wt*Bl[l] +*/ Wz*Bm[m]);
                bm(l,m) = (1.0+bm(l,m)) / (1.0-bm(l,m)); 
             }
        } 
        
        for (unsigned m(1); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0()+1; ++l) {
                SDin.mul_cplx( l,m,bm(l,m) );           
             }
        }

    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
//--------------------------------------------------------------
    void OperatorSpaceU1::Precession2::
    operator() (const valarray<double>& SDin, const valarray<double>& MFin, //--> Input
                      valarray<double>& SDh){                               //--> Output
        // 0--> x, 1--> y, 2--> z
        //  (2,0) -> 0
        // R(2,1) -> 1
        // I(2,1) -> 2
        // R(2,2) -> 3
        // I(2,2) -> 4

        // (2,0)  +=  6.0 * ( Wx*I(2,1) + Wy*R(2,1) ); 
        // R(2,1) +=  Wz*I(2,1)- 0.5*Wy*R(2,0)+ 2.0*(Wy*R(2,2)+Wx*I(2,2));
        // I(2,1) += (-1.0)* Wz*R(2,1)- 0.5*Wx*R(2,0)+ 2.0*(Wy*I(2,2)-Wx*R(2,2));
        // R(2,2) +=   2.0 * Wz*I(2,2) - 0.5*(Wy*R(2,1)-Wx*I(2,1));
        // I(2,2) += (-2.0)* Wz*R(2,2) - 0.5*(Wx*R(2,1)+Wy*I(2,1));
        SDh[0] += 6.0 * ( MFin[0]*SDin[2] + MFin[1]*SDin[1] ); 
        SDh[1] += MFin[2]*SDin[2]- 0.5*MFin[1]*SDin[0]+ 2.0* (MFin[1]*SDin[3]+MFin[0]*SDin[4]);
        SDh[2] += (-1.0)* MFin[2]*SDin[1] - 0.5*MFin[0]*SDin[0]+ 2.0*(MFin[1]*SDin[4]-MFin[0]*SDin[3]);
        SDh[3] +=   2.0 * MFin[2]*SDin[4] - 0.5*(MFin[1]*SDin[1]-MFin[0]*SDin[2]);
        SDh[4] += (-2.0)* MFin[2]*SDin[3] - 0.5*(MFin[0]*SDin[1]+MFin[1]*SDin[2]);
    }

//--------------------------------------------------------------
    void OperatorSpaceU1::Damping::
    operator() (const LocalState_U1& Yin,          //--> Input
                      SpinDistribution_U1& SDh){   //--> Output

    double SH00(0.25/M_PI);

        // Exchange field
        double Wx((-8.0) * M_PI * Yin.kB_Tc_over_mu * (*Yin.SD().sd)[2]); 
        double Wy(  8.0  * M_PI * Yin.kB_Tc_over_mu * (*Yin.SD().sd)[5]); 
        double Wz(  4.0  * M_PI * Yin.kB_Tc_over_mu * (*Yin.SD().sd)[0]); 

        // Add external fields
        Wx += Yin.MF().Hx();    
        Wy += Yin.MF().Hy();    
        Wz += Yin.MF().Hz();    

        Wx *= Yin.lambda;
        Wy *= Yin.lambda;
        Wz *= Yin.lambda;

//      Omega z
        //(1,0) += Wz*[ 2.0*(0,0) - 0.4*(2,0)) ]
        //R(1,1) += (-0.6)*Wz * R(2,1);
        //I(1,1) += (-0.6)*Wz * I(2,1);
        (*SDh.sd)[0] +=   2.0  * Wz * ( SH00 - 0.2*  (*Yin.SD().sd)[1]  );
        (*SDh.sd)[2] += (-0.6) * Wz * (*Yin.SD().sd)[3];
        (*SDh.sd)[5] += (-0.6) * Wz * (*Yin.SD().sd)[6];

        //(2,0) += 2.0*Wz*(1,0);
        //R(2,1) += Wz * R(1,1);
        //I(2,1) += Wz * I(1,1);
        (*SDh.sd)[1] +=   2.0  * Wz * (*Yin.SD().sd)[0];
        (*SDh.sd)[3] +=          Wz * (*Yin.SD().sd)[2];
        (*SDh.sd)[6] +=          Wz * (*Yin.SD().sd)[5]; 

//      Omega x,y: l = 1
        // (1,0) += 1.2*[ Wx*R(2,1) - Wy*I(2,1) ]; 
        // R(1,1)+= 1.2*[ Wx*R(2,2) - Wy*I(2,2) ] - Wx*[ (0,0) + 0.1*(2,0) ]; 
        // I(1,1)+= 1.2*[ Wy*R(2,2) + Wx*I(2,2) ] + Wy*[ (0,0) + 0.1*(2,0) ]; 
        (*SDh.sd)[0] +=   1.2 * ( Wx* (*Yin.SD().sd)[3] - Wy* (*Yin.SD().sd)[6] ); 
        (*SDh.sd)[2] +=   1.2 * ( Wx* (*Yin.SD().sd)[4] - Wy* (*Yin.SD().sd)[7] ) - Wx*(SH00 + 0.1 * (*Yin.SD().sd)[1] ); 
        (*SDh.sd)[5] +=   1.2 * ( Wy* (*Yin.SD().sd)[4] + Wx* (*Yin.SD().sd)[7] ) + Wy*(SH00 + 0.1 * (*Yin.SD().sd)[1] );

//                 l = 2
        //(2,0)  +=   2.0 * ( Wx*R(1,1) - Wy*I(1,1) ); 
    //R(2,1) += (-0.5)*   Wx*(1,0);
        //I(2,1) +=   0.5 *   Wy*(1,0);
        //R(2,2) += (-0.5)* ( Wx*R(1,1) + Wy*I(1,1) );
        //I(2,2) +=   0.5 * ( Wy*R(1,1) - Wx*I(1,1) ); 
        (*SDh.sd)[1] +=   2.0 * ( Wx* (*Yin.SD().sd)[2] - Wy* (*Yin.SD().sd)[5] ); 
    (*SDh.sd)[3] += (-0.5)*Wx* (*Yin.SD().sd)[0];
        (*SDh.sd)[6] +=  0.5 * Wy* (*Yin.SD().sd)[0];
        (*SDh.sd)[4] += (-0.5)* ( Wx* (*Yin.SD().sd)[2] + Wy* (*Yin.SD().sd)[5] );
        (*SDh.sd)[7] +=   0.5 * ( Wy* (*Yin.SD().sd)[2] - Wx* (*Yin.SD().sd)[5] ); 

//      Rotation 
        // (1,0)  +=   2.0*( Wx*I(1,1) +     Wy*R(1,1) ); 
        // R(1,1) +=         Wz*I(1,1) - 0.5*Wy*R(1,0); 
        // I(1,1) += (-1.0)* Wz*R(1,1) - 0.5*Wx*R(1,0); 
        (*SDh.sd)[0] +=   2.0*( Yin.MF().Hx()* (*Yin.SD().sd)[5] +     Yin.MF().Hy()* (*Yin.SD().sd)[2] ); 
        (*SDh.sd)[2] +=         Yin.MF().Hz()* (*Yin.SD().sd)[5] - 0.5*Yin.MF().Hy()* (*Yin.SD().sd)[0]; 
        (*SDh.sd)[5] += (-1.0)* Yin.MF().Hz()* (*Yin.SD().sd)[2] - 0.5*Yin.MF().Hx()* (*Yin.SD().sd)[0]; 

    }

//--------------------------------------------------------------
    void OperatorSpaceU1::RK_Precession::
    ResetFields(const LocalState_U1& Yin){

        Htotal[0] = (-8.0) * M_PI * Yin.kB_Tc_over_mu * (*Yin.SD().sd)[2]; 
        Htotal[1] =   8.0  * M_PI * Yin.kB_Tc_over_mu * (*Yin.SD().sd)[5]; 
        Htotal[2] =   4.0  * M_PI * Yin.kB_Tc_over_mu * (*Yin.SD().sd)[0]; 

        // External fields
        Htotal[0]  += Yin.MF().Hx();    
        Htotal[1]  += Yin.MF().Hy();    
        Htotal[2]  += Yin.MF().Hz();    
    }

    void OperatorSpaceU1::RK_Precession::
    operator()(const valarray<double>& Yin, valarray<double>& Yslope){

        Yslope = 0.0;
        Pfld( Yin, Htotal, Yslope );
    }

//  A generic Runge Kutta Functor initialized with the collection
//  of operators that have been implemented. In a Runge-Kutta method
//  solving df/dx = F(f,x) = F1(f,x)+F2(f,x)+... the RKfunctor is the
//  operator F, where the implementations of the specific F1, F2 etc
//  is described elsewhere.
    OperatorSpaceU1::RKFunctor::
    RKFunctor(const unsigned _Nl, const unsigned _Nm) 
        //:  //Hfld(_Nl,_Nm) , ---> subcycle for rotation
        //  Efld(_Nl,_Nm),
        //  Afld(_Nl,_Nm)
    { }
//--------------------------------------------------------------

//  Collect all of the suboperators as in F = F1 + F2 ...
    void OperatorSpaceU1::RKFunctor::
    operator()(const LocalState_U1& Yin, LocalState_U1& Yslope){

        Yslope = 0.0;
        Dfld( Yin, Yslope.SD() );
    }


//--------------------------------------------------------------
//  The diffusion term is initialized with the damping frequency
//  for each spherical harmonic.
    OperatorSpaceU1::ExactDiffusion::
    ExactDiffusion(const unsigned _Nl) 
        : A1(_Nl){
//      A1 
        for (unsigned l(0); l < _Nl; ++l) {
            double el(l);
            A1[l] = (-1.0)*(el+1.0)*el; 
        }
    }

//  A Crank-Nicholson implementation of the diffusion 
    void OperatorSpaceU1::ExactDiffusion::
    operator()(const double _T, const double _dt, SpinDistribution_U1& SDin, const double _lambda){

        valarray<double> adv(A1); 
        adv *= _lambda * _T *_dt;
        for (unsigned l(0); l < adv.size(); ++l) {
            adv[l] = exp(adv[l]); 
        }

        for (unsigned l(1); l < SDin.l0()+1; ++l) {
            SDin.real( l,0)*=adv[l];        
        }
        for (unsigned m(1); m < SDin.m0()+1; ++m) {
            for (unsigned l(m); l < SDin.l0()+1; ++l) {
                SDin.real( l,m)*=adv[l];        
                SDin.imag( l,m)*=adv[l];        
            }
        }    
    }

    void OperatorSpaceU1::ExactDiffusion::
    operator()(const double _T, const double _dt, SpinDistribution_U1& SDin, const double _lambda, const unsigned fixed_l){
        double decay( exp( _lambda * _T *_dt * A1[fixed_l] ) );
        SDin.real(fixed_l,0)*=decay;        
        for (unsigned m(1); m < min(fixed_l,SDin.m0())+1; ++m) {
            SDin.real( fixed_l,m)*= decay;        
            SDin.imag( fixed_l,m)*= decay;        
        }
    }

