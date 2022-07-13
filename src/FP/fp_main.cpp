///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified: 6 January 2014
///////////////////////////////////////////////////////////

// Standard libraries 
#include <iostream>
#include <stdio.h>
#include <vector>
#include <valarray>
#include <complex>

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stdarg.h>
 
#include <map>
#include <string>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cstring>

// My libraries 
#include "lib-array.h"
#include "lib-algorithms.h"

// Misc Declerations
#include "state.h"
#include "formulary.h"
#include "operators.h"
#include "export.h"


//**************************************************************
// DECLERATIONS
//**************************************************************

//--------------------------------------------------------------
//  Default Parameters used for the Fokker-Planck method
    namespace DefaultParameters {
        const unsigned Nl = 5; 
        const unsigned Nm = 5; 
        const double   lambda = 0.1; 
        const unsigned Ncos8 = 120; 
        const unsigned Nphi = 120; 

        const double CFL(0.001);             // dt*lambda = 0.0001 for limiting the effect of 
                                            // numerical diffusion. Higher values may still be stable
        const double CFL_Rotation(0.0003);  // This is 1/20 of the period for 1000Tesla
                                            // Raising this above 0.0004 is not advisable
    }
//--------------------------------------------------------------

//  Decleration of a Wrapper Class for the Fokker-Planck code
    class WrapperFPE {
    public:
//      Constructors/Destructor
        WrapperFPE();
        WrapperFPE(const WrapperFPE& other);
        ~WrapperFPE(){ }

        void Advance(double& _Hx, double& _Hy, double& _Hz, 
                     double& _kB_T_over_mu, double& _kB_Tc_over_mu, 
                    double& _dt); 
        void Advance(double& _Hx, double& _Hy, double& _Hz, 
                     double& _kB_T_over_mu, double& _kB_Tc_over_mu, 
                    double& _dt, double& _lambda); 
        void Output(const size_t _fileID, const double _time); 
        void WriteRestart(const size_t _fileID); 
        void ReadRestart(const size_t _fileID); 

        unsigned Nl() const    { return __Nl; }
        unsigned Nm() const    { return __Nm; }
        double lambda() const  { return __lambda; }

        unsigned Ncos8() const { return __Ncos8; } 
        unsigned Nphi() const  { return __Nphi; } 

        void Variance(double& var);

    private:
//      Data 
        unsigned __Nl, __Nm;  // Number of harmonics in l, m directions
        LocalState Y;

//      Explicit Algorithm
        Algorithms::RK2<LocalState> RKadv; 
        Algorithms::RK4<LocalState> RKrot; 

//      Operators
        double __lambda;
        OperatorSpace::RKFunctor   rkF;
        OperatorSpace::RK_Rotation rkRotation;
        OperatorSpace::ImplicitDiffusion Diffusion;

//      Output
        unsigned __Ncos8, __Nphi;  // Number of cells for displaying the spin distribution
        Export_Files::Grid_Data gdata;
        Output_Data::Output_Preprocessor output; 

//      Restart facility
        Export_Files::Restart_Facility Re;
    };
//--------------------------------------------------------------

//  A factory for Fokker-Planck Wrapper class
    class FPsolver_Factory {
    public:
        int _number_FPsolver_created;
        int _number_FPsolver_destroyed;
        FPsolver_Factory() : _number_FPsolver_created(0), _number_FPsolver_destroyed(0) {};
        
        void destroy_a_FPsolver(WrapperFPE* fpe) {
            delete fpe;
            _number_FPsolver_destroyed++;
            
        }
     
        WrapperFPE* birth_a_FPsolver (){
            _number_FPsolver_created++;

            return new WrapperFPE();
        }
    };


//--------------------------------------------------------------


//**************************************************************
// DECLERATIONS
//**************************************************************

//--------------------------------------------------------------
//  Default Parameters used for the Fokker-Planck method
    namespace DefaultParameters_U1 {
        const unsigned Nl = 3; 
        const unsigned Nm = 3; 
        const double   lambda = 0.1; 
        const unsigned Ncos8 = 120; 
        const unsigned Nphi = 120; 

        const double CFL(0.003);             // dt*lambda = 0.0001 for limiting the effect of 
                                            // numerical diffusion. This is all about stability of diffusion
        const double CFL_Rotation(0.003);  // This is 1/20 of the period for 1000Tesla
                                            // Raising this above 0.0004 is not advisable
    }
//--------------------------------------------------------------

//  Decleration of a Wrapper Class for the Fokker-Planck code
    class WrapperFPE_U1 {
    public:
//      Constructors/Destructor
        WrapperFPE_U1();
        WrapperFPE_U1(const WrapperFPE_U1& other);
        ~WrapperFPE_U1(){ }

        void Advance(double& _Hx, double& _Hy, double& _Hz, 
                     double& _kB_T_over_mu, double& _kB_Tc_over_mu, 
                    double& _dt); 
        void Advance(double& _Hx, double& _Hy, double& _Hz, 
                     double& _kB_T_over_mu, double& _kB_Tc_over_mu, 
                    double& _dt, double& _lambda); 
        void Output(const size_t _fileID, const double _time); 
        void WriteRestart(const size_t _fileID); 
        void ReadRestart(const size_t _fileID); 

        unsigned Nl() const    { return __Nl; }
        unsigned Nm() const    { return __Nm; }
        double lambda() const  { return __lambda; }

        unsigned Ncos8() const { return __Ncos8; } 
        unsigned Nphi() const  { return __Nphi; } 

        void Variance(double& var);

    private:
//      Data 
        unsigned __Nl, __Nm;  // Number of harmonics in l, m directions
        LocalState_U1 Y;
        valarray<double> m2;

//      Explicit Algorithm
        Algorithms::RK4<LocalState_U1> RKadv; 
        //Algorithms::RK4<LocalState> RKrot; 
        Algorithms::RK4< valarray<double> > RKrot; 

//      Operators
        double __lambda;
        OperatorSpaceU1::RKFunctor      rkF;
        OperatorSpaceU1::RK_Precession  rkRotation;
        OperatorSpaceU1::ExactDiffusion Diffusion;

//      Output
        unsigned __Ncos8, __Nphi;  // Number of cells for displaying the spin distribution
        Export_Files_U1::Grid_Data gdata;
        Output_Data_U1::Output_Preprocessor output; 

//      Restart facility
        Export_Files_U1::Restart_Facility Re;
    };
//--------------------------------------------------------------

//  A factory for Fokker-Planck Wrapper class
    class FPsolver_Factory_U1 {
    public:
        int _number_FPsolver_created;
        int _number_FPsolver_destroyed;
        FPsolver_Factory_U1() : _number_FPsolver_created(0), _number_FPsolver_destroyed(0) {};
        
        void destroy_a_FPsolver(WrapperFPE_U1* fpe) {
            delete fpe;
            _number_FPsolver_destroyed++;
            
        }
     
        WrapperFPE_U1* birth_a_FPsolver (){
            _number_FPsolver_created++;

            return new WrapperFPE_U1();
        }
    };

//--------------------------------------------------------------



//**************************************************************
// DEFINITIONS
//**************************************************************

//**************************************************************
//  This Clock controls an iteration loop, given an initial
//  time tout_start*dt_out it evaluates the number of time-
//  steps it takes to get to (tout_start+1)*dt_out and can 
//  be incremented; 

    class Clock {
    public:
//      Constructor
        Clock(int tout_start, double dt_out, double CFL) {
            _hstep = 0;                        
            _t_start = double(tout_start)*dt_out;
            _numh  = size_t(static_cast<int>(dt_out/CFL))+1;
            _h     = dt_out/static_cast<double>(_numh);
        }

//      Clock readings
        double h()     const {return _h;}                    // Step size
        size_t numh()  const {return _numh;}                 // # steps 
        size_t tick()  const {return _hstep;}                // Current step
        double time()  const {return tick()*h() + _t_start;} // Current time
             
//      Increment time
        Clock& operator++() { ++_hstep; return *this;}  

    private:
        double _t_start;   // Initial output timestep
        size_t _numh;      // # of steps derived from this CFL
        double _h;         // resulting time-step from CFL
        size_t _hstep;
    };
//**************************************************************


//**************************************************************
//---------------------------------------------------------------------------
//  Constructor
    WrapperFPE::WrapperFPE()

//        FPE solver setup
        : __Nl(DefaultParameters::Nl), 
          __Nm(DefaultParameters::Nm), 
          __lambda(DefaultParameters::lambda),
          Y(__Nl,__Nm),
          RKadv(Y),
          RKrot(Y),
          rkF(__Nl,__Nm),
          rkRotation(__Nl,__Nm),
          Diffusion(__Nl),

//        output setup
          __Ncos8(DefaultParameters::Ncos8), 
          __Nphi(DefaultParameters::Nphi),
          gdata( __Nl, __Nm, __Ncos8, __Nphi),
          output(gdata){ 
                    
        Y.SD().assign(0,0,complex<double>(0.25/M_PI,0.0));  
    }

//  Copy constructor
    WrapperFPE::WrapperFPE(const WrapperFPE& other) 
        : __Nl(other.Nl() ), __Nm( other.Nm() ), __lambda( other.lambda() ),
          Y(__Nl,__Nm),
          RKadv(Y),
          RKrot(Y),
          rkF(__Nl,__Nm),
          rkRotation(__Nl,__Nm),
          Diffusion(__Nl),

//        output setup
          __Ncos8( other.Ncos8() ), __Nphi( other.Nphi() ),
          gdata( __Nl, __Nm, __Ncos8, __Nphi),
          output(gdata){ 
                    
        Y.SD().assign(0,0,complex<double>(0.25/M_PI,0.0));  
    }

//     Destructor is in the decleration, simply ~WrapperFPE(){ }
//---------------------------------------------------------------------------

//  Advance method: use default lambda
    void WrapperFPE::Advance(double& _Hx, double& _Hy, double& _Hz, 
                             double& _kB_T_over_mu, double& _kB_Tc_over_mu, 
                             double& _dt) {

        Y.MF().Hx() = _Hx;
        Y.MF().Hy() = _Hy;
        Y.MF().Hz() = _Hz;

        Y.kB_Tc_over_mu = _kB_Tc_over_mu; 
        Y.lambda        = lambda();         // Using default lambda

        double CFL(DefaultParameters::CFL);
        double CFL_Rotation(DefaultParameters::CFL_Rotation);  // This is ~3% of the period for 1000Tesla

//      Repeat until you get to the target step
        for (Clock W(0, _dt, CFL); W.tick() < W.numh(); ++W) { 
            Y = RKadv(Y, W.h(), &rkF);

//          Subcycling for rotation
            for (Clock Wr(0, W.h(), CFL_Rotation); Wr.tick() < Wr.numh(); ++Wr) { 
                Y = RKrot(Y, Wr.h(), &rkRotation); 
            }
            Diffusion(_kB_T_over_mu, W.h(), Y.SD(), Y.lambda);

        }

        MagneticField M( Eval_Magnetization(Y.SD()) ); 
        _Hx = M.Hx();
        _Hy = M.Hy();
        _Hz = M.Hz();
    }
//---------------------------------------------------------------------------

//  Advance method
    void WrapperFPE::Advance(double& _Hx, double& _Hy, double& _Hz, 
                             double& _kB_T_over_mu, double& _kB_Tc_over_mu, 
                             double& _dt, double& _lambda) {

        Y.MF().Hx() = _Hx;
        Y.MF().Hy() = _Hy;
        Y.MF().Hz() = _Hz;

        Y.kB_Tc_over_mu = _kB_Tc_over_mu; 
        Y.lambda        = _lambda; 

        double CFL( (DefaultParameters::CFL)*lambda() / Y.lambda );    // Rescale the CFL condition
        double CFL_Rotation(DefaultParameters::CFL_Rotation);          // This is ~3% of the period for 1000Tesla

//      Repeat until you get to the target step
        for (Clock W(0, _dt, CFL); W.tick() < W.numh(); ++W) { 
            Y = RKadv(Y, W.h(), &rkF);

//          Subcycling for rotation
            for (Clock Wr(0, W.h(), CFL_Rotation); Wr.tick() < Wr.numh(); ++Wr) { 
                Y = RKrot(Y, Wr.h(), &rkRotation); 
            }
            Diffusion(_kB_T_over_mu, W.h(), Y.SD(), Y.lambda);

        }

        MagneticField M( Eval_Magnetization(Y.SD()) ); 
        _Hx = M.Hx();
        _Hy = M.Hy();
        _Hz = M.Hz();

    }
//---------------------------------------------------------------------------

//  export the spin distribution
    void WrapperFPE::Output(const size_t _fileID, const double _time){
        output( Y, _fileID, _time);
    }

    void WrapperFPE::WriteRestart(const size_t _fileID){ 
        Re.Write(0, _fileID, Y);
    }

    void WrapperFPE::ReadRestart(const size_t _fileID){ 
        Re.Read(0, _fileID, Y);
    }
//---------------------------------------------------------------------------

//  Calculate the variance at the current state of the system
    void WrapperFPE::Variance(double& var) {
        double fouroverthree_pi(4.0*M_PI/3.0);
        double fourpi(4.0*M_PI);

        double Nx = (-2.0) * fouroverthree_pi * Y.SD().real(1,1);
        double Ny =   2.0  * fouroverthree_pi * Y.SD().imag(1,1);
        double Nz =          fouroverthree_pi * Y.SD().real(1,0);
        var += fourpi * Y.SD().real(0,0);
        var -= Nx*Nx + Ny*Ny + Nz*Nz;

        var +=  2.0*((-0.8)*fourpi*Y.SD().imag(2,2) - Nx*Ny);
        var +=  2.0*((-0.4)*fourpi*Y.SD().imag(2,1) - Nx*Nz);
        var +=  2.0*(  0.4 *fourpi*Y.SD().real(2,1) - Ny*Nz);
        
    }


//**************************************************************
//---------------------------------------------------------------------------
//  Constructor
    WrapperFPE_U1::WrapperFPE_U1()

//        FPE solver setup
        : __Nl(DefaultParameters_U1::Nl), 
          __Nm(DefaultParameters_U1::Nm), 
          __lambda(DefaultParameters_U1::lambda),
          Y(__Nl,__Nm),
          m2(5),
          RKadv(Y),
          RKrot(m2),
          rkF(__Nl,__Nm),
          //rkRotation(__Nl,__Nm),
          Diffusion(__Nl),

//        output setup
          __Ncos8(DefaultParameters_U1::Ncos8), 
          __Nphi(DefaultParameters_U1::Nphi),
          gdata( __Nl, __Nm, __Ncos8, __Nphi),
          output(gdata){ 
                    
       // Y.SD().assign(0,0,complex<double>(0.25/M_PI,0.0));  
    }

//  Copy constructor
    WrapperFPE_U1::WrapperFPE_U1(const WrapperFPE_U1& other) 
        : __Nl(other.Nl() ), __Nm( other.Nm() ), __lambda( other.lambda() ),
          Y(__Nl,__Nm),
          m2(5),
          RKadv(Y),
          RKrot(m2),
          rkF(__Nl,__Nm),
          //rkRotation(__Nl,__Nm),
          Diffusion(__Nl),

//        output setup
          __Ncos8( other.Ncos8() ), __Nphi( other.Nphi() ),
          gdata( __Nl, __Nm, __Ncos8, __Nphi),
          output(gdata){ 
                    
        //Y.SD().assign(0,0,complex<double>(0.25/M_PI,0.0));  
    }

//     Destructor is in the decleration, simply ~WrapperFPE(){ }
//---------------------------------------------------------------------------

//  Advance method: use default lambda
    void WrapperFPE_U1::Advance(double& _Hx, double& _Hy, double& _Hz, 
                                double& _kB_T_over_mu, double& _kB_Tc_over_mu, 
                                double& _dt) {

        Y.MF().Hx() = _Hx;
        Y.MF().Hy() = _Hy;
        Y.MF().Hz() = _Hz;

        Y.kB_Tc_over_mu = _kB_Tc_over_mu; 
        Y.lambda        = lambda();         // Using default lambda

        double CFL(DefaultParameters_U1::CFL);
        double CFL_Rotation(DefaultParameters_U1::CFL_Rotation);  // This is ~3% of the period for 1000Tesla

//      Repeat until you get to the target step
        for (Clock W(0, _dt, CFL); W.tick() < W.numh(); ++W) { 
            Y = RKadv(Y, W.h(), &rkF);
            Diffusion(_kB_T_over_mu, W.h(), Y.SD(), Y.lambda,1); // for l = 1

//          Subcycling for rotation
        //    m2[0] = (*Y.SD().sd)[1];
        //    m2[1] = (*Y.SD().sd)[3];
        //    m2[2] = (*Y.SD().sd)[6];
        //    m2[3] = (*Y.SD().sd)[4];
        //    m2[4] = (*Y.SD().sd)[7];
        //    rkRotation.ResetFields(Y);

        //    for (Clock Wr(0.0, W.h(), CFL_Rotation); Wr.tick() < Wr.numh(); ++Wr) { 
        //        m2 = RKrot(m2, Wr.h(), &rkRotation); 
        //    }
        //    (*Y.SD().sd)[1] = m2[0]; 
        //    (*Y.SD().sd)[3] = m2[1];
        //    (*Y.SD().sd)[6] = m2[2]; 
        //    (*Y.SD().sd)[4] = m2[3]; 
        //    (*Y.SD().sd)[7] = m2[4]; 

            Diffusion(_kB_T_over_mu, W.h(), Y.SD(), Y.lambda,2); // for l = 2

        }

        MagneticField_U1 M( Eval_Magnetization(Y.SD()) ); 
        _Hx = M.Hx();
        _Hy = M.Hy();
        _Hz = M.Hz();
    }
//---------------------------------------------------------------------------

//  Advance method
    void WrapperFPE_U1::Advance(double& _Hx, double& _Hy, double& _Hz, 
                                double& _kB_T_over_mu, double& _kB_Tc_over_mu, 
                                double& _dt, double& _lambda) {

        Y.MF().Hx() = _Hx;
        Y.MF().Hy() = _Hy;
        Y.MF().Hz() = _Hz;

        Y.kB_Tc_over_mu = _kB_Tc_over_mu; 
        Y.lambda        = _lambda; 

        double CFL( (DefaultParameters_U1::CFL)*lambda() / Y.lambda );    // Rescale the CFL condition
        double CFL_Rotation(DefaultParameters_U1::CFL_Rotation);          // This is ~3% of the period for 1000Tesla

//      Repeat until you get to the target step
        for (Clock W(0, _dt, CFL); W.tick() < W.numh(); ++W) { 
            Y = RKadv(Y, W.h(), &rkF);
            Diffusion(_kB_T_over_mu, W.h(), Y.SD(), Y.lambda,1); // for l = 1

//          Subcycling for rotation
        //    m2[0] = (*Y.SD().sd)[1];
        //    m2[1] = (*Y.SD().sd)[3];
        //    m2[2] = (*Y.SD().sd)[6];
        //    m2[3] = (*Y.SD().sd)[4];
        //    m2[4] = (*Y.SD().sd)[7];
        //    rkRotation.ResetFields(Y);

        //    for (Clock Wr(0.0, W.h(), CFL_Rotation); Wr.tick() < Wr.numh(); ++Wr) { 
        //        m2 = RKrot(m2, Wr.h(), &rkRotation); 
        //    }
        //    (*Y.SD().sd)[1] = m2[0]; 
        //    (*Y.SD().sd)[3] = m2[1];
        //    (*Y.SD().sd)[6] = m2[2]; 
        //    (*Y.SD().sd)[4] = m2[3]; 
        //    (*Y.SD().sd)[7] = m2[4]; 

            Diffusion(_kB_T_over_mu, W.h(), Y.SD(), Y.lambda,2); // for l = 2
        }

        MagneticField_U1 M( Eval_Magnetization(Y.SD()) ); 
        _Hx = M.Hx();
        _Hy = M.Hy();
        _Hz = M.Hz();

    }
//---------------------------------------------------------------------------

//  export the spin distribution
    //void WrapperFPE_U1::Output(const size_t _fileID, const double _time){
        // output( Y, _fileID, _time);
    //}

    void WrapperFPE_U1::WriteRestart(const size_t _fileID){ 
        Re.Write(0, _fileID, Y);
    }

    void WrapperFPE_U1::ReadRestart(const size_t _fileID){ 
        Re.Read(0, _fileID, Y);
    }
//---------------------------------------------------------------------------

//  Calculate the variance at the current state of the system
    void WrapperFPE_U1::Variance(double& var) {
        double fouroverthree_pi(4.0*M_PI/3.0);
        double fourpi(4.0*M_PI);

        double Nx = (-2.0) * fouroverthree_pi * Y.SD().real(1,1);
        double Ny =   2.0  * fouroverthree_pi * Y.SD().imag(1,1);
        double Nz =          fouroverthree_pi * Y.SD().real(1,0);
        var += 1.0;
        var -= Nx*Nx + Ny*Ny + Nz*Nz;

        var +=  2.0*((-0.8)*fourpi*Y.SD().imag(2,2) - Nx*Ny);
        var +=  2.0*((-0.4)*fourpi*Y.SD().imag(2,1) - Nx*Nz);
        var +=  2.0*(  0.4 *fourpi*Y.SD().real(2,1) - Ny*Nz);
        
    }



    //FPsolver_Factory __FPsolver_Factory;
    FPsolver_Factory_U1 __FPsolver_Factory;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
//  C-wrapper for the Fokker-Planck solvers
    extern "C" {

//      Kills a Fokker-Planck Wrapper
        //void CFP_delete(WrapperFPE* This) {
        //    __FPsolver_Factory.destroy_a_FPsolver(This);
        //}
        void CFP_delete(WrapperFPE_U1* This) {
            __FPsolver_Factory.destroy_a_FPsolver(This);
        }

//      Copies a Fokker-Planck Wrapper
        //WrapperFPE* CFP_copy(WrapperFPE This) {
        //    return new WrapperFPE(This);
        //}
        WrapperFPE_U1* CFP_copy(WrapperFPE_U1 This) {
            return new WrapperFPE_U1(This);
        }

//      Generates a new Fokker-Planck Wrapper
        //WrapperFPE* CFP_new() {
        //   return __FPsolver_Factory.birth_a_FPsolver();
        //}
        WrapperFPE_U1* CFP_new() {
           return __FPsolver_Factory.birth_a_FPsolver();
        }

//      Calls the "Advance" method of a given Fokker-Planck solver
//      with default lambda
        //void CFP_advance(WrapperFPE* This, double& Hx, double& Hy, double& Hz,
        //                 double& kB_T_over_mu, double& kB_Tc_over_mu, double& dt) {
        //    double dt_out(formulary().Sec_to_N(dt));
        //    This->Advance(Hx, Hy, Hz, kB_T_over_mu, kB_Tc_over_mu, dt_out);
        //   return;
        //}
        void CFP_advance(WrapperFPE_U1* This, double& Hx, double& Hy, double& Hz,
                         double& kB_T_over_mu, double& kB_Tc_over_mu, double& dt) {
            double dt_out(formulary().Sec_to_N(dt));
            This->Advance(Hx, Hy, Hz, kB_T_over_mu, kB_Tc_over_mu, dt_out);
            return;
        }

//      Calls the "Advance" method of a given Fokker-Planck solver
        //void CFP_advance_l(WrapperFPE* This, double& Hx, double& Hy, double& Hz, 
        //                   double& kB_T_over_mu, double& kB_Tc_over_mu, 
        //                   double& dt, double& lambda) {
        //    double dt_out(formulary().Sec_to_N(dt));
        //    This->Advance(Hx, Hy, Hz, kB_T_over_mu, kB_Tc_over_mu, dt_out, lambda);
            //double v(0.0);
            //This->Variance(v);
        //    return;
        //}
        void CFP_advance_l(WrapperFPE_U1* This, double& Hx, double& Hy, double& Hz, 
                           double& kB_T_over_mu, double& kB_Tc_over_mu, 
                           double& dt, double& lambda) {
            double dt_out(formulary().Sec_to_N(dt));
            This->Advance(Hx, Hy, Hz, kB_T_over_mu, kB_Tc_over_mu, dt_out, lambda);
            //double v(0.0);
            //This->Variance(v);
            return;
        }
    }
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

