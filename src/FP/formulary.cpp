///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	5 December 2013
///////////////////////////////////////////////////////////

//  Standard libraries
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <string>
#include <cstring>
#include <math.h>
#include <map>

//  My libraries
#include "lib-array.h"

//  Declerations
#include "state.h"
#include "formulary.h"


//*******************************************************************
//-------------------------------------------------------------------
//  List all the units as in "label, conversion_factor" 
//-------------------------------------------------------------------
    Formulary::Formulary()  {
        // Physically: label * d = const in any system

        // H-Field
        D["W"]        = units("we", 1.0);
        D["W_si"]     = units("rad/s", Physical_Constants::gm_ratio);

        D["H"]        = units("we/g", 1.0);
        D["H_kOe"]    = units("kOe", 10.0);
        D["H_Oe"]     = units("Oe",  1.0e+4);
        D["H_si"]     = units("A/m", 4.0*M_PI*1.0e+7);

        D["B"]        =  units("Tesla", 1.0);
        D["B_si"]     =  units("Tesla", 1.0);
        D["B_G"]      =  units("Gauss", 1.0e+4);

        // time
        D["Time"]     = units("1/we", 1.0);    
        D["Time_cgs"] = units("sec",  1.0    /(Physical_Constants::gm_ratio) );    
        D["Time_si"]  = units("sec",  1.0    /(Physical_Constants::gm_ratio) );    
        D["Time_fs"]  = units("fsec", 1.0e+15/(Physical_Constants::gm_ratio) );    
        D["Time_ps"]  = units("psec", 1.0e+12/(Physical_Constants::gm_ratio) ); 
        D["Time_ns"]  = units("nsec", 1.0e+9 /(Physical_Constants::gm_ratio) ); 

        // Temperature
        D["T"]      = units("2we*mB/(g*kB)", 1.0 );    
        D["T_Tesla"]= units("Tesla",1.0 );    
        D["T_kOe"]  = units("kOe",  10.0 );
        /* Remove lande    
        D["T_J"]    = units("J",    Physical_Constants::g_lande * Physical_Constants::mB );    
        D["T_eV"]   = units("eV",   Physical_Constants::g_lande * Physical_Constants::mB
			          / Physical_Constants::e_Cb );    
        D["T_K"]    = units("K",    Physical_Constants::g_lande * Physical_Constants::mB
     	   		          / Physical_Constants::kB );    */
    }

    /*double Formulary::N_to_Kelvin(const double _Nor){
        return _Nor*D["T_K"].d; 
    }*/

    /*double Formulary::Kelvin_to_N(const double _kelvin){
        return _kelvin/D["T_K"].d; 
    }*/

    double Formulary::N_to_Sec(const double _Nor){
        return _Nor*D["Time_si"].d; 
    }

    double Formulary::Sec_to_N(const double _sec){
        return _sec/D["Time_si"].d; 
    }

    double Formulary::N_to_kOe(const double _Nor){
        return _Nor*D["H_kOe"].d; 
    }

    double Formulary::kOe_to_N(const double _kOe){
        return _kOe/D["H_kOe"].d; 
    }
//-------------------------------------------------------------------

    Formulary& formulary() {
        static Formulary f;
        return f;
    }
//-------------------------------------------------------------------
//*******************************************************************

//-------------------------------------------------------------------
//  Calculate the Magnetization for a given spin distribution
    MagneticField Eval_Magnetization(const SpinDistribution& SDin){

        double fouroverthree_pi(4.0*M_PI/3.0);
        MagneticField M;
        M.Hx() =  (-2.0) *fouroverthree_pi * SDin.real(1,1);
        M.Hy() = (2.0)*fouroverthree_pi * SDin.imag(1,1);
        M.Hz() =        fouroverthree_pi * SDin(1,0).real();

        return M;
    }
//-------------------------------------------------------------------
    
//-------------------------------------------------------------------
//  Calculate the Magnetization for a given spin distribution
    MagneticField_U1 Eval_Magnetization(const SpinDistribution_U1& SDin){

        double fouroverthree_pi(4.0*M_PI/3.0);
        MagneticField_U1 M;
        M.Hx() =  (-2.0) *fouroverthree_pi * SDin.real(1,1);
        M.Hy() = (2.0)*fouroverthree_pi * SDin.imag(1,1);
        M.Hz() =        fouroverthree_pi * SDin(1,0).real();

        return M;
    }
//-------------------------------------------------------------------
