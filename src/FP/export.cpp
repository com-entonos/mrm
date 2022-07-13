///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	May 17 2013
///////////////////////////////////////////////////////////

//   
//   This cpp file contains the definitions for the functions
//   required to export the data
///////////////////////////////////////////////////////////
//
// 
//   class Export_Formatted_Data::
//
//   This class receives the output matrices and saves the 
//   data in txt files with recognizable name after it attaches
//   a small header with information necessary for plotting. 
//
// 
//   class Restart_Facility::
//
//   This class writes restart files from each node, and 
//   reads restart files for each node.  
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <algorithm>
    #include <fstream>
    #include <iomanip>
    #include <cstdlib>
    #include <sstream>
    #include <string>
    #include <cstring>

    #include <math.h>
    #include <map>

    #include <sys/stat.h>
    #include <sys/types.h>

//  My libraries
    #include "lib-array.h"
    #include "lib-algorithms.h"

//  Declerations
    #include "state.h"
    #include "formulary.h"
    #include "export.h"


//**************************************************************
//---------------------------------------------------------------
//   Create a folder
void Export_Files::Makefolder(string _name){

     mode_t _permissions(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
     char*   foldername = new char [_name.size()];
     strcpy(foldername, _name.c_str());
     int   status(mkdir(foldername,_permissions));

//     if (status != 0) cout << "Warning: Folder "<< _name<< " exists\n";

     delete[] foldername;
        
     return;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
// List of the acceptable Tags 
// Constructor 
Export_Files::DefaultTags::DefaultTags(size_t species){ 

//  Time 
    time.push_back( "Time_cgs");
    time.push_back( "Time_si" );
    time.push_back( "Time_fs" );
    time.push_back( "Time_ps" );
    time.push_back( "Time"    ); // Default tag

//  Space
    space.push_back( "Space_cgs");
    space.push_back( "Space_si" );
    space.push_back( "Space");

//  Fields
    momfld.push_back( "Ex"     ); 
    momfld.push_back( "Ex_cgs" );
    momfld.push_back( "Ex_si"  );
    momfld.push_back( "Ey"     );
    momfld.push_back( "Ey_cgs" );
    momfld.push_back( "Ey_si"  );
    momfld.push_back( "Ez"     );
    momfld.push_back( "Ez_cgs" );
    momfld.push_back( "Ez_si"  );
    momfld.push_back( "Bx"     );
    momfld.push_back( "Bx_cgs" );
    momfld.push_back( "Bx_si"  );
    momfld.push_back( "By"     );
    momfld.push_back( "By_cgs" );
    momfld.push_back( "By_si"  );
    momfld.push_back( "Bz"     );
    momfld.push_back( "Bz_cgs" );
    momfld.push_back( "Bz_si"  );
    momfld.push_back( "Jx"     );
    momfld.push_back( "Jx_cgs" );
    momfld.push_back( "Jx_si"  );
    momfld.push_back( "Jy"     );
    momfld.push_back( "Jy_cgs" );
    momfld.push_back( "Jy_si"  );
    momfld.push_back( "Jz"     );
    momfld.push_back( "Jz_cgs" );
    momfld.push_back( "Jz_si"  );

//  Moments
    momfld.push_back( "P"      );
    momfld.push_back( "P_cgs"  );
    momfld.push_back( "P_si"   );  
    momfld.push_back( "P_Mbar" );
    momfld.push_back( "T"      );
    momfld.push_back( "T_cgs"  );
    momfld.push_back( "T_si"   );  
    momfld.push_back( "T_eV" );
    momfld.push_back( "n"      );
    momfld.push_back( "n_cgs"  );
    momfld.push_back( "n_si"   );  

//  p-x
    for (size_t s(0); s < species; ++s) {
        pvsx.push_back( "px-x_"+stringify(s) ); 
        pvsx.push_back( "px-y_"+stringify(s) ); 
        pvsx.push_back( "px-z_"+stringify(s) ); 
        pvsx.push_back( "py-x_"+stringify(s) ); 
        pvsx.push_back( "py-y_"+stringify(s) ); 
        pvsx.push_back( "py-z_"+stringify(s) ); 
        pvsx.push_back( "pz-x_"+stringify(s) ); 
        pvsx.push_back( "pz-y_"+stringify(s) ); 
        pvsx.push_back( "pz-z_"+stringify(s) ); 
     }

}

// Convert data structure to float structure
valarray<float> Export_Files::vfloat(const valarray<double>& vDouble) {
    valarray<float> vf(vDouble.size());
    for (size_t i(0); i < vf.size(); ++i) {
        vf[i] = static_cast<float>(vDouble[i]);
    }
    return vf;
}
vector<float> Export_Files::vfloat(const vector<double> vDouble) {
    vector<float> vf;
    for (size_t i(0); i < vDouble.size(); ++i) {
        vf.push_back(static_cast<float>(vDouble[i]));
    }
    return vf;
}
//--------------------------------------------------------------

// Definition of the output axis
// Constructor 
Export_Files::oAxis::oAxis() : label(""), min(0.0), max(1.0), sz(3) {}
Export_Files::oAxis::oAxis( const float _m, const float _M, 
                            const size_t _sz) 
                : label(""), min(_m), max(_M), sz(_sz) {}
Export_Files::oAxis::oAxis(const string _l, const float _m, const float _M, 
                   const size_t _sz) 
                : label(_l), min(_m), max(_M), sz(_sz) {}
// Copy constructor 
Export_Files::oAxis::oAxis(const oAxis& other) { 
                label = other.label; 
                min   = other.min; 
                max   = other.max; 
                sz    = other.sz; 
} 
//--------------------------------------------------------------

//--------------------------------------------------------------
// 1D header constructor
Export_Files::Header::Header(oAxis _x, 
                             string _Ql, float _Qc, 
                             string _tl, float _tc,
                             string _oD)  
    : title(_Ql), titleC(_Qc), 
      time(_tl),  timeC(_tc),
      oDir(_oD) {
    xyz_axis.push_back(_x);
}

// 2D header constructor
Export_Files::Header::Header(oAxis _x, oAxis _y, 
                             string _Ql, float _Qc, 
                             string _tl, float _tc,
                             string _oD)  
    : title(_Ql),  time(_tl), 
      titleC(_Qc), timeC(_tc),
      oDir(_oD) {
    xyz_axis.push_back(_x);
    xyz_axis.push_back(_y);
}

// 3D header constructor
Export_Files::Header::Header(oAxis _x, oAxis _y, oAxis _z, 
                             string _Ql, float _Qc, 
                             string _tl, float _tc, 
                             string _oD)  
    : title(_Ql), time(_tl), 
      titleC(_Qc), timeC(_tc), 
      oDir(_oD) {
    xyz_axis.push_back(_x);
    xyz_axis.push_back(_y);
    xyz_axis.push_back(_z);
}

// xD header constructor
Export_Files::Header::Header(vector< oAxis > _xyz,
                             string _Ql, float _Qc, 
                             string _tl, float _tc, 
                             string _oD)  
    : xyz_axis(_xyz), title(_Ql), time(_tl), 
      titleC(_Qc), timeC(_tc), 
      oDir(_oD) {}

// number of header dimensions
size_t Export_Files::Header::dim() { 
    return xyz_axis.size(); 
}     
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

valarray<float> Export_Files::Header::axis(const size_t i) { 
    return Algorithms::MakeAxis(xyz_axis[i].min, xyz_axis[i].max, xyz_axis[i].sz); 
}
string Export_Files::Header::label(const size_t i) { 
    return xyz_axis[i].label; 
}
            
string Export_Files::Header::Title_label() { return title; }
float Export_Files::Header::Title_conv()  { return titleC; }
string Export_Files::Header::Time_label()  { return time; }
float Export_Files::Header::Time_conv()   { return timeC; }
string Export_Files::Header::Directory()   { return oDir; }
//--------------------------------------------------------------


//--------------------------------------------------------------
// Constructor of the export facility for data structures
Export_Files::Xport::Xport(const Grid_Data& _gdata, //const size_t Nphi, const size_t Ncos8,
                           //const vector< string > oTags,
                           string homedir){

    vector< string > oTags;
    oTags.push_back("Time_ps");
    oTags.push_back("Distribution");

    string folder(homedir + "OUTPUT/");
//    Makefolder(folder);
    Export_Files::oAxis phi("phi",0.0+M_PI/float(_gdata.Nphi),2*M_PI-M_PI/float(_gdata.Nphi),size_t(_gdata.Nphi));
    Export_Files::oAxis cos8("cos8",-1.0+1.0/float(_gdata.Ncos8),1.0-1.0/float(_gdata.Ncos8),size_t(_gdata.Ncos8));
    
    DefaultTags dTags(0);

/*    vector< oAxis > xyz, pxyz;    //TODO remove
    xyz.push_back( Export_Files::oAxis(_axis.xmin(0), _axis.xmax(0), _axis.Nx(0)) ); //TODO remove

    for (size_t s(0); s < species; ++s) {
        pxyz.push_back( Export_Files::oAxis( (-1.0)*float(_axis.pmax(s)), float(_axis.pmax(s)), _axis.Npx(s)) );
    }*/

//  Time 
    size_t tloc(0);                  // Find the location of the right tag
    while ( ( tloc < dTags.time.size()-1 ) &&
            ( find(oTags.begin(),oTags.end(), dTags.time[tloc]) == oTags.end() ) ) {
        ++tloc;
    } 
    string tlabel = "t[" +formulary().Label(dTags.time[tloc])+"]";
    float  tconv  =  formulary().Uconv(dTags.time[tloc]);

    Hdr["Distribution"] = Header(phi , cos8, 
	        "Probability Density",
	        1,
                tlabel, tconv, folder);

//  xyz Axis
/*    size_t xloc(0);                  // Find the location of the right tag
    while ( ( xloc < dTags.space.size()-1 ) &&
            ( find(oTags.begin(),oTags.end(), dTags.space[xloc]) == oTags.end() ) ) {
        ++xloc;
    } 
    xyz[0].label = "x["+ formulary().Label(dTags.space[xloc]) +"]";
    if ( xyz.size() > 1 ) xyz[1].label = "y["+ formulary().Label(dTags.space[xloc]) +"]";
    if ( xyz.size() > 2 ) xyz[2].label = "z["+ formulary().Label(dTags.space[xloc]) +"]";*/

//  pxyz Axis
/*    for (size_t i(0); i < species; ++i) {
        pxyz[i].label = "px[mc]";
        if ( pxyz.size()/species > 1 ) pxyz[  species+i].label = "py[mc]";
        if ( pxyz.size()/species > 2 ) pxyz[2*species+i].label = "pz[mc]";
    }

//  Tags for Moments and Fields -->
    for (size_t i(0); i < dTags.momfld.size(); ++i) {

 //     If this tag is an output tag
	if ( find(oTags.begin(),oTags.end(), dTags.momfld[i]) != oTags.end() ) { 

            string nounits = dTags.momfld[i].substr(0, dTags.momfld[i].find("_"));
            string folder = homedir + "OUTPUT/" + nounits + "/";
            Makefolder(folder);
 //         Generate a header file for this tag
            Hdr[dTags.momfld[i]] = Header(xyz, 
	        nounits+"["+formulary().Label(dTags.momfld[i])+"]",
	        formulary().Uconv(dTags.momfld[i]),
                tlabel, tconv, folder);    
        }
    } // <--

//  Tags for p-x -->
    for (size_t i(0); i < dTags.pvsx.size(); ++i) {

 //     If this tag is an output tag
	if ( find(oTags.begin(),oTags.end(), dTags.pvsx[i]) != oTags.end() ) { 

            string folder = homedir + "OUTPUT/" + dTags.pvsx[i] + "/";
            Makefolder(folder);
//          Generate a header file for this tag
//          For each 9 you have a different species 
            Hdr[dTags.pvsx[i]] = Header( pxyz[i/9+((i%9)/3)*species], xyz[i%3], 
                "f"+stringify(i/9), 1.0, tlabel, tconv, folder); 
        }
    }*/ //<--
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//export 1D data structure
void Export_Files::Xport::operator() 
          (const string tag,  valarray<float> data, 
           const size_t step, const double time,
           const size_t spec){

    string      filename(Hdr[tag].Directory());
    ofstream    oFile;

    Makefolder(Hdr[tag].Directory());
//  Check Header file correctness
    if (Hdr[tag].dim() != 1) {
        cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 1D structure\n"; 
        exit(1);
    }
 
//  Open File
    filename.append(tag).append(oFextension(spec,step));
    oFile.open(filename.c_str()); 
    oFile.flush();

//  Export dimensions
    oFile << Hdr[tag].dim() << "\n";

//  Export time
    oFile << Hdr[tag].Time_label() << "\n";
    //oFile << float(step)*Hdr[tag].Time_conv() << "\n";  // 
    oFile << float(time)*Hdr[tag].Time_conv() << "\n";  // 

//  Export all axes
    for (size_t i(0); i< Hdr[tag].dim(); ++i) {
        oFile << Hdr[tag].label(i) << "\n"; 
        oFile << Hdr[tag].axis(i);
    }
//  Renormalize and export data
    oFile << Hdr[tag].Title_label() << "\n"; 
    data *= Hdr[tag].Title_conv();         
    oFile << data; 

    oFile.flush();
    oFile.close();
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//export 2D data structure
void Export_Files::Xport::operator() 
          (const string tag,  Array2D<float> data, 
           const size_t step, const double time, 
           const size_t spec){

    string      filename(Hdr[tag].Directory());
    ofstream    oFile;
     
    Makefolder(Hdr[tag].Directory());
//  Check Header file correctness
    if (Hdr[tag].dim() != 2) {
        cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 2D structure\n"; 
        exit(1);
    }

//  Open File
    filename.append(tag).append(oFextension(spec,step));
    oFile.open(filename.c_str()); 
    oFile.flush();

//  Export dimensions
    oFile << Hdr[tag].dim() << "\n";

//  Export time
    oFile << Hdr[tag].Time_label() << "\n";
    //oFile << float(step)*Hdr[tag].Time_conv() << "\n";  // 
    oFile << float(time)*Hdr[tag].Time_conv() << "\n";  // 

//  Export all axes
    for (size_t i(0); i< Hdr[tag].dim(); ++i) {
        oFile << Hdr[tag].label(i) << "\n"; 
        oFile << Hdr[tag].axis(i);
    }
//  Renormalize and export data
    oFile << Hdr[tag].Title_label() << "\n"; 
    data *= Hdr[tag].Title_conv();         
    oFile << data; 

    oFile.flush();
    oFile.close();
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//export 3D data structure
void Export_Files::Xport::operator() 
          (const string tag,  Array3D<float> data, 
           const size_t step, const double time,
           const size_t spec){

    string      filename(Hdr[tag].Directory());
    ofstream    oFile;
     
    Makefolder(Hdr[tag].Directory());
//  Check Header file correctness
    if (Hdr[tag].dim() != 3) {
        cout << "ERROR "<<tag<<" : "  << Hdr[tag].dim() << " dimensions != 3D structure\n"; 
        exit(1);
    }

//  Open File
    filename.append(tag).append(oFextension(spec,step));
    oFile.open(filename.c_str()); 
    oFile.flush();

//  Export dimensions
    oFile << Hdr[tag].dim() << "\n";

//  Export time
    oFile << Hdr[tag].Time_label() << "\n";
    //oFile << float(step)*Hdr[tag].Time_conv() << "\n";  // 
    oFile << float(time)*Hdr[tag].Time_conv() << "\n";  // 

//  Export all axes
    for (size_t i(0); i< Hdr[tag].dim(); ++i) {
        oFile << Hdr[tag].label(i) << "\n"; 
        oFile << Hdr[tag].axis(i);
    }
//  Renormalize and export data
    oFile << Hdr[tag].Title_label() << "\n"; 
    data *= Hdr[tag].Title_conv();         
    oFile << data; 

    oFile.flush();
    oFile.close();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Adjust filenames with zeros to reach some prescribed length. 
//  Add the filename extension.  
string Export_Files::Xport::oFextension(size_t species, size_t step){
     stringstream sFilename;
     sFilename << "_s" << species <<"_";
        
     // Number of zeros to add to the filename
     int Nzeros(ofconventions::ofile_digits - stringify(step).length()); 
     while (Nzeros-- > 0) { 
         sFilename << "0"; 
     }
        
     sFilename << step << ofconventions::ofile_extension;

     return sFilename.str();
}
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
    Export_Files::Restart_Facility::Restart_Facility(string homedir) {
        hdir = homedir;
//        if (makefolder(hdir+"RESTART/") != 0) cout<<"Warning: Folder "<< hdir+"RESTART/"<<" exists\n";
    }

//--------------------------------------------------------------
//  Read restart file
    void Export_Files::Restart_Facility::Read(const int rank, const size_t re_step, LocalState& Y) {

//      Generate filename 
        string   filename(hdir+"RESTART/re_");
        filename.append(rFextension(rank,re_step));

//      Open file
        ifstream  fin(filename.c_str(), ios::binary);

//      Read spin distribution
	for (unsigned i(0); i < Y.SD().dim(); ++i) {
            fin.read( (char *)(&Y.SD().real(i)), sizeof(Y.SD().real(i)) );
            fin.read( (char *)(&Y.SD().imag(i)), sizeof(Y.SD().imag(i)) );
        }
            
//      Read H
        for(unsigned i(0); i < Y.MF().dim(); ++i){
            fin.read( (char *)(&Y.MF()(i)), sizeof(Y.MF()(i)) );
        }

        fin.close();
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Write restart file
    void Export_Files::Restart_Facility::Write(const int rank, const size_t re_step, LocalState& Y) {

//      Make restart folder if there is none
        Makefolder(hdir+"RESTART/");

//      Generate filename 
        string   filename(hdir+"RESTART/re_");
        filename.append(rFextension(rank,re_step));

//      Open file
        ofstream  fout(filename.c_str(), ios::binary);

//      Write spin distribution
	for (unsigned i(0); i < Y.SD().dim(); ++i) {
            fout.write( (char *)(&Y.SD().real(i)), sizeof(Y.SD().real(i)) );
            fout.write( (char *)(&Y.SD().imag(i)), sizeof(Y.SD().imag(i)) );
        }
            
//      Read H
        for(unsigned i(0); i < Y.MF().dim(); ++i){
            fout.write( (char *)(&Y.MF()(i)), sizeof(Y.MF()(i)) );
        }

        fout.flush();
        fout.close();
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Adjust filenames with zeros to reach some prescribed length. 
//  Add the filename extension. 
string Export_Files::Restart_Facility::rFextension(const int rank, const size_t rstep){
        stringstream sFilename;

        // Number of zeros to add to the filename
        int Nzeros(ofconventions::rank_digits - stringify(rank).length()); 
        while (Nzeros-- > 0) { 
            sFilename << "0"; 
        } 
        
        sFilename << rank << "_";

        // Number of zeros to add to the filename
        Nzeros = ofconventions::rfile_digits - stringify(rstep).length(); 
        while (Nzeros-- > 0) { 
            sFilename << "0";
        } 
        
        sFilename << rstep << ofconventions::rfile_extension;

        return sFilename.str();
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Output_Data::ShellDensity::ShellDensity(const Export_Files::Grid_Data& _G)
    : cos8( size_t(_G.Ncos8) ), phi( size_t(_G.Nphi) ),
      PLegendre( size_t(_G.Ncos8), size_t(_G.Nl), size_t(_G.Nm) ),
      COSmphi(   size_t(_G.Nphi), size_t(_G.Nm)  ),
      SINmphi(   size_t(_G.Nphi), size_t(_G.Nm)  ){

//  Make cos(8) axis
    float cos8min(-1.0+1.0/float(_G.Ncos8));
    float cos8max( 1.0-1.0/float(_G.Ncos8));
    cos8 = Algorithms::MakeAxis(cos8min, cos8max, size_t(_G.Ncos8));

//  Get legendre polynomials for each cosine
    for (size_t i(0); i < PLegendre.dim1(); ++i) {
        Array2D<float> aL( Algorithms::Legendre( cos8[i], size_t(_G.Nl), size_t(_G.Nm) ) );

//  Copy legendre polynomials into the array
        for (size_t l(0); l < PLegendre.dim2(); ++l) {
            for (size_t m(0); m < PLegendre.dim3(); ++m) {
                PLegendre(i,l,m) = aL(l,m); 
            }
        }
    }
    
//  Make phi axis
    float phimin( M_PI/float(_G.Nphi));
    float phimax( 2*M_PI-M_PI/float(_G.Nphi));
    phi =  Algorithms::MakeAxis( phimin,phimax, size_t(_G.Nphi)) ;

    for (size_t i(0); i < COSmphi.dim1(); ++i) {
        for (size_t m(0); m < COSmphi.dim2(); ++m) {
            COSmphi(i,m) = cos(m*phi[i]);
            SINmphi(i,m) = sin(m*phi[i]);
        }
    }
 
} 

float Output_Data::ShellDensity::COS8(const size_t _cos8) const { 
    return cos8[_cos8]; 
}
float Output_Data::ShellDensity::PHI(const size_t _phi) const { 
    return phi[_phi]; 
}
float Output_Data::ShellDensity::PL(const size_t _cos8, const size_t _l, const size_t _m) const { 
    return PLegendre( _cos8, _l, _m); 
}
float Output_Data::ShellDensity::cosMphi (const size_t _phi, const size_t _m) const {
    return COSmphi( _phi, _m);
}
float Output_Data::ShellDensity::sinMphi (const size_t _phi, const size_t _m) const {
    return SINmphi( _phi, _m);
}
//--------------------------------------------------------------
//**************************************************************

Output_Data::Output_Preprocessor::Output_Preprocessor(const Export_Files::Grid_Data& _gdata,  
                            //const vector< string > _oTags, 
                                                      string homedir)
        : expo( _gdata, //_oTags, 
                 homedir),
          sd( _gdata )
          //oTags(_oTags) 
          { }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void Output_Data::Output_Preprocessor::operator()(
    const LocalState& Y, const size_t tout, const double _time) { 

    float YSH_re(0.0), YSH_im(0.0), mphi(0.0);
    Array2D<float> phicos8( size_t(sd.Nphi()), size_t(sd.Ncos8()));

//  Loop for each location in phi-theta
//  Loop for each harmonic
    for (size_t i(0); i < sd.Nphi(); ++i) { 

        for (size_t l(0); l < sd.Nl(); ++l) {
            for (size_t m(0); m < sd.Nm(); ++m) {
                YSH_re = sd.cosMphi(i,m) * float(Y.SD().real(l,m));
                YSH_im = sd.sinMphi(i,m) * float(Y.SD().imag(l,m));
                YSH_re -= YSH_im; 

                for(size_t j(0); j < sd.Ncos8(); ++j) {
                    phicos8(i,j) += float(2.0) * sd.PL(j,l,m) * YSH_re;
                }

            }
        }

    }
    expo("Distribution", phicos8, tout, _time);
}
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//---------------------------------------------------------------
//   Create a folder
void Export_Files_U1::Makefolder(string _name){

     mode_t _permissions(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
     char*   foldername = new char [_name.size()];
     strcpy(foldername, _name.c_str());
     int   status(mkdir(foldername,_permissions));

//     if (status != 0) cout << "Warning: Folder "<< _name<< " exists\n";

     delete[] foldername;
        
     return;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
// List of the acceptable Tags 
// Constructor 
Export_Files_U1::DefaultTags::DefaultTags(size_t species){ 

//  Time 
    time.push_back( "Time_cgs");
    time.push_back( "Time_si" );
    time.push_back( "Time_fs" );
    time.push_back( "Time_ps" );
    time.push_back( "Time"    ); // Default tag

//  Space
    space.push_back( "Space_cgs");
    space.push_back( "Space_si" );
    space.push_back( "Space");

//  Fields
    momfld.push_back( "Ex"     ); 
    momfld.push_back( "Ex_cgs" );
    momfld.push_back( "Ex_si"  );
    momfld.push_back( "Ey"     );
    momfld.push_back( "Ey_cgs" );
    momfld.push_back( "Ey_si"  );
    momfld.push_back( "Ez"     );
    momfld.push_back( "Ez_cgs" );
    momfld.push_back( "Ez_si"  );
    momfld.push_back( "Bx"     );
    momfld.push_back( "Bx_cgs" );
    momfld.push_back( "Bx_si"  );
    momfld.push_back( "By"     );
    momfld.push_back( "By_cgs" );
    momfld.push_back( "By_si"  );
    momfld.push_back( "Bz"     );
    momfld.push_back( "Bz_cgs" );
    momfld.push_back( "Bz_si"  );
    momfld.push_back( "Jx"     );
    momfld.push_back( "Jx_cgs" );
    momfld.push_back( "Jx_si"  );
    momfld.push_back( "Jy"     );
    momfld.push_back( "Jy_cgs" );
    momfld.push_back( "Jy_si"  );
    momfld.push_back( "Jz"     );
    momfld.push_back( "Jz_cgs" );
    momfld.push_back( "Jz_si"  );

//  Moments
    momfld.push_back( "P"      );
    momfld.push_back( "P_cgs"  );
    momfld.push_back( "P_si"   );  
    momfld.push_back( "P_Mbar" );
    momfld.push_back( "T"      );
    momfld.push_back( "T_cgs"  );
    momfld.push_back( "T_si"   );  
    momfld.push_back( "T_eV" );
    momfld.push_back( "n"      );
    momfld.push_back( "n_cgs"  );
    momfld.push_back( "n_si"   );  

//  p-x
    for (size_t s(0); s < species; ++s) {
        pvsx.push_back( "px-x_"+stringify(s) ); 
        pvsx.push_back( "px-y_"+stringify(s) ); 
        pvsx.push_back( "px-z_"+stringify(s) ); 
        pvsx.push_back( "py-x_"+stringify(s) ); 
        pvsx.push_back( "py-y_"+stringify(s) ); 
        pvsx.push_back( "py-z_"+stringify(s) ); 
        pvsx.push_back( "pz-x_"+stringify(s) ); 
        pvsx.push_back( "pz-y_"+stringify(s) ); 
        pvsx.push_back( "pz-z_"+stringify(s) ); 
     }

}

// Convert data structure to float structure
valarray<float> Export_Files_U1::vfloat(const valarray<double>& vDouble) {
    valarray<float> vf(vDouble.size());
    for (size_t i(0); i < vf.size(); ++i) {
        vf[i] = static_cast<float>(vDouble[i]);
    }
    return vf;
}
vector<float> Export_Files_U1::vfloat(const vector<double> vDouble) {
    vector<float> vf;
    for (size_t i(0); i < vDouble.size(); ++i) {
        vf.push_back(static_cast<float>(vDouble[i]));
    }
    return vf;
}
//--------------------------------------------------------------

// Definition of the output axis
// Constructor 
Export_Files_U1::oAxis::oAxis() : label(""), min(0.0), max(1.0), sz(3) {}
Export_Files_U1::oAxis::oAxis( const float _m, const float _M, 
                            const size_t _sz) 
                : label(""), min(_m), max(_M), sz(_sz) {}
Export_Files_U1::oAxis::oAxis(const string _l, const float _m, const float _M, 
                   const size_t _sz) 
                : label(_l), min(_m), max(_M), sz(_sz) {}
// Copy constructor 
Export_Files_U1::oAxis::oAxis(const oAxis& other) { 
                label = other.label; 
                min   = other.min; 
                max   = other.max; 
                sz    = other.sz; 
} 
//--------------------------------------------------------------

//--------------------------------------------------------------
// 1D header constructor
Export_Files_U1::Header::Header(oAxis _x, 
                             string _Ql, float _Qc, 
                             string _tl, float _tc,
                             string _oD)  
    : title(_Ql), titleC(_Qc), 
      time(_tl),  timeC(_tc),
      oDir(_oD) {
    xyz_axis.push_back(_x);
}

// 2D header constructor
Export_Files_U1::Header::Header(oAxis _x, oAxis _y, 
                             string _Ql, float _Qc, 
                             string _tl, float _tc,
                             string _oD)  
    : title(_Ql),  time(_tl), 
      titleC(_Qc), timeC(_tc),
      oDir(_oD) {
    xyz_axis.push_back(_x);
    xyz_axis.push_back(_y);
}

// 3D header constructor
Export_Files_U1::Header::Header(oAxis _x, oAxis _y, oAxis _z, 
                             string _Ql, float _Qc, 
                             string _tl, float _tc, 
                             string _oD)  
    : title(_Ql), time(_tl), 
      titleC(_Qc), timeC(_tc), 
      oDir(_oD) {
    xyz_axis.push_back(_x);
    xyz_axis.push_back(_y);
    xyz_axis.push_back(_z);
}

// xD header constructor
Export_Files_U1::Header::Header(vector< oAxis > _xyz,
                             string _Ql, float _Qc, 
                             string _tl, float _tc, 
                             string _oD)  
    : xyz_axis(_xyz), title(_Ql), time(_tl), 
      titleC(_Qc), timeC(_tc), 
      oDir(_oD) {}

// number of header dimensions
size_t Export_Files_U1::Header::dim() { 
    return xyz_axis.size(); 
}     
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

valarray<float> Export_Files_U1::Header::axis(const size_t i) { 
    return Algorithms::MakeAxis(xyz_axis[i].min, xyz_axis[i].max, xyz_axis[i].sz); 
}
string Export_Files_U1::Header::label(const size_t i) { 
    return xyz_axis[i].label; 
}
            
string Export_Files_U1::Header::Title_label() { return title; }
float Export_Files_U1::Header::Title_conv()  { return titleC; }
string Export_Files_U1::Header::Time_label()  { return time; }
float Export_Files_U1::Header::Time_conv()   { return timeC; }
string Export_Files_U1::Header::Directory()   { return oDir; }
//--------------------------------------------------------------


//--------------------------------------------------------------
// Constructor of the export facility for data structures
Export_Files_U1::Xport::Xport(const Grid_Data& _gdata, //const size_t Nphi, const size_t Ncos8,
                           //const vector< string > oTags,
                           string homedir){

    vector< string > oTags;
    oTags.push_back("Time_ps");
    oTags.push_back("Distribution");

    string folder(homedir + "OUTPUT/");
//    Makefolder(folder);
    Export_Files_U1::oAxis phi("phi",0.0+M_PI/float(_gdata.Nphi),2*M_PI-M_PI/float(_gdata.Nphi),size_t(_gdata.Nphi));
    Export_Files_U1::oAxis cos8("cos8",-1.0+1.0/float(_gdata.Ncos8),1.0-1.0/float(_gdata.Ncos8),size_t(_gdata.Ncos8));
    
    DefaultTags dTags(0);

/*    vector< oAxis > xyz, pxyz;    //TODO remove
    xyz.push_back( Export_Files_U1::oAxis(_axis.xmin(0), _axis.xmax(0), _axis.Nx(0)) ); //TODO remove

    for (size_t s(0); s < species; ++s) {
        pxyz.push_back( Export_Files_U1::oAxis( (-1.0)*float(_axis.pmax(s)), float(_axis.pmax(s)), _axis.Npx(s)) );
    }*/

//  Time 
    size_t tloc(0);                  // Find the location of the right tag
    while ( ( tloc < dTags.time.size()-1 ) &&
            ( find(oTags.begin(),oTags.end(), dTags.time[tloc]) == oTags.end() ) ) {
        ++tloc;
    } 
    string tlabel = "t[" +formulary().Label(dTags.time[tloc])+"]";
    float  tconv  =  formulary().Uconv(dTags.time[tloc]);

    Hdr["Distribution"] = Header(phi , cos8, 
            "Probability Density",
            1,
                tlabel, tconv, folder);

//  xyz Axis
/*    size_t xloc(0);                  // Find the location of the right tag
    while ( ( xloc < dTags.space.size()-1 ) &&
            ( find(oTags.begin(),oTags.end(), dTags.space[xloc]) == oTags.end() ) ) {
        ++xloc;
    } 
    xyz[0].label = "x["+ formulary().Label(dTags.space[xloc]) +"]";
    if ( xyz.size() > 1 ) xyz[1].label = "y["+ formulary().Label(dTags.space[xloc]) +"]";
    if ( xyz.size() > 2 ) xyz[2].label = "z["+ formulary().Label(dTags.space[xloc]) +"]";*/

//  pxyz Axis
/*    for (size_t i(0); i < species; ++i) {
        pxyz[i].label = "px[mc]";
        if ( pxyz.size()/species > 1 ) pxyz[  species+i].label = "py[mc]";
        if ( pxyz.size()/species > 2 ) pxyz[2*species+i].label = "pz[mc]";
    }

//  Tags for Moments and Fields -->
    for (size_t i(0); i < dTags.momfld.size(); ++i) {

 //     If this tag is an output tag
    if ( find(oTags.begin(),oTags.end(), dTags.momfld[i]) != oTags.end() ) { 

            string nounits = dTags.momfld[i].substr(0, dTags.momfld[i].find("_"));
            string folder = homedir + "OUTPUT/" + nounits + "/";
            Makefolder(folder);
 //         Generate a header file for this tag
            Hdr[dTags.momfld[i]] = Header(xyz, 
            nounits+"["+formulary().Label(dTags.momfld[i])+"]",
            formulary().Uconv(dTags.momfld[i]),
                tlabel, tconv, folder);    
        }
    } // <--

//  Tags for p-x -->
    for (size_t i(0); i < dTags.pvsx.size(); ++i) {

 //     If this tag is an output tag
    if ( find(oTags.begin(),oTags.end(), dTags.pvsx[i]) != oTags.end() ) { 

            string folder = homedir + "OUTPUT/" + dTags.pvsx[i] + "/";
            Makefolder(folder);
//          Generate a header file for this tag
//          For each 9 you have a different species 
            Hdr[dTags.pvsx[i]] = Header( pxyz[i/9+((i%9)/3)*species], xyz[i%3], 
                "f"+stringify(i/9), 1.0, tlabel, tconv, folder); 
        }
    }*/ //<--
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//export 1D data structure
void Export_Files_U1::Xport::operator() 
          (const string tag,  valarray<float> data, 
           const size_t step, const double time,
           const size_t spec){

    string      filename(Hdr[tag].Directory());
    ofstream    oFile;

    Makefolder(Hdr[tag].Directory());
//  Check Header file correctness
    if (Hdr[tag].dim() != 1) {
        cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 1D structure\n"; 
        exit(1);
    }
 
//  Open File
    filename.append(tag).append(oFextension(spec,step));
    oFile.open(filename.c_str()); 
    oFile.flush();

//  Export dimensions
    oFile << Hdr[tag].dim() << "\n";

//  Export time
    oFile << Hdr[tag].Time_label() << "\n";
    //oFile << float(step)*Hdr[tag].Time_conv() << "\n";  // 
    oFile << float(time)*Hdr[tag].Time_conv() << "\n";  // 

//  Export all axes
    for (size_t i(0); i< Hdr[tag].dim(); ++i) {
        oFile << Hdr[tag].label(i) << "\n"; 
        oFile << Hdr[tag].axis(i);
    }
//  Renormalize and export data
    oFile << Hdr[tag].Title_label() << "\n"; 
    data *= Hdr[tag].Title_conv();         
    oFile << data; 

    oFile.flush();
    oFile.close();
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//export 2D data structure
void Export_Files_U1::Xport::operator() 
          (const string tag,  Array2D<float> data, 
           const size_t step, const double time, 
           const size_t spec){

    string      filename(Hdr[tag].Directory());
    ofstream    oFile;
     
    Makefolder(Hdr[tag].Directory());
//  Check Header file correctness
    if (Hdr[tag].dim() != 2) {
        cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 2D structure\n"; 
        exit(1);
    }

//  Open File
    filename.append(tag).append(oFextension(spec,step));
    oFile.open(filename.c_str()); 
    oFile.flush();

//  Export dimensions
    oFile << Hdr[tag].dim() << "\n";

//  Export time
    oFile << Hdr[tag].Time_label() << "\n";
    //oFile << float(step)*Hdr[tag].Time_conv() << "\n";  // 
    oFile << float(time)*Hdr[tag].Time_conv() << "\n";  // 

//  Export all axes
    for (size_t i(0); i< Hdr[tag].dim(); ++i) {
        oFile << Hdr[tag].label(i) << "\n"; 
        oFile << Hdr[tag].axis(i);
    }
//  Renormalize and export data
    oFile << Hdr[tag].Title_label() << "\n"; 
    data *= Hdr[tag].Title_conv();         
    oFile << data; 

    oFile.flush();
    oFile.close();
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//export 3D data structure
void Export_Files_U1::Xport::operator() 
          (const string tag,  Array3D<float> data, 
           const size_t step, const double time,
           const size_t spec){

    string      filename(Hdr[tag].Directory());
    ofstream    oFile;
     
    Makefolder(Hdr[tag].Directory());
//  Check Header file correctness
    if (Hdr[tag].dim() != 3) {
        cout << "ERROR "<<tag<<" : "  << Hdr[tag].dim() << " dimensions != 3D structure\n"; 
        exit(1);
    }

//  Open File
    filename.append(tag).append(oFextension(spec,step));
    oFile.open(filename.c_str()); 
    oFile.flush();

//  Export dimensions
    oFile << Hdr[tag].dim() << "\n";

//  Export time
    oFile << Hdr[tag].Time_label() << "\n";
    //oFile << float(step)*Hdr[tag].Time_conv() << "\n";  // 
    oFile << float(time)*Hdr[tag].Time_conv() << "\n";  // 

//  Export all axes
    for (size_t i(0); i< Hdr[tag].dim(); ++i) {
        oFile << Hdr[tag].label(i) << "\n"; 
        oFile << Hdr[tag].axis(i);
    }
//  Renormalize and export data
    oFile << Hdr[tag].Title_label() << "\n"; 
    data *= Hdr[tag].Title_conv();         
    oFile << data; 

    oFile.flush();
    oFile.close();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Adjust filenames with zeros to reach some prescribed length. 
//  Add the filename extension.  
string Export_Files_U1::Xport::oFextension(size_t species, size_t step){
     stringstream sFilename;
     sFilename << "_s" << species <<"_";
        
     // Number of zeros to add to the filename
     int Nzeros(ofconventions::ofile_digits - stringify(step).length()); 
     while (Nzeros-- > 0) { 
         sFilename << "0"; 
     }
        
     sFilename << step << ofconventions::ofile_extension;

     return sFilename.str();
}
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
    Export_Files_U1::Restart_Facility::Restart_Facility(string homedir) {
        hdir = homedir;
//        if (makefolder(hdir+"RESTART/") != 0) cout<<"Warning: Folder "<< hdir+"RESTART/"<<" exists\n";
    }

//--------------------------------------------------------------
//  Read restart file
    void Export_Files_U1::Restart_Facility::Read(const int rank, const size_t re_step, LocalState_U1& Y) {

//      Generate filename 
        string   filename(hdir+"RESTART/re_");
        filename.append(rFextension(rank,re_step));

//      Open file
        ifstream  fin(filename.c_str(), ios::binary);

//      Read spin distribution
    for (unsigned i(0); i < Y.SD().size(); ++i) {
            fin.read( (char *)(&(*Y.SD().sd)[i]), sizeof((*Y.SD().sd)[i]) );
        }
            
    //for (unsigned i(0); i < Y.SD().dim(); ++i) {
        //    fin.read( (char *)(&Y.SD().real(i)), sizeof(Y.SD().real(i)) );
        //    fin.read( (char *)(&Y.SD().imag(i)), sizeof(Y.SD().imag(i)) );
        //}
            
//      Read H
        for(unsigned i(0); i < Y.MF().dim(); ++i){
            fin.read( (char *)(&Y.MF()(i)), sizeof(Y.MF()(i)) );
        }

        fin.close();
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Write restart file
    void Export_Files_U1::Restart_Facility::Write(const int rank, const size_t re_step, LocalState_U1& Y) {

//      Make restart folder if there is none
        Makefolder(hdir+"RESTART/");

//      Generate filename 
        string   filename(hdir+"RESTART/re_");
        filename.append(rFextension(rank,re_step));

//      Open file
        ofstream  fout(filename.c_str(), ios::binary);

//      Write spin distribution
    for (unsigned i(0); i < Y.SD().size(); ++i) {
            fout.write( (char *)(&(*Y.SD().sd)[i]), sizeof((*Y.SD().sd)[i]) );
            //fout.write( (char *)(&Y.SD().real(i)), sizeof(Y.SD().real(i)) );
            //fout.write( (char *)(&Y.SD().imag(i)), sizeof(Y.SD().imag(i)) );
        }
            
//      Read H
        for(unsigned i(0); i < Y.MF().dim(); ++i){
            fout.write( (char *)(&Y.MF()(i)), sizeof(Y.MF()(i)) );
        }

        fout.flush();
        fout.close();
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Adjust filenames with zeros to reach some prescribed length. 
//  Add the filename extension. 
string Export_Files_U1::Restart_Facility::rFextension(const int rank, const size_t rstep){
        stringstream sFilename;

        // Number of zeros to add to the filename
        int Nzeros(ofconventions::rank_digits - stringify(rank).length()); 
        while (Nzeros-- > 0) { 
            sFilename << "0"; 
        } 
        
        sFilename << rank << "_";

        // Number of zeros to add to the filename
        Nzeros = ofconventions::rfile_digits - stringify(rstep).length(); 
        while (Nzeros-- > 0) { 
            sFilename << "0";
        } 
        
        sFilename << rstep << ofconventions::rfile_extension;

        return sFilename.str();
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Output_Data_U1::ShellDensity::ShellDensity(const Export_Files_U1::Grid_Data& _G)
    : cos8( size_t(_G.Ncos8) ), phi( size_t(_G.Nphi) ),
      PLegendre( size_t(_G.Ncos8), size_t(_G.Nl), size_t(_G.Nm) ),
      COSmphi(   size_t(_G.Nphi), size_t(_G.Nm)  ),
      SINmphi(   size_t(_G.Nphi), size_t(_G.Nm)  ){

//  Make cos(8) axis
    float cos8min(-1.0+1.0/float(_G.Ncos8));
    float cos8max( 1.0-1.0/float(_G.Ncos8));
    cos8 = Algorithms::MakeAxis(cos8min, cos8max, size_t(_G.Ncos8));

//  Get legendre polynomials for each cosine
    for (size_t i(0); i < PLegendre.dim1(); ++i) {
        Array2D<float> aL( Algorithms::Legendre( cos8[i], size_t(_G.Nl), size_t(_G.Nm) ) );

//  Copy legendre polynomials into the array
        for (size_t l(0); l < PLegendre.dim2(); ++l) {
            for (size_t m(0); m < PLegendre.dim3(); ++m) {
                PLegendre(i,l,m) = aL(l,m); 
            }
        }
    }
    
//  Make phi axis
    float phimin( M_PI/float(_G.Nphi));
    float phimax( 2*M_PI-M_PI/float(_G.Nphi));
    phi =  Algorithms::MakeAxis( phimin,phimax, size_t(_G.Nphi)) ;

    for (size_t i(0); i < COSmphi.dim1(); ++i) {
        for (size_t m(0); m < COSmphi.dim2(); ++m) {
            COSmphi(i,m) = cos(m*phi[i]);
            SINmphi(i,m) = sin(m*phi[i]);
        }
    }
 
} 

float Output_Data_U1::ShellDensity::COS8(const size_t _cos8) const { 
    return cos8[_cos8]; 
}
float Output_Data_U1::ShellDensity::PHI(const size_t _phi) const { 
    return phi[_phi]; 
}
float Output_Data_U1::ShellDensity::PL(const size_t _cos8, const size_t _l, const size_t _m) const { 
    return PLegendre( _cos8, _l, _m); 
}
float Output_Data_U1::ShellDensity::cosMphi (const size_t _phi, const size_t _m) const {
    return COSmphi( _phi, _m);
}
float Output_Data_U1::ShellDensity::sinMphi (const size_t _phi, const size_t _m) const {
    return SINmphi( _phi, _m);
}
//--------------------------------------------------------------
//**************************************************************

Output_Data_U1::Output_Preprocessor::Output_Preprocessor(const Export_Files_U1::Grid_Data& _gdata,  
                            //const vector< string > _oTags, 
                                                      string homedir)
        : expo( _gdata, //_oTags, 
                 homedir),
          sd( _gdata )
          //oTags(_oTags) 
          { }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void Output_Data_U1::Output_Preprocessor::operator()(
    const LocalState_U1& Y, const size_t tout, const double _time) { 

    float YSH_re(0.0), YSH_im(0.0), mphi(0.0);
    Array2D<float> phicos8( size_t(sd.Nphi()), size_t(sd.Ncos8()));

//  Loop for each location in phi-theta
//  Loop for each harmonic
    for (size_t i(0); i < sd.Nphi(); ++i) { 

        for (size_t l(0); l < sd.Nl(); ++l) {
            for (size_t m(0); m < sd.Nm(); ++m) {
                YSH_re = sd.cosMphi(i,m) * float(Y.SD().real(l,m));
                YSH_im = sd.sinMphi(i,m) * float(Y.SD().imag(l,m));
                YSH_re -= YSH_im; 

                for(size_t j(0); j < sd.Ncos8(); ++j) {
                    phicos8(i,j) += float(2.0) * sd.PL(j,l,m) * YSH_re;
                }

            }
        }

    }
    expo("Distribution", phicos8, tout, _time);
}
//--------------------------------------------------------------
//**************************************************************
