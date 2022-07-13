///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Mar 15, 2013
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the 
//   structures that are required for exporting data:
//
//   1. Xport:
//     Output data
//
//   2. Restart Facility:
//     Restart data
///////////////////////////////////////////////////////////

    #ifndef DECL_EXPORT_H
    #define DECL_EXPORT_H



//**************************************************************
//--------------------------------------------------------------
namespace Export_Files{

    namespace ofconventions {
        const int    ofile_digits = 5;
        const string ofile_extension = ".txt";
        const int    ofile_precision = 6; 

        const int    rfile_digits = 3;
        const int    rank_digits = 6;
        const string rfile_extension = ".dat";
    }
//--------------------------------------------------------------

//  Folder and filename functions
    template<typename T> inline std::string stringify(T const& x) {
        std::ostringstream out;
        out << x;
        return out.str();
    }

    void Makefolder(string _name);
//--------------------------------------------------------------

//  Contains all the info you need to construct an output axis
    class DefaultTags {
        public:
//          Contents
        vector<string> time;
        vector<string> space;
        vector<string> momfld;
        vector<string> pvsx;
        
        DefaultTags(size_t species); 
        ~DefaultTags(){}
    };
//--------------------------------------------------------------

//  redefinitions of "<<" operator for data containers
    template <class T> 
    ofstream& operator<<(ofstream& s, const vector<T>& v) {
        s << setprecision(ofconventions::ofile_precision);
        s << 1 <<"\n";
        s << v.size()<<"\n";
        for (size_t i(0); i < v.size(); ++i) {
            s << v[i]<<"\n";
        }
        return s;
    }

    template <class T> 
    ofstream& operator<<(ofstream& s, const valarray<T>& v) {
        s << setprecision(ofconventions::ofile_precision);
        s << 1 <<"\n";
        s << v.size()<<"\n";
        for (size_t i(0); i < v.size(); ++i) {
            s << v[i]<<"\n";
        }
        return s;
    }

    template <class T> 
    ofstream& operator<<(ofstream& s, const Array2D<T>& array2D) {
        s << setprecision(ofconventions::ofile_precision);
        s << 2 <<"\n";
        s << array2D.dim1()<<"\n";
        s << array2D.dim2()<<"\n";
        for (size_t i(0); i < array2D.dim(); ++i) {
            s << array2D(i)<<"\n";
        }
        return s;
    }

    template <class T> 
    ofstream& operator<<(ofstream& s, const Array3D<T>& array3D) {
        s << setprecision(ofconventions::ofile_precision);
        s << 3 <<"\n";
        s << array3D.dim1()<<"\n";
        s << array3D.dim2()<<"\n";
        s << array3D.dim3()<<"\n";
        for (size_t i(0); i < array3D.dim(); ++i) {
            s << array3D(i)<<"\n";
        }
        return s;
    }
//--------------------------------------------------------------

//  Convert data structure to float structure
    valarray<float> vfloat(const valarray<double>& vDouble); 
    vector<float>   vfloat(const vector<double> vDouble); 


//--------------------------------------------------------------

//  Contains all the info you need to construct an output axis
    class oAxis {
        public:
//          Contents
            string label;           // e.g. label = cm 
            float min, max;  
            size_t sz;

//          Constructors
            oAxis();
            oAxis(const float _m, const float _M, const size_t _sz); 
            oAxis(const string _l, const float _m, const float _M, 
                  const size_t _sz);

//          Copy constructor 
            oAxis(const oAxis& other);
            ~oAxis(){}
    };
//--------------------------------------------------------------

//  Facilitates the generation of a header
    class Header {
        public:
//          Constructor
            Header() { };
            Header(oAxis _x,                                        // 1D
                   string _Ql, float _Qc, string _tl, float _tc, string _oD);
            Header(oAxis _x, oAxis _y,                             // 2D 
                   string _Ql, float _Qc,  string _tl, float _tc, string _oD);
            Header(oAxis _x, oAxis _y, oAxis _z,                  // 3D
                   string _Ql, float _Qc, string _tl, float _tc, string _oD);
            Header(vector< oAxis > _xyz,                            // xD
                   string _Ql, float _Qc, string _tl, float _tc, string _oD);
            size_t dim();    

            valarray<float> axis(const size_t i); // this is axis 0, 1, 2 
            string          label(const size_t i);
            float           conv(const size_t i);
            
            string  Title_label(); 
            float   Title_conv(); 
            string  Time_label();
            float   Time_conv();
            string  Directory();
        
        private:
            vector< oAxis >    xyz_axis;         // axis 0, 1, 2  : size 0-2
            string             title,  time;
            float              titleC, timeC;
            string 	       oDir;
    };
//--------------------------------------------------------------

class Grid_Data {
public:
//  Constructor
    Grid_Data(const unsigned _Nl,    const unsigned _Nm, 
	      const unsigned _Ncos8, const unsigned _Nphi) 
	      : Nl(_Nl), Nm(_Nm), Ncos8(_Ncos8), Nphi(_Nphi){ }

//  Copy constructor
    Grid_Data(const Grid_Data& other): Nl(other.Nl),       Nm(other.Nm), 
                                       Ncos8(other.Ncos8), Nphi(other.Nphi) { }
    ~Grid_Data(){};

    const unsigned Nl;
    const unsigned Nm;
    const unsigned Ncos8;
    const unsigned Nphi;
};
//--------------------------------------------------------------

//  Main facility for exporting data 
    class Xport {
        public:
//          Constructor
            Xport(const Grid_Data& _gdata, //const size_t Nphi, const size_t Ncos8, 
                  // const vector< string > oTags,
                  string homedir=""); 

//	    export 1D data structure
            void operator() (const string tag,  valarray<float> data,
                             const size_t step, const double time, 
                             const size_t spec = 0);  
//	           2D data structure
            void operator() (const string tag,  Array2D<float> data, 
                             const size_t step, const double time, 
                             const size_t spec = 0); 
//	           3D data structure
            void operator() (const string tag,  Array3D<float> data, 
                             const size_t step, const double time,
                             const size_t spec = 0);
        
        private:
            map< string, Header > Hdr; // Dictionary of headers
            string oFextension(size_t species, size_t step);
    };
//--------------------------------------------------------------

    class Restart_Facility {

    public:
        Restart_Facility(string homedir="");
        
        void Read(const int rank, const size_t re_step, LocalState& Y);
        void Write(const int rank, const size_t re_step, LocalState& Y);

    private:
        string hdir;
        string rFextension(const int rank, const size_t rstep);
    }; 
//--------------------------------------------------------------
}
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
namespace Output_Data{


//--------------------------------------------------------------
    class ShellDensity {
    public:
//      Constructor/Destructor
        ShellDensity(const Export_Files::Grid_Data& _G);
        ShellDensity(const ShellDensity& other);
        ~ShellDensity(){ }

//      Axis
        float COS8(const size_t _cos8) const;
        float PHI(const size_t _phi) const;

//      Functions 
        float PL(const size_t _cos8, const size_t _l, const size_t _m) const; 
        float cosMphi (const size_t _phi, const size_t _m) const;
        float sinMphi (const size_t _phi, const size_t _m) const;

//      Shape
	unsigned Ncos8() const { return unsigned( PLegendre.dim1() ); }
	unsigned Nl()    const { return unsigned( PLegendre.dim2() ); }
	unsigned Nm()    const { return unsigned( PLegendre.dim3() ); }
	unsigned Nphi()  const { return unsigned( COSmphi.dim1()   ); }

    private:
	Array3D<float> PLegendre;
	Array2D<float> COSmphi, SINmphi;

        valarray<float> cos8, phi;
    };
//--------------------------------------------------------------

//  Output Functor
    class Output_Preprocessor {
    public:
//      Constructor
        Output_Preprocessor(const Export_Files::Grid_Data& _gdata,  
                          //  const vector< string > _oTags, 
                            string homedir="");  
//      Functor
        void operator()(const LocalState& Y, const size_t tout, const double _time);

    private:
        Export_Files::Xport expo;
        ShellDensity        sd;
        // vector< string >    oTags;
    };
//--------------------------------------------------------------

}
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
namespace Export_Files_U1{

    namespace ofconventions {
        const int    ofile_digits = 5;
        const string ofile_extension = ".txt";
        const int    ofile_precision = 6; 

        const int    rfile_digits = 3;
        const int    rank_digits = 6;
        const string rfile_extension = ".dat";
    }
//--------------------------------------------------------------

//  Folder and filename functions
    template<typename T> inline std::string stringify(T const& x) {
        std::ostringstream out;
        out << x;
        return out.str();
    }

    void Makefolder(string _name);
//--------------------------------------------------------------

//  Contains all the info you need to construct an output axis
    class DefaultTags {
        public:
//          Contents
        vector<string> time;
        vector<string> space;
        vector<string> momfld;
        vector<string> pvsx;
        
        DefaultTags(size_t species); 
        ~DefaultTags(){}
    };
//--------------------------------------------------------------

//  redefinitions of "<<" operator for data containers
    template <class T> 
    ofstream& operator<<(ofstream& s, const vector<T>& v) {
        s << setprecision(ofconventions::ofile_precision);
        s << 1 <<"\n";
        s << v.size()<<"\n";
        for (size_t i(0); i < v.size(); ++i) {
            s << v[i]<<"\n";
        }
        return s;
    }

    template <class T> 
    ofstream& operator<<(ofstream& s, const valarray<T>& v) {
        s << setprecision(ofconventions::ofile_precision);
        s << 1 <<"\n";
        s << v.size()<<"\n";
        for (size_t i(0); i < v.size(); ++i) {
            s << v[i]<<"\n";
        }
        return s;
    }

    template <class T> 
    ofstream& operator<<(ofstream& s, const Array2D<T>& array2D) {
        s << setprecision(ofconventions::ofile_precision);
        s << 2 <<"\n";
        s << array2D.dim1()<<"\n";
        s << array2D.dim2()<<"\n";
        for (size_t i(0); i < array2D.dim(); ++i) {
            s << array2D(i)<<"\n";
        }
        return s;
    }

    template <class T> 
    ofstream& operator<<(ofstream& s, const Array3D<T>& array3D) {
        s << setprecision(ofconventions::ofile_precision);
        s << 3 <<"\n";
        s << array3D.dim1()<<"\n";
        s << array3D.dim2()<<"\n";
        s << array3D.dim3()<<"\n";
        for (size_t i(0); i < array3D.dim(); ++i) {
            s << array3D(i)<<"\n";
        }
        return s;
    }
//--------------------------------------------------------------

//  Convert data structure to float structure
    valarray<float> vfloat(const valarray<double>& vDouble); 
    vector<float>   vfloat(const vector<double> vDouble); 


//--------------------------------------------------------------

//  Contains all the info you need to construct an output axis
    class oAxis {
        public:
//          Contents
            string label;           // e.g. label = cm 
            float min, max;  
            size_t sz;

//          Constructors
            oAxis();
            oAxis(const float _m, const float _M, const size_t _sz); 
            oAxis(const string _l, const float _m, const float _M, 
                  const size_t _sz);

//          Copy constructor 
            oAxis(const oAxis& other);
            ~oAxis(){}
    };
//--------------------------------------------------------------

//  Facilitates the generation of a header
    class Header {
        public:
//          Constructor
            Header() { };
            Header(oAxis _x,                                        // 1D
                   string _Ql, float _Qc, string _tl, float _tc, string _oD);
            Header(oAxis _x, oAxis _y,                             // 2D 
                   string _Ql, float _Qc,  string _tl, float _tc, string _oD);
            Header(oAxis _x, oAxis _y, oAxis _z,                  // 3D
                   string _Ql, float _Qc, string _tl, float _tc, string _oD);
            Header(vector< oAxis > _xyz,                            // xD
                   string _Ql, float _Qc, string _tl, float _tc, string _oD);
            size_t dim();    

            valarray<float> axis(const size_t i); // this is axis 0, 1, 2 
            string          label(const size_t i);
            float           conv(const size_t i);
            
            string  Title_label(); 
            float   Title_conv(); 
            string  Time_label();
            float   Time_conv();
            string  Directory();
        
        private:
            vector< oAxis >    xyz_axis;         // axis 0, 1, 2  : size 0-2
            string             title,  time;
            float              titleC, timeC;
            string         oDir;
    };
//--------------------------------------------------------------

    class Grid_Data {
    public:
//      Constructor
        Grid_Data(const unsigned _Nl,    const unsigned _Nm, 
          const unsigned _Ncos8, const unsigned _Nphi) 
          : Nl(_Nl), Nm(_Nm), Ncos8(_Ncos8), Nphi(_Nphi){ }

//  Copy constructor
        Grid_Data(const Grid_Data& other): Nl(other.Nl),       Nm(other.Nm), 
                                       Ncos8(other.Ncos8), Nphi(other.Nphi) { }
        ~Grid_Data(){};

        const unsigned Nl;
        const unsigned Nm;
        const unsigned Ncos8;
        const unsigned Nphi;
    };
//--------------------------------------------------------------

//  Main facility for exporting data 
    class Xport {
        public:
//          Constructor
            Xport(const Grid_Data& _gdata, //const size_t Nphi, const size_t Ncos8, 
                  // const vector< string > oTags,
                  string homedir=""); 

//      export 1D data structure
            void operator() (const string tag,  valarray<float> data,
                             const size_t step, const double time, 
                             const size_t spec = 0);  
//             2D data structure
            void operator() (const string tag,  Array2D<float> data, 
                             const size_t step, const double time, 
                             const size_t spec = 0); 
//             3D data structure
            void operator() (const string tag,  Array3D<float> data, 
                             const size_t step, const double time,
                             const size_t spec = 0);
        
        private:
            map< string, Header > Hdr; // Dictionary of headers
            string oFextension(size_t species, size_t step);
    };
//--------------------------------------------------------------

    class Restart_Facility {

    public:
        Restart_Facility(string homedir="");
        
        void Read(const int rank, const size_t re_step, LocalState_U1& Y);
        void Write(const int rank, const size_t re_step, LocalState_U1& Y);

    private:
        string hdir;
        string rFextension(const int rank, const size_t rstep);
    }; 
//--------------------------------------------------------------
}
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
namespace Output_Data_U1{


//--------------------------------------------------------------
    class ShellDensity {
    public:
//      Constructor/Destructor
        ShellDensity(const Export_Files_U1::Grid_Data& _G);
        ShellDensity(const ShellDensity& other);
        ~ShellDensity(){ }

//      Axis
        float COS8(const size_t _cos8) const;
        float PHI(const size_t _phi) const;

//      Functions 
        float PL(const size_t _cos8, const size_t _l, const size_t _m) const; 
        float cosMphi (const size_t _phi, const size_t _m) const;
        float sinMphi (const size_t _phi, const size_t _m) const;

//      Shape
    unsigned Ncos8() const { return unsigned( PLegendre.dim1() ); }
    unsigned Nl()    const { return unsigned( PLegendre.dim2() ); }
    unsigned Nm()    const { return unsigned( PLegendre.dim3() ); }
    unsigned Nphi()  const { return unsigned( COSmphi.dim1()   ); }

    private:
    Array3D<float> PLegendre;
    Array2D<float> COSmphi, SINmphi;

        valarray<float> cos8, phi;
    };
//--------------------------------------------------------------

//  Output Functor
    class Output_Preprocessor {
    public:
//      Constructor
        Output_Preprocessor(const Export_Files_U1::Grid_Data& _gdata,  
                          //  const vector< string > _oTags, 
                            string homedir="");  
//      Functor
        void operator()(const LocalState_U1& Y, const size_t tout, const double _time);

    private:
        Export_Files_U1::Xport expo;
        ShellDensity        sd;
        // vector< string >    oTags;
    };
//--------------------------------------------------------------

}
//**************************************************************


    #endif
