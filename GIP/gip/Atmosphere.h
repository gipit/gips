#ifndef ATMOSPHERE_INCLUDED
#define ATMOSPHERE_INCLUDED

namespace gip {

    //! Represents an atmosphere holding various parameters (transmission, radiance, etc)
    class Atmosphere {
    public:
        //! Default constructor
        Atmosphere() :
            Xa(0.0), Xb(0.0), Xc(0.0),
            _Valid(false), _t(0), _Lu(0), _Ld(0) {}
        //! Create new atmosphere with parameters
        Atmosphere(double t, double Lu, double Ld)
            : _Valid(true), _t(t), _Lu(Lu), _Ld(Ld) {}
        //! Copy constructor
        Atmosphere(const Atmosphere& atm)
            : Xa(atm.Xa), Xb(atm.Xb), Xc(atm.Xc),
            _Valid(atm._Valid), _t(atm._t), _Lu(atm._Lu), _Ld(atm._Ld) {}
        ~Atmosphere() {};

        bool Valid() const { return _Valid; }
        double t() const { return _t; }
        double Lu() const { return _Lu; }
        double Ld() const { return _Ld; }

        // Some temp stuff
        bool Coef() const {
            return (Xa*Xb*Xc == 0.0) ? false : true;
        }
        double Xa, Xb, Xc;

    private:
        bool _Valid;
        double _t;
        double _Lu;
        double _Ld;
        //double _e;
    };

}

#endif // ATMOSPHERE_INCLUDED
