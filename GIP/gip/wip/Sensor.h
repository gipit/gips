#ifndef SENSOR_INCLUDED
#define SENSOR_INCLUDED

namespace gip {

    //! Represents a sensor and holds sensor specific parameters used in processing
    class Sensor {
    public:
        //! Default constructor
        Sensor()
            //: _minDC(0), _maxDC(255), _RefCoef(1.0), _K1(0), _K2(0) {}
            : _minDC(0), _maxDC(255), _K1(0), _K2(0) {}
        //! Copy constructor
        Sensor(const Sensor& s)
            //: _minDC(s._minDC), _maxDC(s._maxDC), _Radiance(s._Radiance), _RefCoef(s._RefCoef), _K1(s._K1), _K2(s._K2) {}
            : _minDC(s._minDC), _maxDC(s._maxDC), _Radiance(s._Radiance), _K1(s._K1), _K2(s._K2) {}
        ~Sensor() {};

        // Get functions
        int minDC() const { return _minDC; }
        int maxDC() const { return _maxDC; }
        double K1() const { return _K1; }
        double K2() const { return _K2; }
        //double RefCoef() const { return _RefCoef; }
        double Radiance() const { return _Radiance; }

        //! Sets dyanmic range of sensor (min to max digital counts)
        void SetDynamicRange(int min, int max) {
            _minDC = min;
            _maxDC = max;
        }

        //! Set total, either exoatmospheric or downward on target
        void SetTotalRadiance(float rad) {
            _Radiance = rad;
        }



		//! Is this a thermal sensor band?
		bool Thermal() const { if (_K1*_K2 == 0) return false; else return true; }
		//! Set thermal band
		void SetThermal(float k1=666.09, float k2=1282.71) { _K1=k1; _K2=k2; }

    private:
        int _minDC;
        int _maxDC;
        double _Radiance;
        //double _RefCoef;
        double _K1;
        double _K2;
    };

}



#endif // SENSOR_INCLUDED
