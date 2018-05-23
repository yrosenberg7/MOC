#ifndef THROAT_H
#define THROAT_H

#include <iostream>

#include "Point.h"
#include "Gas_Model.h"
/** WARNING: UNTESTED*/
class Throat
{
    public:
        Throat();
        Throat(double, double, double, double, bool, bool);
        //virtual ~Throat(){std::cout << "Throat Deconstructed" << std::endl;};
        virtual ~Throat(){};

        double fwall(double);
        double fpwall(double);
        double wall_end();

        void print();

        double get_height();
        double get_max_angle();
        double get_upstream_rth();
        double get_downstream_rth();
        bool get_axi();
        bool get_ysmooth();
        int get_delta();

    protected:
        double upstream_rth;
        double downstream_rth;
        double height;
        double max_angle; // in radians
        double xbreak;
        int delta;
        bool axi;
        bool ysmooth;

        double wall_help(double); // used in fwall function
        double slope_help(double);
        double wall_slope(double, double, double);

    private:


        //double wall_shape(double, double, double, double, bool);

};

#endif // THROAT_H
