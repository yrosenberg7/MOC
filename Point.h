#ifndef POINT_H
#define POINT_H

#include "Gas_Model.h"
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>

class Point
{
    public:
        Point();
        Point(double,double); /// create point at location
        Point(double,double, double, double); /// create point at location, with velocity
        Point(const Point& other); /// copy

        virtual ~Point(){}

        void print(); /// to cout
        void print_xm();
        void print_xu();
        void print(std::ostream*); /// to file

        /// COMPARES ONLY ID
        bool operator==(const Point&)const;
        bool operator!=(const Point&)const;
        bool operator>(const Point&)const;
        bool operator>=(const Point&)const;
        bool operator<(const Point&)const;
        bool operator<=(const Point&)const;

        void finish_pt(Gas_Model*);

        void set_x(double);
        void set_y(double);
        void set_u(double);
        void set_v(double);
        void set_a(double);
        void set_T(double);
        void set_p(double);
        void set_rho(double);
        void set_M(double);
        void set_uv(double, double);

        double get_x()const;
        double get_y()const;
        double get_u()const;
        double get_v()const;
        double get_a()const;
        double get_T()const;
        double get_p()const;
        double get_rho()const;
        double get_M()const;
        int get_id() const;

    protected:

    private:
        static int id;
        int cid;

        double x;
        double y;
        double u;
        double v;
        double a;
        double T;
        double p;
        double rho;
        double M;
};

#endif // POINT_H
