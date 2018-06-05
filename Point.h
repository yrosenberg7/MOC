#ifndef POINT_H
#define POINT_H

#include "Gas_Model.h"
#include <string>
#include <stdio.h>
#include <iostream>
/** WARNING: INCOMPLETE*/
class Point
{
    public:
        Point();
        Point(double,double); /// create point at location
        Point(double,double, double, double); /// create point at location, with velocity
        Point(const Point& other);

        virtual ~Point(){}

        void print();
        void print_xm();
        void print(FILE*);
        ///void print(&file); // prints to file directly?

        /// COMPARES ONLY ID!
        bool operator==(const Point&)const;
        bool operator!=(const Point&)const;
        bool operator>(const Point&)const;
        bool operator>=(const Point&)const;
        bool operator<(const Point&)const;
        bool operator<=(const Point&)const;


        void finish_pt(Gas_Model*);

        enum LOC
        {
            NA = 0,
            THROAT = 1,
            IDL = 2,
            IDL_R = 3,
            FIXED_W = 4,
            VAR_W = 5,
        };

        void set_x(double); /// change to address?
        void set_y(double);
        void set_u(double);
        void set_v(double);
        void set_a(double);
        void set_T(double);
        void set_p(double);
        void set_rho(double);
        void set_M(double);
        void set_uv(double, double);
        ///void set_loc();

        double get_x()const; /// change to pointers?
        double get_y()const;
        double get_u()const;
        double get_v()const;
        double get_a()const;
        double get_T()const;
        double get_p()const;
        double get_rho()const;
        double get_M()const;
        int get_id() const;
        ///Point::LOC get_loc() const;

        //bool operator==(const Point& b) const;
        //bool operator==(const Point& a, const Point& b);

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

        ///Point::LOC tag;
        ///string tags[6] = {"NA","Throat","IDL","IDL_R","Fixed_W","Var_W"}

};

/*
class Tri: public Point
{
    public:
        std:string print();
        ///void print(&file)
};
*/
#endif // POINT_H
