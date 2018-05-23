#include "Throat.h"
//#include "y_math.h"
#include <cmath>

#include <iostream> ///REMOVE?

#define PI 3.14159265

///TODO: add variable wall capability
//TODO: change print to take a file


Throat::Throat()
{
    this->upstream_rth = 2.0;
    this->downstream_rth = 0.5;
    this->height = 1;
    this->max_angle = 13.65*PI/180.0;
    this->xbreak = downstream_rth*sin(max_angle);
    this->axi = true;
    this->ysmooth = true;


    if (axi)
    {
        this->delta = 1;
    }else{
        this->delta = 0;
    }
}

// constructor, with max angle in degrees
Throat::Throat(double upstream_rth, double downstream_rth, double height, double max_angle, bool axi, bool ysmooth)
{
    this->upstream_rth = upstream_rth;
    this->downstream_rth = downstream_rth;
    this->height = height;
    this->max_angle = max_angle*PI/180.0;
    this->axi = axi; ///FIX!
    this->ysmooth = ysmooth; ///FIX!

    this->xbreak = downstream_rth*sin(max_angle);

    if (axi)
    {
        this->delta = 1;
    }else{
        this->delta = 0;
    }
}

//Throat::~Throat(){}

void Throat::print()
{
    std::cout << "THROAT GEOM: upstream rad: " << this->get_upstream_rth() << ",\t"
              << "downstream rad: " << this->get_downstream_rth() << ",\t"
              << "height: " << this->get_height() << ",\t"
              << "max angle: " << this->get_max_angle() << ",\t"
              << "axi: " << this->get_axi() << ",\t"
              << "ymooth: " << this->get_ysmooth() << std::endl;
}

///public
double Throat::wall_end()
{
    return this->xbreak;
}

///public
double Throat::fwall(double x)
{
    double ybreak = wall_help(this->xbreak);
    double mbreak = wall_slope(this->xbreak, this->max_angle, this->downstream_rth);

    if(x<=this->xbreak)
    {
        return wall_help(x);
    }else
    {
        return ybreak+mbreak*(x - this->xbreak);
    }
}

///public
double Throat::fpwall(double x)
{
    return wall_slope(x, this->max_angle, this->downstream_rth);
}

///private
double Throat::wall_help(double x)
{
    return this->height + this->downstream_rth - sqrt(pow(this->downstream_rth,2) - pow(x,2));
}

///private
double Throat::slope_help(double x)
{
    double x2 = pow(x,2);
    double down2 = pow(this->downstream_rth, 2);
    double diff = down2 - x2;
    double sqr = sqrt(diff);
    double result = x/sqr;

/*
    printf("\n#############################\nslope help: \nx^2: %f"
            "\ndownstream rad^2: %f"
            "\ndiff: %f, \nsqrt: %f,"
            "\nresult: %f", x2,down2,diff,sqr,result);
*/
    return result;
    //return x/sqrt(pow(this->downstream_rth,2) - pow(x,2));
}

///private
double Throat::wall_slope(double x, double max_ange, double downstream_rth)
{
    //double xbreak = downstream_rth*sin(max_angle);
    if(x<=this->xbreak)
    {
        return slope_help(x);
    }else
    {
        return slope_help(this->xbreak);
    }

}

double Throat::get_height(){return this->height;}
double Throat::get_max_angle(){return this->max_angle;}
double Throat::get_upstream_rth(){return this->upstream_rth;}
double Throat::get_downstream_rth(){return this->downstream_rth;}
bool Throat::get_axi(){return this->axi;}
bool Throat::get_ysmooth(){return this->ysmooth;}
int Throat::get_delta(){return this->delta;}
