#include "Point.h"
#include <cmath>
#include <string>
#include <iostream>
/** WARNING: INCOMPLETE*/


int Point::id = 0; //init object counter

Point::Point()
{
    this->x = 0.0;
    this->y = 0.0;
    this->u = 0.0;
    this->v = 0.0;
    this->a = 0.0;
    this->T = 0.0;
    this->p = 0.0;
    this->rho = 0.0;
    this->M = 0.0;

    this->cid = 0;
    //this->tag = NA;
}
/// create point at location
Point::Point(double xin,double yin)
{
    this->x = xin;
    this->y = yin;
    this->u = 0;
    this->v = 0;
    this->a = 0;
    this->T = 0;
    this->p = 0;
    this->rho = 0;
    this->M = 0;

    this->cid = id++;
    //this->tag = NA;
}

Point::Point(double xin,double yin, double uin, double vin)
{
    this->x = xin;
    this->y = yin;
    this->u = uin;
    this->v = vin;
    this->a = 0;
    this->T = 0;
    this->p = 0;
    this->rho = 0;
    this->M = 0;

    this->cid = id++;
    //this->tag = NA;
}

Point::Point(const Point& other)
{
//    std::cout << "COPY POINT!" << std::endl;
    this->x = other.get_x();
    this->y = other.get_y();
    this->u = other.get_u();
    this->v = other.get_v();
    this->a = other.get_a();
    this->T = other.get_T();
    this->p = other.get_p();
    this->rho = other.get_rho();
    this->M = other.get_M();

    //this->cid = id++;
    this->cid = other.get_id();
}

//Point::~Point()
//{
//    dtor
//}

//std::string Point::print()
void Point::print()
{
    std::cout << ""
    << this->get_id() <<"\t"
    << this->get_x() <<"\t"
    << this->get_y() <<"\t"
    << this->get_u() <<"\t"
    << this->get_v() <<"\t"
    << this->get_a() <<"\t"
    << this->get_T() <<"\t"
    << this->get_p() <<"\t"
    << this->get_rho() <<"\t"
    << this->get_M() <<"\t"
    << std::endl;
}

void Point::print_xm()
{
    std::cout << ""
    <<"x: " << this->get_x() <<"\t\t"
    <<"y: "<< this->get_y() <<"\t\t"
    <<"M: "<< this->get_M() <<"\t\t"
    << std::endl;
}
void Point::print_xu()
{
    std::cout <<""
    <<"x: "<< this->get_x() <<"\t\t"
    <<"y: "<< this->get_y() <<"\t\t"
    <<"u: "<< this->get_u() <<"\t\t"
    <<"v: "<< this->get_v() <<"\t\t"
    << std::endl;
}
void Point::print(FILE* fout)
{
    fprintf(fout, "%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n",this->get_id(),this->get_x(),
            this->get_y(),this->get_u(),this->get_v(),
            this->get_a(),this->get_T(),this->get_p(),this->get_rho(),this->get_M());
}



void Point::finish_pt(Gas_Model* GM)
{
    //GM->print();

//    #warning this just fills them with NAN!
    double U0 = sqrt(pow(this->get_u(),2) + pow(this->get_v(),2));

    //double aa = GM->a_from_u(U0);
    this->set_a(GM->a_from_u(U0));

    //std::cout << "a: " << aa << std::endl; /// why is a defined ?!?!?!?!?!?!?!
    //this->set_a(aa);
    this->set_T(GM->T_from_a(this->get_a()));
    this->set_p(GM->p_from_T(this->get_T()));
    this->set_rho(GM->rho_from_pT(this->get_p(),this->get_T()));
    this->set_M(U0/this->get_a());

}

//bool Point::operator==(const Point& a, const Point& b)
//bool Point::operator==(const Point& b) const
//{
//    //return ((b.get_x() == this->get_x()) && (b.get_y() == this->get_y()));
//    return (this->get_x() == b.get_x()) && (this->get_y() == b.get_y());
//}



/// setters
void Point::set_x(double xin){this->x = xin;}
void Point::set_y(double yin){this->y = yin;}
void Point::set_u(double uin){this->u = uin;}
void Point::set_v(double vin){this->v = vin;}
void Point::set_a(double ain){this->a = ain;}
void Point::set_T(double Tin){this->T = Tin;}
void Point::set_p(double pin){this->p = pin;}
void Point::set_rho(double rhoin){this->rho = rhoin;}
void Point::set_M(double Min){this->M = Min;}
void Point::set_uv(double uin, double vin)
{
    this->u = uin;
    this->v = vin;
}

///getters
double Point::get_x()const {return this->x;}
double Point::get_y()const {return this->y;}
double Point::get_u()const {return this->u;}
double Point::get_v()const {return this->v;}
double Point::get_a()const {return this->a;}
double Point::get_T()const {return this->T;}
double Point::get_p()const {return this->p;}
double Point::get_rho()const {return this->rho;}
double Point::get_M()const {return this->M;}
int Point::get_id()const {return this->cid;}
/// TODO: get loc

bool Point::operator==(const Point& other)const
{
    return this->get_id() == other.get_id();
}
bool Point::operator!=(const Point&other)const
{
    return this->get_id() != other.get_id();
}
bool Point::operator>(const Point&other)const
{
    return this->get_id() > other.get_id();
}
bool Point::operator>=(const Point&other)const
{
    return this->get_id() >= other.get_id();
}
bool Point::operator<(const Point&other)const
{
    return this->get_id() < other.get_id();
}
bool Point::operator<=(const Point&other)const
{
    return this->get_id() <= other.get_id();
}
