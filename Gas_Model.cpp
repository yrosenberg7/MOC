#include "Gas_Model.h"
#include <cmath>
#include <iostream>
/** Initial Tests Passed */
// Default Constructor
Gas_Model::Gas_Model()
{
    this->gamma = 1.4;
    this->R = 287.04;
    this->T0 = 300; ///
    this->p0 = 30e5;

    this->a0 = sqrt(gamma*R*T0);
    this->Tstar = 2*T0/(gamma+1);
    this->astar = sqrt(gamma*R*Tstar);
}

// Specifying Gas and stagnation properties
Gas_Model::Gas_Model(double gamma, double R, double T0, double p0)
{
    this->gamma = gamma;
    this->R = R;
    this->T0 = T0;
    this->p0 = p0;

    this->a0 = sqrt(gamma*R*T0);
    this->Tstar = 2*T0/(gamma+1);
    this->astar = sqrt(gamma*R*Tstar);
}

void Gas_Model::print()
{
    std::cout << "GAS MODEL:  gamma: " << gamma << ",\t"
         << "R: " << R << ",\t"
         << "T0: " << T0 << ",\t"
         << "p0: " << p0 << ",\t"
         << "a0: " << a0 << ",\t"
         << "Tstar: " << Tstar << ",\t"
         << "astar: " << astar << std::endl;
}


double Gas_Model::a_from_u(double u){

        //std::cout << "is u here?: " << u << std::endl; /// yes it is! WHY IS IT NOT WORKING?!?!?!

        ///sqrt(abs(gas_model.a0^2 -0.5*(gas_model.gamma-1)*u.^2)
        double term1 = pow(this->get_a0(),2.0);
        double term2 = 0.5*(this->get_gamma()-1)*pow(u,2);

        double a = sqrt(fabs(term1-term2)); // a
//       return sqrt(pow(a0,2.0)-0.5*(gamma-1)*pow(u,2)); // a
        //std::cout << "is a here?: " << a << std::endl;
        //std::cout <<"\nu: "<< u << ", term1: " << term1 << ", term2: " << term2 << ", a: " << a <<std::endl;
        return a;
}

double Gas_Model::p_from_T(double T){
    //gas_model.p_from_T=@(T)(gas_model.p0*(T/gas_model.T0).^(gas_model.gamma/(gas_model.gamma-1)));

    //Returns p0 WHHHHYYYYY:(
    //return ((p0)*(pow(T/T0,gamma/gamma-1)));

    // THIS ONE WORKS
    //return this->get_p0()*pow(T/this->get_T0(),this->get_gamma()/(this->get_gamma()-1));
    //return get_p0()*pow(T/get_T0(),get_gamma()/(get_gamma()-1));
    return (p0*pow(T/T0,gamma/(gamma-1)));
}
double Gas_Model::T_from_a(double a){
    ///gas_model.T_from_a=@(a)(a^2/(gas_model.gamma*gas_model.R));

    return pow(a,2)/(gamma*R); //a
}
double Gas_Model::rho_from_pT(double p,double T){
    ///gas_model.rho_from_p_T=@(p, T)(p/(gas_model.R*T));
    return p/(R*T); //rho
}

/** Getter methods */
double Gas_Model::get_gamma(){return this->gamma;}
double Gas_Model::get_R(){return this->R;}
double Gas_Model::get_T0(){return this->T0;}
double Gas_Model::get_p0(){return this->p0;}
double Gas_Model::get_a0(){return this->a0;}
double Gas_Model::get_Tstar(){return this->Tstar;}
double Gas_Model::get_astar(){return this->astar;}
