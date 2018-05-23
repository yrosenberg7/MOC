#ifndef GAS_MODEL_H
#define GAS_MODEL_H


#include <iostream>
class Gas_Model
{
    public:

        Gas_Model();
        Gas_Model(double, double, double, double);
//        virtual ~Gas_Model(){std::cout << "Gas Model deconstructed" << std::endl;}
        virtual ~Gas_Model(){};

        void print();

        double a_from_u(double);
        double p_from_T(double);
        double T_from_a(double);
        double rho_from_pT(double,double);

        double get_gamma();
        double get_R();
        double get_T0();
        double get_p0();
        double get_a0();
        double get_Tstar();
        double get_astar();

    protected:

    private:
        double gamma;
        double R;
        double T0;
        double p0;
        double a0; // v- must be solved for -v
        double Tstar;
        double astar;
};

#endif // GAS_MODEL_H
