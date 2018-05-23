#ifndef THROAT_SOLVER_H
#define THROAT_SOLVER_H

#include <list>


#include "Point.h"
#include "Gas_Model.h"
#include "Throat.h"


class Throat_Solver
{
    public:
        Throat_Solver();
        Throat_Solver(Gas_Model*,Throat*);
        virtual ~Throat_Solver();

        std::list<Point> Compute_IDL(int npts);
//        std::array<Point> Compute_IDL(int npts);
        std::list<std::list<Point>> Compute_THROAT(std::list<Point>* idl);
//        std::list<std::list<Point>> Compute_THROAT(std::array<Point>* idl);


    protected:
        Point Transonic_Velocity(double x, double y);

    private:
        Gas_Model* GM;
        Throat* TH;

        double astar;
        double alpha;
        double eps;
        //double xbar;

        double get_astar();
        double get_alpha();
        double get_eps();
        //double get_xbar();
};
#endif // THROAT_SOLVER_H


//std::list<Point> Compute_IDL(Gas_Model* GM,Throat* TH,int npts);
//std::list<Point> Compute_THROAT(Gas_Model* GM,Throat* TH,std::list<Point>);
//Point Compute_Throat((Gas_Model* GM,Throat*, double, double);
