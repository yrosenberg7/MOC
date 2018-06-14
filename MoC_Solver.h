#ifndef MOC_SOLVER_H
#define MOC_SOLVER_H

#include <vector>
#include "Point.h"
#include "Gas_Model.h"
#include "Throat.h"
#include <Eigen/Dense>
class MoC_Solver
{
    public:
        MoC_Solver();
        MoC_Solver(Gas_Model*, Throat*, double,bool,int);
        virtual ~MoC_Solver();

        Point Symmetry_Pt(Point*);
        Point Interior_Pt(Point*,Point*);
        Point Fixed_Wall_Pt(Point*);
        Point Fixed_Wall_Pt(Point* pt1, double (*fw)(double), double (*fp)(double));
        Point Variable_Wall_Pt(Point*,Point*);

        double get_eps()const{return this->eps;}
        bool get_disp() const{return this->disp;}
        int get_iter() const{return this->iter;}

    protected:
        double calc_lam_p(double, double, double);
        double calc_lam_n(double, double, double);
        //double get_max(std::vector<double> X);
        double mean(double a, double b);
        double max_err(Eigen::Vector4d x, Eigen::Vector4d x_old);
        void print_Vec(Eigen::Vector4d x);
    private:
        Gas_Model* GM;
        Throat* TH;
        double eps;
        bool disp;
        int iter;
};

#endif // MOC_SOLVER_H
