#ifndef MOC_SOLVER_H
#define MOC_SOLVER_H

#include "Point.h"
#include "Gas_Model.h"
#include "Throat.h"

class MoC_Solver
{
    public:
        MoC_Solver();
        Moc_Solver(Gas_Model*, Throat*);
        virtual ~MoC_Solver();

        Point Symmetry_Pt(Point*);
        Point Interior_Pt(Point*,Point*);
        Point Fixed_Wall_Pt(Point*,Point*);
        Point Variable_Wall_Pt(Point*,Point*);

    protected:
        double calc_lam_p(double, double, double);
        double calc_lam_n(double, double, double);
    private:
        Gas_Model* GM;
        Throat* TH;
};

#endif // MOC_SOLVER_H
