#include "MoC_Solver.h"
//#include <cmath> // sqrt
#include <math.h>
#include <vector>

#include "Point.h"
#include "Gas_Model.h"
#include "Throat.h"

#include <Eigen/Dense> // A\b

MoC_Solver::MoC_Solver(){}
///
MoC_Solver::MoC_Solver(Gas_Model* GM_in, Throat* TH_in, double eps_in, bool disp_in, int iter_in)
{
    this->GM = GM_in;
    this->TH = TH_in;
    this->eps = eps_in;
    this->disp = disp_in;
    this->iter = iter_in;
}
MoC_Solver::~MoC_Solver()
{
    //dtor
}

Point MoC_Solver::Symmetry_Pt(Point* pt2)
{
    double a0 = this->GM->get_a0();
    double gamma = this->GM->get_gamma();
    int delta = this->TH->get_delta();
    double eps = this->get_eps();

    double x2 = pt2->get_x();
    double y2 = pt2->get_y();
    double u2 = pt2->get_u();
    double v2 = pt2->get_v();
    double a2 = this->GM->a_from_u(sqrt(u2*u2 + v2*v2));

    double lam2_p = calc_lam_p(u2,v2,a2);
    double lam2_n = calc_lam_n(u2,v2,a2);

    Eigen::Matrix4d A;
    Eigen::Vector4d b;
    Eigen::Vector4d x;
    Eigen::Vector4d x_old;

    bool done = false;
    double err = 1;
    int counter = 0;
    double y23,u23,v23,lam23,x3,y3,u3,v3,a3,lam3_p,lam3_n,a23_2,Q23,R23,S23;

    while(!done)
    {
        if(counter == 0) /// set initial state
        {
            if(y2 != 0 ){y23 = y2;}
            else{y23 = 1;}

            u23 = u2;
            v23 = v2;
            lam23 = lam2_n; /// notes say positive but doesn't work?
        }
        else /// update guess
        {
            x_old = x;

            x3 = x_old(0);
            y3 = x_old(1);
            u3 = x_old(2);
            v3 = x_old(3);
            a3 = this->GM->a_from_u(sqrt(u3*u3 + v3*v3));
            y23 = mean(y2,y3);

            lam3_p = calc_lam_p(u3,v3,a3);
            lam3_n = calc_lam_n(u3,v3,a3);
            lam23 = mean(lam2_n,lam3_n);

            u23 = mean(u2,u3);
            v23 = mean(v2,v3);
        }

        a23_2 = a0*a0 - (gamma/2 - 0.5)*(u23*u23 + v23*v23);
        Q23 = u23*u23 - a23_2;
        R23 = 2*u23*v23 - Q23*lam23;
        S23 = delta*(a23_2*v23/y23);

        A << lam23,-1,0,0,
             0,1,0,0,
             -1*S23, 0, Q23, R23,
             0,0,0,1;
        b << lam23*x2-y2,
             0,
             -1*S23*x2+Q23*u2+R23*v2,
             0;

        /// X = x3,y3,u3,u4'
        x = A.colPivHouseholderQr().solve(b); /// A\b

        if(counter != 0) // skip otherwise
        {
            err = max_err(x, x_old);
            if(this->get_disp()){std::cout << "at counter " << counter << " error: "  << err << std::endl;}
            /// find err
            if (err < eps)
            {
                done = true;
            }else if(counter >= this->get_iter())
            {
                std::cout << "WARNING: MoC Irrotational Symmetry Point Solve did not converge in "<< this->get_iter() <<" iterations. " << std::endl;
                break;
            }
        }

        counter++;
    } /// convergence loop

    return Point(x(0),x(1),x(2),x(3));
}

Point MoC_Solver::Interior_Pt(Point* pt1,Point* pt2)
{
    double a0 = this->GM->get_a0();
    double gamma = this->GM->get_gamma();
    int delta = this->TH->get_delta();
    double eps = this->get_eps();

    double x1 = pt1->get_x();
    double y1 = pt1->get_y();
    double u1 = pt1->get_u();
    double v1 = pt1->get_v();
    double a1 = this->GM->a_from_u(sqrt(u1*u1 + v1*v1));

    double x2 = pt2->get_x();
    double y2 = pt2->get_y();
    double u2 = pt2->get_u();
    double v2 = pt2->get_v();
    double a2 = this->GM->a_from_u(sqrt(u2*u2 + v2*v2));

    double lam1_p = calc_lam_p(u1,v1,a1);
    double lam1_n = calc_lam_n(u1,v1,a1);
    double lam2_p = calc_lam_p(u2,v2,a2);
    double lam2_n = calc_lam_n(u2,v2,a2);

    Eigen::Matrix4d A;
    Eigen::Vector4d b;
    Eigen::Vector4d x;
    Eigen::Vector4d x_old;

    bool done = false;
    double err = 1;
    int counter = 0;
    double y23,u23,v23,lam23,x3,y3,u3,v3,a3,lam3_p,lam3_n,a23_2,Q23,R23,S23;
    double u13,v13,y13,lam13,a13_2,Q13,R13,S13;

    while(!done)
    {
        if(counter == 0) /// set initial state
        {
            /// ###
            u13 =  u1;
            v13 = v1;

            if(y1 != 0){y13 = y1;}
            else{y13 = 1;}

            if(y2 != 0){y23 = y2;}
            else{y23 = 1;}

            u23 = u2;
            v23 = v2;
            lam13 = lam1_p;
            lam23 = lam2_n; /// notes say positive but doesn't work?
            //lam23 = lam2_p; /// notes say positive but doesn't work?

        }
        else /// update guess
        {
            x_old = x;

            x3 = x_old(0);
            y3 = x_old(1);
            u3 = x_old(2);
            v3 = x_old(3);
            a3 = this->GM->a_from_u(sqrt(u3*u3 + v3*v3));

            y13 = mean(y1,y3);
            y23 = mean(y2,y3);

            lam3_p = calc_lam_p(u3,v3,a3);
            lam3_n = calc_lam_n(u3,v3,a3);
            lam13 = mean(lam1_p,lam3_p);
            lam23 = mean(lam2_n,lam3_n);

            u13 = mean(u1,u3);
            v13 = mean(v1,v3);
            u23 = mean(u2,u3);
            v23 = mean(v2,v3);
        }

        /// positive characteristic 13
//    a13_2 = a0^2 - (gamma/2 - 1/2)*(u13^2 + v13^2); % = a13^2
//    Q13 = u13^2 - a13_2;
//    R13 = 2*u13*v13 - Q13*lam13;
//    S13 = delta*(a13_2 * v13/y13);
        a13_2 = a0*a0 - (gamma/2 - 0.5)*(u13*u13 + v13*v13);
        Q13 = u13*u13 - a13_2;
        R13 = 2*u13*v13 - Q13*lam13;
        S13 = delta*(a13_2*v13/y13);

        /// check! negative characteristic 23
        a23_2 = a0*a0 - (gamma/2 - 0.5)*(u23*u23 + v23*v23);
        Q23 = u23*u23 - a23_2;
        R23 = 2*u23*v23 - Q23*lam23;
        S23 = delta*(a23_2*v23/y23);

        A << lam13,-1,0,0,
             lam23,-1,0,0,
             -1*S13, 0, Q13, R13,
             -1*S23, 0, Q23, R23;
        b << lam13*x1-y1,
             lam23*x2-y2,
             (-1*S13*x1)+(Q13*u1)+(R13*v1),
             (-1*S23*x2)+(Q23*u2)+(R23*v2);

        /// X = x3,y3,u3,u4'
        x = A.colPivHouseholderQr().solve(b); /// A\b

        if(counter != 0) // skip otherwise
        {
            err = max_err(x, x_old);
            if(this->get_disp()){std::cout << "at counter " << counter << " error: "  << err << std::endl;}
            /// find err
            if (err < eps)
            {
                done = true;
            }else if(counter >= this->get_iter())
            {
                std::cout << "WARNING: MoC Irrotational Symmetry Point Solve did not converge in "<< this->get_iter() <<" iterations. " << std::endl;
                break;
            }
        }

        counter++;
    } /// convergence loop

    return Point(x(0),x(1),x(2),x(3));

}

Point MoC_Solver::ISEN_Interior_Pt(Point* pt1, Point* pt2)
{
    ///make sure a,p,T,rho,M have been calculated
    pt1->finish_pt(this->GM);
    pt2->finish_pt(this->GM);
    double a0 = this->GM->get_a0();
    double gamma = this->GM->get_gamma();
    int delta = this->TH->get_delta();
    double eps = this->get_eps();

    double x1 = pt1->get_x();
    double y1 = pt1->get_y();
    double u1 = pt1->get_u();
    double v1 = pt1->get_v();
    double a1 = pt1->get_a();
    double p1 = pt1->get_p();
    double T1 = pt1->get_T();
    double rh1 = pt1->get_rho();
    double M1 = pt1->get_M();

    double x2 = pt2->get_x();
    double y2 = pt2->get_y();
    double u2 = pt2->get_u();
    double v2 = pt2->get_v();
    double a2 = pt2->get_a();
    double p2 = pt2->get_p();
    double T2 = pt2->get_T();
    double rh2 = pt2->get_rho();
    double M2 = pt2->get_M();

    double V2 = sqrt(u2*u2 + v2*v2);
    double V1 = sqrt(u1*u1 + v1*v1);
    double th2 = atan2(v2,u2);
    double th1 = atan2(v1,u1);
    double alpha2 = asin(1/M2);
    double alpha1 = asin(1/M1);



    /// CREATE INTERPOLATION VARIABLES
    Eigen::Matrix2d A_interp;

    A_interp << x1, y1, /// same for all
                x2, y2;

    Eigen::Vector2d x_interpv;
    x_interpv = Interp_Vars(A_interp,V1,V2);
    double av = x_interpv(0);
    double bv = x_interpv(1);

    Eigen::Vector2d x_interpth ;
    x_interpth = Interp_Vars(A_interp,th1,th2);
    double ath = x_interpth(0);
    double bth = x_interpth(1);

    Eigen::Vector2d x_interpp;
    x_interpp = Interp_Vars(A_interp,p1,p2);
    double ap = x_interpp(0);
    double bp = x_interpp(1);

    Eigen::Vector2d x_interprh;
    x_interprh = Interp_Vars(A_interp,rh1,rh2);
    double arh = x_interprh(0);
    double brh = x_interprh(1);


    /// Solver matrices
    Eigen::MatrixXd A(12,12); /// 12 vars!
    Eigen::VectorXd b(12);
    Eigen::VectorXd x(12);
    Eigen::VectorXd x_old(12);

    bool done = false;
    double err = 1;
    int counter = 0;

    double x3,y3,V3,p3,th3,rh3,x4,y4,p4,V4,th4,rh4;
    double lam_0, lam_12, P_0, A_0, R_0;
    double lam_p, lam_n, Q_p, Q_n,S_p,S_n;
    double a3,alpha3,a4,T4;

    while(!done)
    {
        if(counter == 0) /// set initial state
        {
            lam_p = tan(th2 + alpha2);
            lam_n = tan(th1 + alpha1);

            Q_p = sqrt(M2*M2 - 1)/(rh2*V2*V2);
            Q_n = sqrt(M1*M1 - 1)/(rh1*V1*V1);
            S_p = delta*sin(th2)/(y2*M2*cos(th2+alpha2));
            S_n = delta*sin(th1)/(y1*M1*cos(th1+alpha1));

            /// point 3 values: avg of 1,2
            p3 = mean(p1,p2);
            a3 = mean(a1,a2);
            rh3 = mean(rh1,rh2);
            V3 = mean(V1,V2);
            th3 = mean(th1,th2);
            alpha3 = mean(alpha1,alpha2);

            /// need p4 to solve for
            a4 = a2;
            V4 = V1;
            T4 = this->GM->T_from_a(a4);
            p4 = this->GM->p_from_T(T4);
            rh4 = this->GM->rho_from_pT(p4,T4);

        }
        else /// update guess
        {
            x_old = x;

            x3 = x_old(0);
            y3 = x_old(1);
            V3 = x_old(2);
            th3 = x_old(3);
            p3 = x_old(4);
            rh3 = x_old(5);
            x4 = x_old(6);
            y4 = x_old(7);
            V4 = x_old(8);
            th4 = x_old(9);
            p4 = x_old(10);
            rh4 = x_old(11);

            /// UPDATE THESE: AVG OF 1or2 and 4
            lam_p = tan(th2 + alpha2);
            lam_n = tan(th1 + alpha1);
            Q_p = sqrt(M2*M2 - 1)/(rh2*V2*V2);
            Q_n = sqrt(M1*M1 - 1)/(rh1*V1*V1);
            S_p = delta*sin(th2)/(y2*M2*cos(th2+alpha2));
            S_n = delta*sin(th1)/(y1*M1*cos(th1+alpha1));

            a3 = this->GM->a_from_u(V3);
            alpha3 = asin(a3/V3);

        }

        R_0 = rh3*V3;
        A_0 = a3*a3*a3;
        lam_0 = tan(th3+alpha3);

        A << -1*lam_0, 1, 0, 0, 0, 0, lam_0, -1, 0, 0, 0, 0,
             0,       0, 0, 0, 0, 0, lam_p, 1, 0, 0, 0, 0,
             0,       0, 0, 0, 0, 0, -1*lam_n, 1, 0, 0, 0, 0,
             -1*lam_0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             av, bv, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             ath, bth, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
             ap, bp, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0,
             arh, brh, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, S_p, 0, 0, 1, Q_p, 0,
             0, 0, 0, 0, 0, 0, S_n, 0, 0, -1, Q_n, 0,
             0, 0, 0, 0, 0, 0, 0, 0, p3,  0,   1, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0,  0,   1, -1*A_0;

        b << 0,
            -lam_p*x2 + y2,
            -lam_n*x1 + y1,
            -lam_12*x2 + y1,
            0,
            0,
            0,
            0,
            S_p*x2 + Q_p*p2 + th2,
            S_n*x1 + Q_n*p1 + th1,
            R_0*V3 + p3,
            p3 - A_0*rh4;

        /// X = x3,y3,V3,Th3,p3,rh3,x4,y4,V4,Th4,p4,rh4
        x = A.colPivHouseholderQr().solve(b); /// A\b

        if(counter != 0) // skip otherwise
        {
            err = max_err(x, x_old);
            if(this->get_disp()){std::cout << "at counter " << counter << " error: "  << err << std::endl;}
            /// find err
            if (err < eps)
            {
                done = true;
            }else if(counter >= this->get_iter())
            {
                std::cout << "WARNING: MoC Irrotational Symmetry Point Solve did not converge in "<< this->get_iter() <<" iterations. " << std::endl;
                break;
            }
        }

        counter++;
    } /// convergence loop

    return Point();
}



Point MoC_Solver::Fixed_Wall_Pt(Point* pt1, double (*fw)(double), double (*fp)(double))
{
    double a0 = this->GM->get_a0();
    double gamma = this->GM->get_gamma();
    int delta = this->TH->get_delta();
    double eps = this->get_eps();

    double x1 = pt1->get_x();
    double y1 = pt1->get_y();
    double u1 = pt1->get_u();
    double v1 = pt1->get_v();
    double a1 = this->GM->a_from_u(sqrt(u1*u1 + v1*v1));

    double lam1_p = calc_lam_p(u1,v1,a1);
    double lam1_n = calc_lam_n(u1,v1,a1);

    Eigen::Matrix4d A;
    Eigen::Vector4d b;
    Eigen::Vector4d x;
    Eigen::Vector4d x_old;

    bool done = false;
    double err = 1;
    int counter = 0;
    double y13,u13,v13,lam13,x3,y3,u3,v3,a3,lam3_p,lam3_n,a13_2,Q13,R13,S13,alpha;

    while(!done)
    {
        if(counter == 0) /// set initial state
        {
            u13 = u1;
            v13 = v1;
            lam13 = lam1_p; /// notes say positive but doesn't work?
            x3 = x1;
            //y3 = this->TH->fwall(x3);
            y3 = (*fw)(x3);
            y13 = mean(y1,y3);
        }
        else /// update guess
        {
            x_old = x;

            x3 = x_old(0);
            //y3 = this->TH->fwall(x3);
            y3 = (*fw)(x3);
            u3 = x_old(2);
            v3 = x_old(3);
            a3 = this->GM->a_from_u(sqrt(u3*u3 + v3*v3));
            y13 = mean(y1,y3);

            lam3_p = calc_lam_p(u3,v3,a3);
            lam3_n = calc_lam_n(u3,v3,a3);
            lam13 = mean(lam1_p,lam3_p);

            u13 = mean(u1,u3);
            v13 = mean(v1,v3);
        }

        a13_2 = a0*a0 - (gamma/2 - 0.5)*(u13*u13 + v13*v13);
        Q13 = u13*u13 - a13_2;
        R13 = 2*u13*v13 - Q13*lam13;
        S13 = delta*(a13_2*v13/y13);
        //alpha = this->TH->fpwall(x3);
        alpha = (*fp)(x3);

        A << lam13,-1,0,0,
             alpha,-1,0,0,
             -1*S13, 0, Q13, R13,
             0,0,alpha,-1;
        b << lam13*x1-y1,
             alpha*x3-y3,
             -1*S13*x1+Q13*u1+R13*v1,
             0;

        /// X = x3,y3,u3,u4'
        x = A.colPivHouseholderQr().solve(b); /// A\b

        if(counter != 0) // skip otherwise
        {
            err = max_err(x, x_old);
            if(this->get_disp()){std::cout << "at counter " << counter << " error: "  << err << std::endl;}
            /// find err
            if (err < eps)
            {
                done = true;
            }else if(counter >= this->get_iter())
            {
                std::cout << "WARNING: MoC Irrotational Symmetry Point Solve did not converge in "<< this->get_iter() <<" iterations. " << std::endl;
                break;
            }
        }

        counter++;
    } /// convergence loop

    return Point(x(0),x(1),x(2),x(3));
}
//Point MoC_Solver::Variable_Wall_Pt(Point*,Point*);

double MoC_Solver::calc_lam_p(double u, double v, double a)
{
    return (u*v + a*sqrt(u*u + v*v - a*a))/(u*u - a*a);
}
double MoC_Solver::calc_lam_n(double u, double v, double a)
{
    return (u*v - a*sqrt(u*u + v*v - a*a))/(u*u - a*a);
}

double MoC_Solver::mean(double a, double b)
{
    return (a+b)/2;
}

double MoC_Solver::max_err(Eigen::Vector4d x, Eigen::Vector4d x_old) /// assumes the length of the vector :(
{
    if(this->get_disp())
    {
        std::cout << std::endl;
        print_Vec(x_old);
        print_Vec(x);
    }

    int n = x.size();
    double vars[n];
    long double tmp;

    for(int i=0; i<n; i++)
    {
        tmp =  x(i)-x_old(i);
        if(this->get_disp()){std::cout << "comparison "<< i << ": " << tmp << std::endl;}
        vars[i] = fabs(tmp);

    }

    return *std::max_element(vars,vars+n);
}

void MoC_Solver::print_Vec(Eigen::Vector4d x)
{
    std::cout << "x(0): " << x(0) << ", x(1): " << x(1) << ", x(2): " << x(2) << ", x(3): " << x(3) << std::endl;
}

Eigen::Vector2d MoC_Solver::Interp_Vars(Eigen::Matrix2d A, double var1, double var2)
{
    Eigen::Vector2d b_interp;
    b_interp << var1,var2;
    Eigen::Vector2d x;
    x = A.colPivHouseholderQr().solve(b_interp); /// A\b
    return x;
}
