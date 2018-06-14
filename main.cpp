#include <iostream>
#include <cmath>
#include <list>
#include <array>
#include <fstream>
#include <iterator>
#include <stdio.h>
#include <sstream>
#include <vector>

#include "Gas_Model.h"
#include "Throat.h"
#include "Point.h"
#include "Throat_Solver.h"
#include "Tris.h"


/// TESTING GITHUB VERSION CONTROL, DID THIS UPLOAD??????

///TODO: ***complete throat region solve***
///         - add astar, alpha, eps to throat solver constructor, with reference to GM and TH
///         - create helper func for IDL solve
///         - change idl solve to use updated class
///         - change transonic throat solve to use updated class
///         - change THROAT SOLVE to use list of points
///         - make debug version of throat solve for only x,y
///         - test transonic throat solve for x,y locations
///         - add copy constructor to point class, (copy id!)
///         - output is list of throat points, including idl
///     *** throat mesh ***
///         - create tris class, with print to file(s)
///         -

///      Add location tags to points



/** INPUT FILE FORMAT
npts	11
upstream_rth	2.0
downstream_rth	0.5
height	1
max_angle	13.65
axi     1
ymooth  1
gamma	1.4
R	287.04
T0	300
p0	30e5
*/

///ver 0.0.1 idl
///ver 0.0.2 throat solve
///ver 0.0.3 throat tris <- HERE...
///ver 0.1.1 throat display
///ver 0.2.2 idl R
///ver 0.2.3 idl R tris
///ver 0.2.4 idl R display
///ver 0.3.1 fixed wall R
///ver 0.3.2 fixed wall R tris
///ver 0.3.3 variable wall
///ver 0.3.4 variable wall tris

///ver 1.0.0 AERO 400 project complete!

int main()
{
    ///header
    std::cout << "\n" << "MoC Supersonic Nozzle Analysis\n" << "C++ Implementation\n"<< "Written by Yoav Rosenberg\n" << "ver 0.0.2" << "\n" << std::endl;

    ///declare inputs
    int npts;
    double upstream_rth,downstream_rth,height,max_angle,gamma,R,T0,p0;
    bool axi,ysmooth;

    std::list<Point> POINTS;

    std::string input_file_name = "MoC_input2.txt";
    std::string output_pts_file_name = "MoC_Points.txt";
    std::string output_tris_file_name = "MoC_Tris.txt";
    std::string inputs[12];

    //FILE* input_file;
    //input_file = fopen("MoC_input2.txt","r");

    /// file io
    FILE* output_file_pts;
    FILE*output_file_tris;
    output_file_pts = fopen(output_pts_file_name.c_str(),"w");
    output_file_tris = fopen(output_tris_file_name.c_str(),"w");
    std::ifstream input_file(input_file_name);

    ///extract inputs
    std::cout<< "Extracting inputs from file: "  << input_file_name << std::endl;
    std::cout<< "Output point info to file: "  << output_pts_file_name << std::endl;
    std::cout<< "Output Tris to file: "  << output_tris_file_name << std::endl << std::endl;
    if(input_file.is_open())
    {
        for(int i = 0; i < 11; ++i)
        {
            input_file >> inputs[i];
        }
    }
    input_file.close();

    ///extract / cast inputs
    std::istringstream ( inputs[0] ) >> npts; /// idl pts
    std::istringstream ( inputs[1] ) >> upstream_rth; /// throat geom
    std::istringstream ( inputs[2] ) >> downstream_rth;
    std::istringstream ( inputs[3] ) >> height;
    std::istringstream ( inputs[4] ) >> max_angle;
    std::istringstream ( inputs[5] ) >> axi;
    std::istringstream ( inputs[6] ) >> ysmooth;
    std::istringstream ( inputs[7] ) >> gamma; /// gas model
    std::istringstream ( inputs[8] ) >> R;
    std::istringstream ( inputs[9] ) >> T0;
    std::istringstream ( inputs[10] ) >> p0;

    ///instantiate models
    Throat TH = Throat(upstream_rth,downstream_rth,height,max_angle,axi,ysmooth);//Throat TH = Throat();
    Gas_Model GM = Gas_Model(gamma,R,T0,p0);// GM = Gas_Model();

    ///display models
    GM.print();
    TH.print();

    Throat_Solver TS = Throat_Solver(&GM, &TH);
    std::cout << "\nComputing IDL..." << std::endl;
    std::vector<Point> idl = TS.Compute_IDL(npts); /// compute IDL

    /// $$$ DISPLAT IDL $$$
    for(int i = 0; i < npts; i++)
    {
        idl[i].print();
        POINTS.emplace_back(idl[i]); /// copy
    }


    /// ### COMPUTE THROAT REGION ###
    std::cout << "\nComputing Throat Region..." << std::endl;
    std::vector<std::vector<Point>> throat_pts(npts);
    std::vector<std::vector<Point> >::iterator il = throat_pts.begin();
    /// throat region looks like this:
    /// with * = idl pt, and . = transonic solve
    /// *
    /// . *
    /// . . *
    /// . . . *

    Point pt;
    for(int i = 0; i <npts; i++)
    {
        for(int j = i; j<npts; j++)
        {
            if(i==j)
            {
                il->emplace_back(idl[i]); /// copy point from the idl
                continue;
            }

            double x = idl[i].get_x();
            double y = idl[j].get_y();
            pt = TS.Transonic_Velocity(x,y);
            //pt = Point(x,y); /// ### solve at throat point ###
            il->emplace_back(pt);
            POINTS.emplace_back(pt);
        }
        il++;
    }

    /// $$$ DISPLAT IDL $$$
    std::vector<std::vector<Point> >::iterator print_it = throat_pts.begin();
    std::vector<Point>::iterator print_il;
    for(;print_it != throat_pts.end();print_it++)
    {
        print_il = print_it->begin();
        for(; print_il != print_it->end(); print_il++)
        {
            print_il->print();
        }
    }

/// ### COMPUT TRIS ###
    std::cout<< std::endl << std::endl;
    Tris tr = Tris(); /// used to call tri1 and tri2 funcs... it really should be a static member

    /// create tris from throat region
    std::list<Tris> TRIS;
    for(int i = 0; i < npts-1; i++)
    {
        for(int j = 0; j < npts-i-1; j++)
        {
            if(j == 0)
            {
                Tris tri1 = tr.tri(i,j,throat_pts);
                TRIS.push_back(tri1);
            }else
            {
                Tris tri1 = tr.tri(i,j,throat_pts);
                Tris tri2 = tr.tri2(i,j,throat_pts);
                TRIS.push_back(tri1);
                TRIS.push_back(tri2);
            }
        }
    }
    /// compute tris from idl_r

        ///point values format / order
    std::cout << "\nDisplaing Points: ";
    std::cout << "\n" << "id" <<"\t"  << "x" <<"\t"
                      << "y" <<"\t" << "u" <<"\t"
                      << "v" <<"\t" << "a" <<"\t"
                      << "T" <<"\t"  << "p" <<"\t"
                      << "rho" <<"\t" << "M" <<"\t" << std::endl;




    ///print throat
//    il = throat_pts.begin();
//    std::vector<Point>::iterator ip;
//    for(; il != throat_pts.end(); il++)
//    {
//        for(ip = il->begin(); ip != il->end(); ip++)
//        {
//            ip->print();
//        }
//    }

    /// finish pts, sort , print ### IMPLEMENT SORT ### /// ### Change to full point list ###
    POINTS.sort();
    for (std::list<Point>::iterator pts =POINTS.begin(); pts != POINTS.end(); ++pts)
    {
        pts->finish_pt(&GM);
        pts->print();
    }

    /// print tris
//    std::cout<< "printing tris: " << std::endl;
//    for(std::list<Tris>::iterator it = TRIS.begin(); it != TRIS.end(); it++)
//    {
//        it->print();
//    }

    ///close output file
    fclose(output_file_pts);
    fclose(output_file_tris);
    std::cout << "\n\n###### END OF PROGRAM ######" << std::endl;
    return 0;
}
