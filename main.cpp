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

int main()
{
    ///header
    std::cout << "MoC Supersonic Nozzle Analysis\n" << "C++ Implementation\n"<< "Written by Yoav Rosenberg\n" << std::endl;

    ///declare inputs
    int npts;
    double upstream_rth,downstream_rth,height,max_angle,gamma,R,T0,p0;
    bool axi,ysmooth;

    std::list<Point> POINTS;

    std::string input_file_name = "MoC_input.txt";
    std::string output_pts_file_name = "MoC_Points.txt";
    std::string output_tris_file_name = "MoC_Tris.txt";
    std::string inputs[12];

    std::cout<< "Extracting inputs from file: "  << input_file_name << std::endl;
    std::cout<< "Output point info to file: "  << output_pts_file_name << std::endl;
    std::cout<< "Output Tris to file: "  << output_tris_file_name << std::endl << std::endl;

    /// file io
    std::ofstream output_file_pts(output_pts_file_name);
    std::ofstream output_file_tris(output_tris_file_name);
    std::ifstream input_file(input_file_name);

    ///extract inputs
    if(input_file.is_open())
    {
        for(int i = 0; i < 11; ++i)
        {
            input_file >> inputs[i];
        }
    }
    else
    {
        throw "Failed to open input file!";
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
    std::cout << "Displaying analysis model: " << std::endl;
    Throat TH = Throat(upstream_rth,downstream_rth,height,max_angle,axi,ysmooth);//Throat TH = Throat();
    Gas_Model GM = Gas_Model(gamma,R,T0,p0);// GM = Gas_Model();

    ///display models
    GM.print();
    TH.print();

    Throat_Solver TS = Throat_Solver(&GM, &TH);
    std::cout << "\nComputing IDL..." << std::endl;
    std::vector<Point> idl = TS.Compute_IDL(npts); /// compute IDL

    for(int i = 0; i < npts; i++)
    {
        POINTS.emplace_back(idl[i]); /// copy to full points list
    }

    /// ### COMPUTE THROAT REGION ###
    std::cout << "Computing Throat Region..." << std::endl;
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
            il->emplace_back(pt);
            POINTS.emplace_back(pt);
        }
        il++;
    }

    /// ### COMPUT TRIS ###
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

    /// finish pts, sort , print
    POINTS.sort();
    for (std::list<Point>::iterator pts =POINTS.begin(); pts != POINTS.end(); ++pts)
    {
        pts->finish_pt(&GM);
        pts->print(&output_file_pts);
    }

    /// print tris
    for(std::list<Tris>::iterator it = TRIS.begin(); it != TRIS.end(); it++)
    {
        it->print(&output_file_tris);
    }

    ///close output file
    output_file_pts.close();
    output_file_tris.close();
    std::cout << "\n###### END OF PROGRAM ######" << std::endl;
    return 0;
}
