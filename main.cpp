#include <iostream>
#include <cmath>
#include <list>
#include <array>
#include <fstream>
#include <iterator>
#include <stdio.h>
#include <sstream>

#include "Gas_Model.h"
#include "Throat.h"
#include "Point.h"
#include "Throat_Solver.h"

//using namespace std;


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
///ver 0.0.2 throat solve <- HERE...
///ver 0.0.3 throat tris
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
    std::cout << "\n" << "MoC Supersonic Nozzle Analysis\n" << "C++ Implementation\n"<< "Written by Yoav Rosenberg\n" << "ver 0.0.1" << "\n" << std::endl;

    ///declare inputs
    int npts;
    double upstream_rth,downstream_rth,height,max_angle,gamma,R,T0,p0;
    bool axi,ysmooth;
    std::string inputs[12];

    //FILE* input_file;
    //input_file = fopen("MoC_input2.txt","r");

    /// file io
    FILE* output_file;
    output_file = fopen("MoC_output.txt","w");
    std::ifstream input_file("MoC_input2.txt");

    ///extract inputs
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
    std::list<Point> idl = TS.Compute_IDL(npts); /// compute IDL
//    array<Point> idl = TS.Compute_IDL(npts); /// compute IDL

    /// ### COMPUTE THROAT REGION ###
    std::list<std::list<Point>> Throat_pts = TS.Compute_THROAT(&idl); /// compute throat region

    ///point values format / order
    std::cout << "\n" << "id" <<"\t"  << "x" <<"\t"
                      << "y" <<"\t" << "u" <<"\t"
                      << "v" <<"\t" << "a" <<"\t"
                      << "T" <<"\t"  << "p" <<"\t"
                      << "rho" <<"\t" << "M" <<"\t" << std::endl;
    /// finish pts, print
    for (std::list<Point>::iterator it=idl.begin(); it != idl.end(); ++it)
        {
            it->finish_pt(&GM);
//            it->print_xm();
            it->print();
        }

    ///close output file
    fclose(output_file);


    std::cout << "\n\n###### END OF PROGRAM ######" << std::endl;
    return 0;
}
