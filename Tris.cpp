#include "Tris.h"

#include <stdio.h>
#include <iostream>

unsigned int Tris::id = 0;

Tris::Tris()
{
    this->pt1 = nullptr;
    this->pt2 = nullptr;
    this->pt3 = nullptr;
    this->cid = 0;
}

Tris::Tris(Point* p1, Point* p2, Point* p3)
{
    this->pt1 = p1;
    this->pt2 = p2;
    this->pt3 = p3;
    this->cid = id++;
}

Tris::Tris(const Tris& other)
{
    this->pt1 = other.get_p1();
    this->pt2 = other.get_p2();
    this->pt3 = other.get_p3();
    this->cid = other.get_cid();
}

Tris::~Tris(){}

void Tris::print()
{
    std::cout << (this->pt1)->get_id() << ", ";
    std::cout << (this->pt2)->get_id() << ", ";
    std::cout << (this->pt3)->get_id() << std::endl;
}

///down left trsi
Tris Tris::tri(int i, int j, std::vector<std::vector<Point> >& grid)
{
    /// requires indecis passed in are for the TOP LEFT of a triangle

    /**
    p1
    |  \
    p2 - p3
    */
    Point* p1 = &(grid[i][j]);
    Point* p2 = &(grid[i][j+1]);
    Point* p3 = &(grid[i+1][j]);

    Tris tris = Tris(p1,p2,p3);
    return tris;
}

/// up right tris
Tris Tris::tri2(int i, int j, std::vector<std::vector<Point> >& grid)
{
    /// requires indecis passed in are for the TOP LEFT of a triangle

    /**
     p1 - p3
       \  |
          p2
    */
    Point* p4 = &(grid[i][j]);
    Point* p5 = &(grid[i+1][j]);
    Point* p6 = &(grid[i+1][j-1]);

    Tris tris = Tris(p4,p5,p6);
    return tris;

}
