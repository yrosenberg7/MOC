#ifndef TRIS_H
#define TRIS_H

#include "Point.h"
#include <list>
#include <vector>
class Tris
{
    public:
        Tris();
        Tris(Point*,Point*,Point*);
        Tris(const Tris& other);
        virtual ~Tris();

        Tris tri(int i, int j, std::vector<std::vector<Point> >& grid);
        Tris tri2(int i, int j, std::vector<std::vector<Point> >& grid);

        void print(); //prints to std out

        Point* get_p1()const{return this->pt1;}
        Point* get_p2()const{return this->pt2;}
        Point* get_p3()const{return this->pt3;}
        int get_cid()const{return this->cid;}


    protected:
        Point* pt1;
        Point* pt2;
        Point* pt3;
        unsigned int cid;
        static unsigned int id;


    private:
};

#endif // TRIS_H
