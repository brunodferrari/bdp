/*
    Bipartite Graph Class
    Author: Bruno Ferrari
*/

#include <iostream>
#include <limits>
#include <vector>
#include "bgraph.h"

using namespace std;

Vertex::Vertex(){
    this->v = -1;
    this->pos = new int;
    (*this->pos) = 0;
}


Vertex::Vertex(int i, int p){
    this->v = i;
    this->pos = new int;
    (*this->pos) = p;
}

Vertex::~Vertex(){
    delete this->pos;
}

int Vertex::get_vertex() const {
    return v;
}

int Vertex::get_pos() const {
    return (*pos);
}

