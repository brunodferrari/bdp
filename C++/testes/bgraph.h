/*
    Bipartite Graph Class
    Author: Bruno Ferrari
*/

#ifndef BGRAPH_H
#define BGRAPH_H

#define pi(v) v.get_pos()

class Vertex{

    private:
        int v;
        int *pos;

    public:
        Vertex(int i, int p);
        ~Vertex();

        int get_vertex() const;
        int get_pos() const;
};


class BGraph{
    private:
        int n_v1;
        int n_v2;
        int n_edges;
        Vertex *v1;
        Vertex *v2;


    public:
        BGraph(int *v1, int *v2);

        ~BGraph();


};


#endif // BGRAPH_H