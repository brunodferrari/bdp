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
        //int *pos;

    public:
        int pos;
        Vertex();
        Vertex(int i, int p);
        ~Vertex();

        void setVertex(int i, int p);
        void setPos(int p);

        int get_vertex() const;
        int get_pos() const;

};

typedef unordered_map<int, unordered_map<int,int>> AdjMatrix;
typedef unordered_map<int, vector<int> > Adj;
typedef unordered_map<int, Vertex> Dict;
typedef map<int, Vertex> VSorted;

typedef struct Edge{
    int u;
    int v;
} edge;

typedef priority_queue<pair<Dict,int>, vector<pair<Dict,int>>, decltype(mcmp)> MinQueue;

class BGraph{
    private:
        int n_edges;
        Vertex *v1;  //list -> v in V_1  dim:  n_v1x1
        Vertex *v2;  //list -> v in V_2   dim:  n_v2x1
        Edge *edgelist;
        AdjMatrix edges;

    public:
        int n_v1;
        int n_v2;
        VSorted map_v1;
        VSorted map_v2;
        Dict pi_1;
        Dict pi_2;
        Adj adj_list;
        unordered_map<int,int> layer;
        int maped;

        BGraph(int n_v1_, int n_v2_, int n_edges_, Vertex *v1_, Vertex *v2_, Edge *e);

        unordered_map<int, int> degrees();

        unordered_map<int, int> degrees(unordered_set<int> V, unordered_set<int> V_sub);

        int greedy_selection(float alpha, unordered_set<int> V, unordered_set<int> V_sub, int seed = -1);

        float bc(int v, int k);

        int n_cross();

        void move_vertex(int v, int to);

        void move_vertex(int v, int to; int inplace);

        void re_map(int k);

        void printBGraph();

        void printGraph();

        void printMatrix()

        ~BGraph();
};


#endif // BGRAPH_H
