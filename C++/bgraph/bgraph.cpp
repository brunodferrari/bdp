/*
    Bipartite Graph Class
    Author: Bruno Ferrari
*/

#include "bgraph.h"

using namespace std;

Vertex::Vertex(){
    this->v = -1;
    this->pos = 0;
}

Vertex::Vertex(int i, int p){
    this->v = i;
    this->pos = p;
}

Vertex::~Vertex(){
}

void Vertex::setVertex(int i, int p){
    this->v = i;
    this->pos = p;
}

void Vertex::setPos(int p){
    //(*this->pos) = p;
    this->pos = p;
}

int Vertex::get_vertex() const {
    return v;
}

int Vertex::get_pos() const {
    return pos;
}

BGraph::BGraph(int n_v1_, int n_v2_, int n_edges_, Vertex *v1_, Vertex *v2_, Edge *e){
    this->n_v1 = n_v1_;
    this->n_v2 = n_v2_;
    this->n_edges = n_edges_;

    this->v1 = v1_;
    this->v2 = v2_;
    this->edgelist = e;

    int i,j,m = 0;
    int n = n_v1_ + n_v2_;

    for(m = 0; m < n_edges_; m++){
        i = e[m].u;
        j = e[m].v;
        this->edges[i][j] = 1;
        this->edges[j][i] = 1;

        this->adj_list[i].push_back(j);
        this->adj_list[j].push_back(i);

        this->edgelist[m].u = i;
        this->edgelist[m].v = j;
    }
    for (i=1;i<=n_v1_;i++){
        this->map_v1[v1_[i-1].get_pos()] = v1_[i-1];
        this->pi_1[v1_[i-1].get_vertex()] = v1_[i-1];
        this->layer[v1_[i-1].get_vertex()] = 1;
    }
    for (i=1;i<=n_v2_;i++){
        this->map_v2[v2_[i-1].get_pos()] = v2_[i-1];
        this->pi_2[v2_[i-1].get_vertex()] = v2_[i-1];
        this->layer[v2_[i-1].get_vertex()] = 2;
    }
    this->maped = 0;
}

unordered_map<int, int> BGraph::degrees(){
    int i,j;
    unordered_map<int, int> deg;
    for(int m = 0; m < n_edges; m++){
        i = this->edgelist[m].u;// - shift;
        j = this->edgelist[m].v;// - shift;

        deg[i]++;
        deg[j]++;
    }
    return deg;
}

unordered_map<int, int> BGraph::degrees(unordered_set<int> V, unordered_set<int> V_sub){
    //int i,j;
    unordered_map<int, int> deg;

    for(auto i = V.begin(); i != V.end(); i++){
        deg[*i] = 0;
        for(auto j = V_sub.begin(); j != V_sub.end(); j++){
            if (edges[*i][*j]) { deg[*i]++;}
        }
    }
    return deg;
}

int BGraph::greedy_selection(float alpha, unordered_set<int> V, unordered_set<int> V_sub, int seed = -1){
    unordered_map<int, int> deg_dict;
    if(!(V.size() && V_sub.size())){
        deg_dict = degrees();
    } else {
        deg_dict = degrees(V, V_sub);
    }

    //Max Heap
    auto cmp = [](pair<int,int> left, pair<int,int> right) { return left.second < right.second; }; // Funcao de Comparacao
    priority_queue<pair<int,int>, vector<pair<int,int>>, decltype(cmp)> max_q(cmp);
    for(auto v = deg_dict.begin(); v != deg_dict.end(); v++)
        max_q.push(*v);

    int deg_max = max_q.top().second * alpha;

    int max_idx = 0;
    vector<int> rcl;

    while( max_q.top().second >= deg_max && !max_q.empty()){
        rcl.push_back(max_q.top().first);
        max_idx++;
        max_q.pop();
    }
    if (seed!=-1) { srand(seed); }
    int idx = rand() % max_idx;

    return rcl[idx];
}

float BGraph::bc(int v, int k){ //barycenter
    float b=0;
    int K=0;
    Dict *pi_k;
    switch(k){
        case 1:
            pi_k = &this->pi_2;
            break;
        case 2:
            pi_k = &this->pi_1;
            break;
    }

    for( auto u = this->adj_list[v].begin(); u != this->adj_list[v].end(); u++){
        if((*pi_k)[*u].get_pos()){
            b+=(*pi_k)[*u].get_pos();
            K++;
        }
    }
    b/= K > 0 ? K : 1;

    return b;
}

int BGraph::n_cross(){
    int c = 0;
    int i,j,k,l;
    i=j=k=l=0;
    int aloc;

    int p_i, p_j, p_k, p_l;

    for (int m=0; m<this->n_edges; m++){
        i = this->edgelist[m].u;
        k = this->edgelist[m].v;
        for (int m_=m+1; m_<this->n_edges; m_++){
            j = this->edgelist[m_].u;
            l = this->edgelist[m_].v;


            p_i = this->pi_1[i].get_pos();
            p_j = this->pi_1[j].get_pos();
            p_k = this->pi_2[k].get_pos();
            p_l = this->pi_2[l].get_pos();

            aloc = ( p_i && p_j && p_k && p_l) ? 1 : 0 ;  // Verifying if the all position is allocated ie != 0.

            if (aloc){
                if (
                    (p_i < p_j) &&
                    (p_k > p_l)
                ) {c++;} else if (
                    (p_i > p_j) &&
                    (p_k < p_l)
                ) {c++;}
            }
        }
    }
    return c;
}


void BGraph::move_vertex(int v, int to){

    int k;
    int from;

    VSorted *pos_assing;
    Dict *pi;

    k = this->layer[v];

    switch(k){
        case 1:
            from = this->pi_1[v].get_pos();
            pos_assing = &this->map_v1;
            pi = &this->pi_1;
        break;

        case 2:
            from = this->pi_2[v].get_pos();
            pos_assing = &this->map_v2;
            pi = &this->pi_2;
        break;
    }

    int pos, u, c;

    c = to < from ? -1 : 1;

    for (pos = from; pos != to; pos += c){
        u = (*pos_assing)[pos + c].get_vertex();
        (*pi)[u].setPos(pos);
        (*pos_assing)[pos] = (*pi)[u];
    }
    (*pi)[v].setPos(to);
    (*pos_assing)[to] = (*pi)[v];
}

Dict BGraph::move_vertex(int v, int to, int inplace){
    int k;
    int from;

    VSorted pos_assing;
    Dict pi;

    k = this->layer[v];

    switch(k){
        case 1:
            from = this->pi_1[v].get_pos();
            pos_assing = this->map_v1;
            pi = this->pi_1;
        break;

        case 2:
            from = this->pi_2[v].get_pos();
            pos_assing = this->map_v2;
            pi = this->pi_2;
        break;
    }

    int pos, u, c;

    c = to < from ? -1 : 1;

    for (pos = from; pos != to; pos += c){
        u = (pos_assing)[pos + c].get_vertex();
        (pi)[u].setPos(pos);
        //pi[u] = *(new Vertex(u, pos));
        (pos_assing)[pos] = (pi)[u];
    }
    (pi)[v].setPos(to);
    //pi[v] = *(new Vertex(v, to));
    (pos_assing)[to] = (pi)[v];

    if (inplace) {
        switch(k){
            case 1:
                this->map_v1 = pos_assing;
                this->pi_1 = pi;
            break;

            case 2:
                this->map_v2 = pos_assing;
                this->pi_2 = pi;
            break;
        }
    } else {
        this->maped = k;
    }
    return pi;
}

void BGraph::re_map(int k){

    VSorted new_assing;
    VSorted *pos_assing;
    Dict *pi;


    switch(k){
        case 1:
            pos_assing = &this->map_v1;
            pi = &this->pi_1;
        break;

        case 2:
            pos_assing = &this->map_v2;
            pi = &this->pi_2;
        break;
    }

    for(auto i = (*pi).begin(); i != (*pi).end(); i++){
        new_assing[i->second.get_pos()] = i->second;
    }
    (*pos_assing) = new_assing;
    this->maped = 0;
}

void BGraph::printBGraph(){
    if (this->maped) { re_map(this->maped); }

    string u;
    string v;

    int ut;
    int vt;
    int n_max = (this->n_v1 > this->n_v2 ? this->n_v1 : this->n_v2);
    for (int i=1;i <= n_max ; i++){
        //u = this->map_v1[i].get_vertex() !=-1  ? (char) this->map_v1[i].get_vertex() + '0' : (char) 0 ;
        //v = this->map_v2[i].get_vertex() !=-1  ? (char) this->map_v2[i].get_vertex() + '0' : (char) 0 ;

        ut = (const) this->map_v1[i].get_vertex();
        vt = (const) this->map_v2[i].get_vertex();

        u = (const) this->map_v1[i].get_vertex() !=-1  ? to_string((const) this->map_v1[i].get_vertex()) : "" ;
        v = (const) this->map_v2[i].get_vertex() !=-1  ? to_string((const) this->map_v2[i].get_vertex()) : "" ;
        cout <<  u << " "
             <<  v << endl;
    }
}

void BGraph::printGraph(){
    int n = n_v1 + n_v2;
    for (auto i = this->adj_list.begin() ; i != this->adj_list.end(); i++) {
        auto lst = i->second;

        cout << endl << "Adjacency list of vertex  : " << i->first << "deg: " << lst.size()  << endl;

        for (auto itr = lst.begin(); itr != lst.end(); itr++)
            cout << *itr << " ";
        cout << endl;
    }
}

void BGraph::printMatrix(){
    int m = this->n_edges;
    int n = this->n_v1+this->n_v2;
    for (int i=1; i<=n; i++){
        for (int j=1; j<=n; j++ ){
            cout << this->edges[i][j] << ",";
        }
    cout << endl;
    }
}

BGraph::~BGraph(){

}
