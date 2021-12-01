#include <iostream>
#include <stdlib.h>

#include <bits/stdc++.h>

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>

#include <random>

#define print(txt) cout << txt << endl;
using namespace std;

class Vertex{

    private:
        int v;

    public:
        int pos;

        Vertex(){
            this->v = -1;
            //this->pos = new int;
            //(*this->pos) = 0;
            this->pos = 0;
        }

        Vertex(int i, int p){
            this->v = i;
            //this->pos = new int;
            //(*this->pos) = p;
            this->pos = p;
        }

        //~Vertex(){
        //    delete this->pos;
        //}

        void setVertex(int i, int p){
            this->v = i;
            //this->pos = new int;
            //(*this->pos) = p;
            this->pos = p;
        }

        void setPos(int p){
            //(*this->pos) = p;
            this->pos = p;
        }

        int get_vertex() const {
            return v;
        }

        int get_pos() const {
            //return (*pos);
            return pos;
        }


};

class Solution {
   public:
   int n;
   vector <int> v;
   Solution(vector<int> &w) {
      srand(time(NULL));
      n = w[0];
      for(int i = 1; i < w.size(); i++){
         w[i] += w[i - 1];
         n = w[i];
      }
      v = w;
   }
   int pickIndex() {
      //cout << rand() % 9 << endl;
     // cout << rand() % v.back() << endl;
      //cout << "v back: " << v.back() << endl;
      return upper_bound(v.begin(), v.end(), rand() % v.back()) - v.begin();
   }
};

typedef unordered_map<int, unordered_map<int,int>> AdjMatrix;

typedef unordered_set<int> Adj;
typedef unordered_map<int, vector<int> > Adj2;

typedef unordered_map<int, Vertex> Dict;
typedef map<int, Vertex> VSorted;



typedef struct Edge{
    int u;
    int v;
} edge;

auto mcmp = [](pair<Dict,int> left, pair<Dict,int> right) { return left.second > right.second; }; // Funcao de Comparacao
typedef priority_queue<pair<Dict,int>, vector<pair<Dict,int>>, decltype(mcmp)> MinQueue;

void print_map(const VSorted mapping) { // NB: pass by value so the print uses a copy
    for(auto i = mapping.begin(); i != mapping.end(); i++){
        cout << i->first << " - " << i->second.get_vertex() << " : "<< i->second.get_pos() << '\n'  ;

    }
}

class BGraph{
    private:
        int n_edges;
        Vertex *v1;  //list -> v in V_1  dim:  n_v1x1
        Vertex *v2;  //list -> v in V_2   dim:  n_v2x1


        Edge *edgelist;
        AdjMatrix edges;
        //int **edges; //edges -> adj matriz dim: (n_v1 + n_v2)x(n_v1 + n_v2)


    public:
        int n_v1;
        int n_v2;
        VSorted map_v1;
        VSorted map_v2;
        Dict pi_1;
        Dict pi_2;
        Adj *adj_list;
        Adj2 adj_list2;
        unordered_map<int,int> layer;
        int maped;

        BGraph(int n_v1_, int n_v2_, int n_edges_, Vertex *v1_, Vertex *v2_, Edge *e, int shift = 0){
            this->n_v1 = n_v1_;
            this->n_v2 = n_v2_;
            this->n_edges = n_edges_;

            this->v1 = v1_;
            this->v2 = v2_;
            this->edgelist = e;

            int i,j,m = 0;
            int n = n_v1_ + n_v2_;

            this->adj_list = new Adj[n+shift]; // USANDO HASH SET
            //this->edges = new int* [n];
            //for(i=0; i < n; i++){
            //    this->edges[i] = (int*) calloc(n, sizeof(int));
            //}


            for(m = 0; m < n_edges_; m++){
                i = e[m].u;// - shift;
                j = e[m].v;// - shift;
                this->edges[i][j] = 1;
                this->edges[j][i] = 1;

                //this->adj_list[i].insert(j);
                //this->adj_list[j].insert(i);

                this->adj_list2[i].push_back(j);
                this->adj_list2[j].push_back(i);

               // cout << *(this->adj_list[i].begin()) << endl;


                this->edgelist[m].u = i;// + shift;
                this->edgelist[m].v = j;// + shift;
            }


           // if (this->test.find(10) == this->test.end()){
             //  cout<< "nao existe" <<endl;
            //}


           // for (auto i = this->test.begin(); i != this->test.end(); i++)
             //   cout << i->first << " : " << i->second.size() << endl;



            for (i=1;i<=n_v1_;i++){
                this->map_v1[v1_[i-1].get_pos()] = v1_[i-1];
                cout << i << "(" << v1_[i-1].get_vertex() << ")" << " ~ ";
                this->pi_1[v1_[i-1].get_vertex()] = v1_[i-1];
                this->layer[v1_[i-1].get_vertex()] = 1;
            }
            print("")
            for (i=1;i<=n_v2_;i++){
                this->map_v2[v2_[i-1].get_pos()] = v2_[i-1];
                cout << i << "(" << v2_[i-1].get_vertex() << ")" << " ~ ";
                this->pi_2[v2_[i-1].get_vertex()] = v2_[i-1];
                this->layer[v2_[i-1].get_vertex()] = 2;
            }
            this->maped = 0;
        }

        unordered_map<int, int> degrees(){
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

        unordered_map<int, int> degrees(unordered_set<int> V, unordered_set<int> V_sub){
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

        int greedy_selection(float alpha, unordered_set<int> V, unordered_set<int> V_sub, int seed = -1){
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

        float bc(int v, int k){ //barycenter
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

            for( auto u = this->adj_list2[v].begin(); u != this->adj_list2[v].end(); u++){
                if((*pi_k)[*u].get_pos()){
                    b+=(*pi_k)[*u].get_pos();
                    K++;
                }
            }
            b/= K > 0 ? K : 1;

            return b;
        }

        //~BGraph();

        int n_cross(){
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


        void move_vertex(int v, int to){

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

        Dict move_vertex(int v, int to, int inplace){
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

        void re_map(int k){

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

        void printBGraph(){
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

        void printGraph(){
            int n = n_v1 + n_v2;
            for (int i = 0; i <= n; i++) {
                Adj lst = this->adj_list[i];
                cout << endl << "Adjacency list of vertex  : " << i << " deg: " << lst.size()  << endl;

                for (auto itr = lst.begin(); itr != lst.end(); itr++)
                    cout << *itr << " ";
                cout << endl;
            }
        }

        void printGraph2(){
            int n = n_v1 + n_v2;
            for (auto i = this->adj_list2.begin() ; i != this->adj_list2.end(); i++) {
                auto lst = i->second;

                cout << endl << "Adjacency list of vertex  : " << i->first << "deg: " << lst.size()  << endl;

                for (auto itr = lst.begin(); itr != lst.end(); itr++)
                    cout << *itr << " ";
                cout << endl;
            }
        }

        void printMatrix_2(){
            int m = this->n_edges;
            int n = this->n_v1+this->n_v2;
            for (int i=1; i<=n; i++){
                for (int j=1; j<=n; j++ ){
                    cout << this->edges[i][j] << ",";
                }
            cout << endl;
            }
        }

        void printMatrix(){
            int n = this->n_v1+this->n_v2;
            for(int i=0; i < n; i++){
                for(int j=0; j < n; j++){
                    cout << this->edges[i+1][j+1] << ",";
                }
            }
            cout << endl;
        }

};

void construction_phase(BGraph& G, unordered_set<int> U_1, unordered_set<int> U_2, float alpha, int seed=-1){

    int v, k;

    float bc_v;
    int bc_v_plus;
    int bc_v_minus;
    int pos;

    unsigned int n_cross_plus=-1;
    unsigned int n_cross_minus=-1;

    auto isInteger = [](float N){ return (N - (int) N) > 0 ? 0 : 1; };

    VSorted *pos_assing;
    Dict *pi;

    unordered_set<int> V;
    unordered_set<int> U;

    for (auto i = U_2.begin(); i != U_2.end(); i++){

        U.insert(*i);

        pos = G.pi_2[*i].get_pos();
        G.pi_2[*i].setPos(0);

        G.map_v2[pos] = *(new Vertex());
    }


    for (auto i = U_1.begin(); i != U_1.end(); i++){

        U.insert(*i);

        pos = G.pi_1[*i].get_pos();
        G.pi_1[*i].setPos(0);

        G.map_v1[pos] = *(new Vertex());
    }

    v = G.greedy_selection(alpha, U, U);
    k = G.layer[v];

    switch(k){
        case 1:
            U_1.erase(v);
            G.pi_1[v].setPos(1);
            G.map_v1[1] = G.pi_1[v];
        break;

        case 2:
            U_2.erase(v);
            G.pi_2[v].setPos(1);
            G.map_v2[1] = G.pi_2[v];
        break;
    }

    U.erase(v);
    V.insert(v);
    while (U.size()){

        v = G.greedy_selection(alpha, U, V, seed);
        k = G.layer[v];

        switch(k){
            case 1:
                pos_assing = &G.map_v1;
                pi = &G.pi_1;
            break;

            case 2:
                pos_assing = &G.map_v2;
                pi = &G.pi_2;
            break;
        }

        bc_v = G.bc(v, k);
        bc_v = bc_v > 0 ? bc_v : 1; // Pos check for start at 1st

        if (isInteger(bc_v)){
            bc_v_minus = bc_v;
            bc_v_plus = bc_v;
        } else {
            bc_v_minus = bc_v;
            bc_v_plus = bc_v+1;
        }

        while((n_cross_plus==-1) && (n_cross_minus==-1)){

            if (!((*pos_assing)[bc_v_minus].get_pos()) && bc_v_minus ){
                (*pi)[v].setPos(bc_v_minus);
                n_cross_minus = G.n_cross();
                (*pi)[v].setPos(0);
            }

            if (!((*pos_assing)[bc_v_plus].get_pos()) && (bc_v_plus  <= pos_assing->size()) ){
                (*pi)[v].setPos(bc_v_plus);
                n_cross_plus = G.n_cross();
                (*pi)[v].setPos(0);
            }

            if( !( (n_cross_minus == n_cross_plus) && (n_cross_plus==-1)) ){
                if(n_cross_minus == -1){
                    (*pi)[v].setPos(bc_v_plus);
                    (*pos_assing)[bc_v_plus] = (*pi)[v];
                    pos=bc_v_plus;
                    break;
                } else if ( n_cross_plus <= n_cross_minus  ){
                    (*pi)[v].setPos(bc_v_plus);
                    (*pos_assing)[bc_v_plus] = (*pi)[v];
                    pos=bc_v_plus;
                    break;
                } else {
                    (*pi)[v].setPos(bc_v_minus);
                    (*pos_assing)[bc_v_minus] = (*pi)[v];
                    pos=bc_v_minus;
                    break;
                }
            }
            bc_v_minus = (bc_v_minus-1) > 0 ? --bc_v_minus : 1;
            bc_v_plus  = (bc_v_plus+1) <= pos_assing->size() ? ++bc_v_plus : 1;
         }
        n_cross_plus = -1;
        n_cross_minus = -1;
        pos_assing=nullptr;
        U.erase(v);
        V.insert(v);
    }
}

vector<int> random_choice(BGraph& G, unordered_map<int,int> Degrees, int seed = -1){

    vector<int> vector_list;
    vector<double> pr_list;
    vector<int> sample;

    unordered_set<int> restrict_set;



    for(auto i = Degrees.begin(); i != Degrees.end() ; i++ ){

        vector_list.push_back(i->first);
        pr_list.push_back(i->second);

    }

    //random_device rd;
    mt19937 gen(seed);

    discrete_distribution<int> dist(pr_list.begin(),pr_list.end());

    vector<double> prob;

    int v;
    int s;
    int i=0;
    while( sample.size() < vector_list.size() ){
        //dist = discrete_distribution<int>(pr_list.begin(),pr_list.end());

        /*
        cout << '\n';
        prob = dist.probabilities();
        for(auto n : prob)
            cout << n << ' ';
        cout << '\n';*/

        s = dist(gen);
        pr_list[s] = 0;

        v = vector_list[s];
        auto p = restrict_set.insert(v);

        if( p.second ){
            //cout //<< "(" << s << ")"
            //<< v << " / ";
            sample.push_back(v);
        }
        i++;
    }
    //cout << endl;
    //cout << i ;

    return sample;
    //while(Degrees.size()){
    //}

}

template<typename T>
void print_queue(T q) { // NB: pass by value so the print uses a copy
    while(!q.empty()) {
        cout << q.top().second << ' ';
        q.pop();
    }
    cout << '\n';
}

void improvement_phase(BGraph& G, int seed =-1, unordered_set<int> U_1 = unordered_set<int>() , unordered_set<int> U_2 =  unordered_set<int>()){

    auto isInteger = [](float N){ return (N - (int) N) > 0 ? 0 : 1; };

    unordered_map<int,int> deg;

    unordered_set<int> U;

    vector<int> sample;

    if (U_1.size() && U_2.size()){
        for (auto i = U_2.begin(); i != U_2.end(); i++){
            U.insert(*i);
        }
        for (auto i = U_1.begin(); i != U_1.end(); i++){
            U.insert(*i);
        }
        deg = G.degrees(U, U);
    } else {
        deg = G.degrees();
    }

    int k,v;

    Dict pi_aux;
    Dict pi_plus;
    Dict pi_minus;
    Dict pi_bc;

    unsigned int n_cross_aux;
    unsigned int n_cross_bc;
    unsigned int n_cross_plus;
    unsigned int n_cross_minus;


    float bc;
    int v_bc;
    int v_plus;
    int v_minus;

    // MIN HEAP
    MinQueue *min_q;

    sample = random_choice(G, deg, seed);

    for (auto i = sample.begin(); i != sample.end(); i++){

        n_cross_aux=-1;
        n_cross_bc=-1;
        n_cross_plus=-1;
        n_cross_minus=-1;

        v_bc = 0;
        v_plus = 0;
        v_minus = 0;

        v = *i;
        k = G.layer[v];
        bc = G.bc(v, k);

        min_q = new MinQueue(mcmp);
        switch(k){
            case 1:

                if (  isInteger(bc)   ){
                    v_bc = bc;
                    v_minus = v_bc - 1;
                    v_plus =  v_bc + 1 > G.n_v2 ? 0 : v_bc + 1 ;
                } else if ( ((int) bc == 1) ){
                    v_bc = bc + 1;
                    v_minus = v_bc - 1;
                    v_plus =  v_bc + 1 > G.n_v2 ? 0 : v_bc + 1 ;
                } else if ( (int) bc == G.n_v2 ) {
                    v_bc = bc;
                    v_minus = v_bc - 1 > 0 ;
                    v_plus =  v_bc + 1 > G.n_v2 ? 0 : v_bc + 1 ;
                } // Think Better this else

                //int_bc = int_bc > G.n_v2 - 1
                pi_aux = G.pi_1;
                n_cross_aux = G.n_cross();
                min_q->push({pi_aux, n_cross_aux});

                if (v_bc){
                    G.pi_1 = pi_aux;
                    pi_bc = G.move_vertex(v, v_bc, 0);
                    G.pi_1 = pi_bc;
                    n_cross_bc = G.n_cross();
                    min_q->push({pi_bc, n_cross_bc});
                }

                if (v_plus){
                    G.pi_1 = pi_aux;
                    pi_plus = G.move_vertex(v, v_plus, 0);
                    G.pi_1 = pi_plus;
                    n_cross_plus = G.n_cross();
                    min_q->push({pi_plus, n_cross_plus});
                }

                if (v_minus){
                    G.pi_1 = pi_aux;
                    pi_minus = G.move_vertex(v, v_minus, 0);
                    G.pi_1 = pi_minus;
                    n_cross_minus = G.n_cross();
                    min_q->push({pi_minus, n_cross_minus});
                }

                //print_queue(*min_q);

                G.pi_1 = min_q->top().first;
                G.re_map(1);

            break;

            case 2:

                if (  isInteger(bc)   ){
                    v_bc = bc;
                    v_minus = v_bc - 1;
                    v_plus =  v_bc + 1 > G.n_v1 ? 0 : v_bc + 1 ;
                } else if ( ((int) bc == 1) ){
                    v_bc = bc + 1;
                    v_minus = v_bc - 1;
                    v_plus =  v_bc + 1 > G.n_v1 ? 0 : v_bc + 1 ;
                } else if ( (int) bc == G.n_v1 ) {
                    v_bc = bc;
                    v_minus = v_bc - 1 > 0 ;
                    v_plus =  v_bc + 1 > G.n_v1 ? 0 : v_bc + 1 ;
                } // Think Better this else

                pi_aux = G.pi_2;
                n_cross_aux = G.n_cross();
                min_q->push({pi_aux, n_cross_aux});

                if (v_bc){
                    G.pi_2 = pi_aux;
                    pi_bc = G.move_vertex(v, v_bc, 0);
                    G.pi_2 = pi_bc;
                    n_cross_bc = G.n_cross();
                    min_q->push({pi_bc, n_cross_bc});
                }

                if (v_plus){
                    G.pi_2 = pi_aux;
                    pi_plus = G.move_vertex(v, v_plus, 0);
                    G.pi_2 = pi_plus;
                    n_cross_plus = G.n_cross();
                    min_q->push({pi_plus, n_cross_plus});
                }

                if (v_minus){
                    G.pi_2 = pi_aux;
                    pi_minus = G.move_vertex(v, v_minus, 0);
                    G.pi_2 = pi_minus;
                    n_cross_minus = G.n_cross();
                    min_q->push({pi_minus, n_cross_minus});
                }

                //print_queue(*min_q);

                G.pi_2 = min_q->top().first;
                G.re_map(2);


            break;
        }
        if ( !min_q->top().second ) { return ;}
    }
}

int main(){

//    Edge e[] = {{1,6}, {2,6}, {3,7}, {3,8}, {4,7}, {5,6}, {5,9}, {2,9},{4,6}};
//    Edge e[] = {{0, 49}, {0, 41}, {0, 42}, {1, 57}, {1, 44}, {1, 60}, {2, 53}, {2, 42}, {2, 60}, {3, 36}, {3, 42}, {3, 59}, {3, 44}, {4, 35}, {4, 47}, {5, 35}, {5, 46}, {6, 51}, {6, 40}, {7, 37}, {8, 34}, {9, 54}, {9, 43}, {9, 31}, {10, 50}, {11, 54}, {11, 55}, {11, 56}, {12, 48}, {12, 39}, {13, 57}, {13, 44}, {14, 43}, {15, 40}, {16, 48}, {17, 37}, {17, 58}, {17, 43}, {18, 52}, {18, 58}, {19, 59}, {19, 45}, {20, 32}, {21, 38}, {22, 41}, {22, 60}, {23, 33}, {23, 55}, {23, 57}, {23, 42}, {24, 55}, {24, 42}, {25, 35}, {25, 58}, {26, 35}, {26, 59}, {27, 51}, {27, 55}, {28, 40}, {29, 57}, {29, 61}, {30, 54}, {30, 47}};
    Edge e[] = {{0, 185}, {0, 186}, {1, 187}, {1, 188}, {2, 189}, {3, 190}, {4, 189}, {5, 191}, {6, 192}, {7, 193}, {8, 194}, {9, 195}, {10, 196}, {11, 197}, {11, 198}, {12, 192}, {13, 199}, {13, 200}, {14, 185}, {15, 187}, {16, 201}, {17, 198}, {18, 202}, {19, 199}, {20, 203}, {21, 185}, {21, 203}, {22, 191}, {22, 198}, {23, 202}, {23, 204}, {24, 205}, {25, 206}, {26, 186}, {27, 207}, {27, 208}, {28, 196}, {29, 193}, {30, 209}, {31, 210}, {32, 211}, {33, 212}, {34, 213}, {35, 214}, {36, 202}, {37, 208}, {38, 192}, {38, 193}, {39, 215}, {40, 211}, {41, 210}, {42, 216}, {43, 186}, {44, 205}, {45, 217}, {46, 206}, {47, 218}, {47, 192}, {48, 219}, {49, 220}, {50, 221}, {51, 208}, {52, 222}, {52, 208}, {53, 223}, {54, 205}, {55, 219}, {55, 195}, {56, 201}, {56, 200}, {57, 211}, {58, 186}, {59, 218}, {59, 224}, {60, 186}, {60, 225}, {61, 194}, {62, 220}, {63, 203}, {64, 194}, {64, 212}, {65, 204}, {66, 212}, {66, 214}, {67, 197}, {68, 200}, {69, 204}, {70, 200}, {71, 209}, {72, 193}, {72, 211}, {73, 201}, {73, 194}, {74, 226}, {75, 200}, {76, 222}, {77, 190}, {77, 210}, {78, 223}, {79, 224}, {80, 216}, {81, 218}, {82, 225}, {83, 210}, {84, 188}, {85, 221}, {86, 212}, {87, 214}, {88, 226}, {88, 195}, {89, 207}, {90, 195}, {91, 220}, {92, 219}, {93, 187}, {93, 207}, {94, 226}, {95, 202}, {95, 225}, {96, 224}, {97, 215}, {98, 197}, {98, 199}, {99, 223}, {100, 189}, {100, 216}, {101, 222}, {102, 206}, {102, 227}, {103, 198}, {104, 217}, {104, 221}, {105, 209}, {106, 228}, {107, 190}, {108, 225}, {109, 191}, {110, 213}, {111, 227}, {112, 227}, {112, 215}, {113, 187}, {114, 223}, {114, 214}, {115, 185}, {116, 229}, {117, 207}, {118, 220}, {118, 205}, {119, 185}, {120, 225}, {121, 218}, {122, 219}, {122, 228}, {123, 192}, {124, 188}, {125, 209}, {125, 228}, {126, 194}, {127, 202}, {128, 197}, {129, 204}, {129, 230}, {130, 196}, {131, 222}, {132, 221}, {133, 229}, {134, 189}, {134, 210}, {135, 203}, {136, 201}, {137, 224}, {138, 226}, {138, 229}, {139, 224}, {139, 216}, {140, 227}, {141, 199}, {142, 230}, {143, 213}, {144, 214}, {145, 215}, {146, 199}, {147, 230}, {148, 206}, {149, 213}, {150, 217}, {151, 203}, {151, 211}, {152, 188}, {153, 221}, {153, 206}, {154, 193}, {155, 228}, {156, 226}, {157, 191}, {157, 196}, {158, 201}, {159, 229}, {160, 216}, {161, 198}, {162, 213}, {162, 230}, {163, 219}, {164, 188}, {164, 229}, {165, 207}
};

   /*
    for (auto i = order.begin();
         i != order.end(); i++)
    {
        cout << i->first
                  << " : "
                  << i->second << '\n';
    }*/

    /*
    vector<int> c = {1,3,1000};
    Solution ob(c);
    while (1){
    cout << ob.pickIndex() << " ";
    ci*n.ignore();
    }*/

    Vertex *v1;//v1[5];
    Vertex *v2;//v2[4];
    v1 = new Vertex[166];
    v2 = new Vertex[46];

    /*
    v1[0].setVertex(1, 1);
    v1[1].setVertex(2, 3);
    v1[2].setVertex(3, 2);
    v1[3].setVertex(4, 4);
    v1[4].setVertex(5, 5);


    v2[0].setVertex(6, 1);
    v2[1].setVertex(7, 2);
    v2[2].setVertex(8, 3);
    v2[3].setVertex(9, 4);


    BGraph New(5, 4, 7,
              v1, v2, e,
              1);

    v1[0].setVertex(0,1);	v2[0].setVertex(31,1);
    v1[1].setVertex(1,2);	v2[1].setVertex(32,2);
    v1[2].setVertex(2,3);	v2[2].setVertex(33,3);
    v1[3].setVertex(3,4);	v2[3].setVertex(34,4);
    v1[4].setVertex(4,5);	v2[4].setVertex(35,5);
    v1[5].setVertex(5,6);	v2[5].setVertex(36,6);
    v1[6].setVertex(6,7);	v2[6].setVertex(37,7);
    v1[7].setVertex(7,8);	v2[7].setVertex(38,8);
    v1[8].setVertex(8,9);	v2[8].setVertex(39,9);
    v1[9].setVertex(9,10);	v2[9].setVertex(40,10);
    v1[10].setVertex(10,11);	v2[10].setVertex(41,11);
    v1[11].setVertex(11,12);	v2[11].setVertex(42,12);
    v1[12].setVertex(12,13);	v2[12].setVertex(43,13);
    v1[13].setVertex(13,14);	v2[13].setVertex(44,14);
    v1[14].setVertex(14,15);	v2[14].setVertex(45,15);
    v1[15].setVertex(15,16);	v2[15].setVertex(46,16);
    v1[16].setVertex(16,17);	v2[16].setVertex(47,17);
    v1[17].setVertex(17,18);	v2[17].setVertex(48,18);
    v1[18].setVertex(18,19);	v2[18].setVertex(49,19);
    v1[19].setVertex(19,20);	v2[19].setVertex(50,20);
    v1[20].setVertex(20,21);	v2[20].setVertex(51,21);
    v1[21].setVertex(21,22);	v2[21].setVertex(52,22);
    v1[22].setVertex(22,23);	v2[22].setVertex(53,23);
    v1[23].setVertex(23,24);	v2[23].setVertex(54,24);
    v1[24].setVertex(24,25);	v2[24].setVertex(55,25);
    v1[25].setVertex(25,26);	v2[25].setVertex(56,26);
    v1[26].setVertex(26,27);	v2[26].setVertex(57,27);
    v1[27].setVertex(27,28);	v2[27].setVertex(58,28);
    v1[28].setVertex(28,29);	v2[28].setVertex(59,29);
    v1[29].setVertex(29,30);	v2[29].setVertex(60,30);
    v1[30].setVertex(30,31);	v2[30].setVertex(61,31);
*/
    v1[0].setVertex(0,1);	v2[0].setVertex(185,1);
    v1[1].setVertex(1,2);	v2[1].setVertex(186,2);
    v1[2].setVertex(2,3);	v2[2].setVertex(187,3);
    v1[3].setVertex(3,4);	v2[3].setVertex(188,4);
    v1[4].setVertex(4,5);	v2[4].setVertex(189,5);
    v1[5].setVertex(5,6);	v2[5].setVertex(190,6);
    v1[6].setVertex(6,7);	v2[6].setVertex(191,7);
    v1[7].setVertex(7,8);	v2[7].setVertex(192,8);
    v1[8].setVertex(8,9);	v2[8].setVertex(193,9);
    v1[9].setVertex(9,10);	v2[9].setVertex(194,10);
    v1[10].setVertex(10,11);	v2[10].setVertex(195,11);
    v1[11].setVertex(11,12);	v2[11].setVertex(196,12);
    v1[12].setVertex(12,13);	v2[12].setVertex(197,13);
    v1[13].setVertex(13,14);	v2[13].setVertex(198,14);
    v1[14].setVertex(14,15);	v2[14].setVertex(199,15);
    v1[15].setVertex(15,16);	v2[15].setVertex(200,16);
    v1[16].setVertex(16,17);	v2[16].setVertex(201,17);
    v1[17].setVertex(17,18);	v2[17].setVertex(202,18);
    v1[18].setVertex(18,19);	v2[18].setVertex(203,19);
    v1[19].setVertex(19,20);	v2[19].setVertex(204,20);
    v1[20].setVertex(20,21);	v2[20].setVertex(205,21);
    v1[21].setVertex(21,22);	v2[21].setVertex(206,22);
    v1[22].setVertex(22,23);	v2[22].setVertex(207,23);
    v1[23].setVertex(23,24);	v2[23].setVertex(208,24);
    v1[24].setVertex(24,25);	v2[24].setVertex(209,25);
    v1[25].setVertex(25,26);	v2[25].setVertex(210,26);
    v1[26].setVertex(26,27);	v2[26].setVertex(211,27);
    v1[27].setVertex(27,28);	v2[27].setVertex(212,28);
    v1[28].setVertex(28,29);	v2[28].setVertex(213,29);
    v1[29].setVertex(29,30);	v2[29].setVertex(214,30);
    v1[30].setVertex(30,31);	v2[30].setVertex(215,31);
    v1[31].setVertex(31,32);	v2[31].setVertex(216,32);
    v1[32].setVertex(32,33);	v2[32].setVertex(217,33);
    v1[33].setVertex(33,34);	v2[33].setVertex(218,34);
    v1[34].setVertex(34,35);	v2[34].setVertex(219,35);
    v1[35].setVertex(35,36);	v2[35].setVertex(220,36);
    v1[36].setVertex(36,37);	v2[36].setVertex(221,37);
    v1[37].setVertex(37,38);	v2[37].setVertex(222,38);
    v1[38].setVertex(38,39);	v2[38].setVertex(223,39);
    v1[39].setVertex(39,40);	v2[39].setVertex(224,40);
    v1[40].setVertex(40,41);	v2[40].setVertex(225,41);
    v1[41].setVertex(41,42);	v2[41].setVertex(226,42);
    v1[42].setVertex(42,43);	v2[42].setVertex(227,43);
    v1[43].setVertex(43,44);	v2[43].setVertex(228,44);
    v1[44].setVertex(44,45);	v2[44].setVertex(229,45);
    v1[45].setVertex(45,46);	v2[45].setVertex(230,46);
    v1[46].setVertex(46,47);
    v1[47].setVertex(47,48);
    v1[48].setVertex(48,49);
    v1[49].setVertex(49,50);
    v1[50].setVertex(50,51);
    v1[51].setVertex(51,52);
    v1[52].setVertex(52,53);
    v1[53].setVertex(53,54);
    v1[54].setVertex(54,55);
    v1[55].setVertex(55,56);
    v1[56].setVertex(56,57);
    v1[57].setVertex(57,58);
    v1[58].setVertex(58,59);
    v1[59].setVertex(59,60);
    v1[60].setVertex(60,61);
    v1[61].setVertex(61,62);
    v1[62].setVertex(62,63);
    v1[63].setVertex(63,64);
    v1[64].setVertex(64,65);
    v1[65].setVertex(65,66);
    v1[66].setVertex(66,67);
    v1[67].setVertex(67,68);
    v1[68].setVertex(68,69);
    v1[69].setVertex(69,70);
    v1[70].setVertex(70,71);
    v1[71].setVertex(71,72);
    v1[72].setVertex(72,73);
    v1[73].setVertex(73,74);
    v1[74].setVertex(74,75);
    v1[75].setVertex(75,76);
    v1[76].setVertex(76,77);
    v1[77].setVertex(77,78);
    v1[78].setVertex(78,79);
    v1[79].setVertex(79,80);
    v1[80].setVertex(80,81);
    v1[81].setVertex(81,82);
    v1[82].setVertex(82,83);
    v1[83].setVertex(83,84);
    v1[84].setVertex(84,85);
    v1[85].setVertex(85,86);
    v1[86].setVertex(86,87);
    v1[87].setVertex(87,88);
    v1[88].setVertex(88,89);
    v1[89].setVertex(89,90);
    v1[90].setVertex(90,91);
    v1[91].setVertex(91,92);
    v1[92].setVertex(92,93);
    v1[93].setVertex(93,94);
    v1[94].setVertex(94,95);
    v1[95].setVertex(95,96);
    v1[96].setVertex(96,97);
    v1[97].setVertex(97,98);
    v1[98].setVertex(98,99);
    v1[99].setVertex(99,100);
    v1[100].setVertex(100,101);
    v1[101].setVertex(101,102);
    v1[102].setVertex(102,103);
    v1[103].setVertex(103,104);
    v1[104].setVertex(104,105);
    v1[105].setVertex(105,106);
    v1[106].setVertex(106,107);
    v1[107].setVertex(107,108);
    v1[108].setVertex(108,109);
    v1[109].setVertex(109,110);
    v1[110].setVertex(110,111);
    v1[111].setVertex(111,112);
    v1[112].setVertex(112,113);
    v1[113].setVertex(113,114);
    v1[114].setVertex(114,115);
    v1[115].setVertex(115,116);
    v1[116].setVertex(116,117);
    v1[117].setVertex(117,118);
    v1[118].setVertex(118,119);
    v1[119].setVertex(119,120);
    v1[120].setVertex(120,121);
    v1[121].setVertex(121,122);
    v1[122].setVertex(122,123);
    v1[123].setVertex(123,124);
    v1[124].setVertex(124,125);
    v1[125].setVertex(125,126);
    v1[126].setVertex(126,127);
    v1[127].setVertex(127,128);
    v1[128].setVertex(128,129);
    v1[129].setVertex(129,130);
    v1[130].setVertex(130,131);
    v1[131].setVertex(131,132);
    v1[132].setVertex(132,133);
    v1[133].setVertex(133,134);
    v1[134].setVertex(134,135);
    v1[135].setVertex(135,136);
    v1[136].setVertex(136,137);
    v1[137].setVertex(137,138);
    v1[138].setVertex(138,139);
    v1[139].setVertex(139,140);
    v1[140].setVertex(140,141);
    v1[141].setVertex(141,142);
    v1[142].setVertex(142,143);
    v1[143].setVertex(143,144);
    v1[144].setVertex(144,145);
    v1[145].setVertex(145,146);
    v1[146].setVertex(146,147);
    v1[147].setVertex(147,148);
    v1[148].setVertex(148,149);
    v1[149].setVertex(149,150);
    v1[150].setVertex(150,151);
    v1[151].setVertex(151,152);
    v1[152].setVertex(152,153);
    v1[153].setVertex(153,154);
    v1[154].setVertex(154,155);
    v1[155].setVertex(155,156);
    v1[156].setVertex(156,157);
    v1[157].setVertex(157,158);
    v1[158].setVertex(158,159);
    v1[159].setVertex(159,160);
    v1[160].setVertex(160,161);
    v1[161].setVertex(161,162);
    v1[162].setVertex(162,163);
    v1[163].setVertex(163,164);
    v1[164].setVertex(164,165);
    v1[165].setVertex(165,166);


    BGraph New(166, 46, 207,
              v1, v2, e,
              0);

    cout << "Crossing : " << New.n_cross();
    /*print("matrix")
    New.printMatrix_2();

    print("TESTE <SADASDADADASDASD")
    print_map(New.map_v1);
    //print("test_deg");

    //New.printGraph();

    print("test_deg_2");

    New.printGraph2();
*/
    print("deg_test");
    unordered_map<int,int> aux;

/*
    unordered_set<int> U_1({1,2,3,4,5});
    unordered_set<int> U_2({6,7,8,9});

    unordered_set<int> V({3,4,7,8});
    unordered_set<int> V_sub({1,2,5,6,9});

    unordered_set<int> U_1({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30});
    unordered_set<int> U_2({31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61});
*/
    unordered_set<int> U_1, U_2;

    for (auto u = New.pi_1.begin(); u != New.pi_1.end(); u++ ){
        U_1.insert(u->first);
    }
    for (auto u_ = New.pi_2.begin(); u_ != New.pi_2.end(); u_++ ){
        U_2.insert(u_->first);
    }

    unordered_set<int> V({3,4,7,8});
    unordered_set<int> V_sub({1,2,5,6,9});

    print("union")
    srand(1);

    //print(New.greedy_selection(1, V, V_sub));

    print("################");
    print("MOVE");

    //auto mcmp = [](pair<Dict,int> left, pair<Dict,int> right) { return left.second > right.second; }; // Funcao de Comparacao
    //priority_queue<pair<Dict,int>, vector<pair<Dict,int>>, decltype(mcmp)> min_q(mcmp);

    MinQueue *test_pointer_queue;
    MinQueue min_q(mcmp);

    test_pointer_queue = &min_q;

    New.printBGraph();
    print("Crossing " << New.n_cross());


    print_queue(*test_pointer_queue);
    print(test_pointer_queue->empty());

    test_pointer_queue = new MinQueue(mcmp);
    print_queue(*test_pointer_queue);

    print(test_pointer_queue->empty());
    print("################");
    print("#####Construct######");
    //construction_phase(New, U_1, U_2, 1);
    New.printBGraph();
    print("Crossing:" << New.n_cross());


    clock_t t_start, t_end;
    aux = New.degrees();


    random_choice(New, aux, 42);

    t_start = clock();
/*    construction_phase(New, U_1, U_2, 0.8, 42);
    improvement_phase(New, 42);

    construction_phase(New, U_1, U_2, 0.8, 50);
    improvement_phase(New, 50);

    construction_phase(New, U_1, U_2, 0.8, 47);
    improvement_phase(New, 47);

    construction_phase(New, U_1, U_2, 0.8, 88);
    improvement_phase(New, 88);

    construction_phase(New, U_1, U_2, 0.8, 99);
    improvement_phase(New, 99);
*/

    for(int i=0; i<30; i++){
        construction_phase(New, U_1, U_2, 0.8, i);
        improvement_phase(New, i);
        print("Crossing:" << New.n_cross());
    }

    t_end = clock();






    New.printBGraph();
    print("Crossing:" << New.n_cross());



    print("time")
    double time_taken = double(t_end - t_start) / double(CLOCKS_PER_SEC);;
    cout << "Time taken by program is : " << fixed
         << time_taken << setprecision(16);
    cout << " sec " << endl;

    auto cmp = [](pair<int,int> left, pair<int,int> right) { return left.second < right.second; };

    priority_queue<pair<int,int>, vector<pair<int,int>>, decltype(cmp)> q3(cmp);

    for(auto n = aux.begin(); n!=aux.end();n++){
        q3.push(*n);
        cout << n->first << ":";
        std::cout << q3.top().first<<endl;}
        //q3.push(make_pair(n->first,n->second));

    while(!q3.empty()) {
        std::cout << q3.top().first << ' ';
        q3.pop();
    }


    delete v1;
    delete v2;
    return 0;
}
