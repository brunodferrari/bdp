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

                this->adj_list[i].insert(j);
                this->adj_list[j].insert(i);

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
            int i,j;
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

            //cout << " Deg size " << deg_dict.size() << endl ;

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
                    //cout << "   TOP "<<max_q.top().first << " ";
                max_idx++;
                max_q.pop();
            }
            //cout << endl;
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
            //print(" asj barycenter");
            //cout << adj_list2[v].size();
            //print("");

            //cout << v << " : " << (*this->pi_1[v].pos) << endl;
            //cout << v << " : " << (*this->pi_2[v].pos) << endl;


            for( auto u = this->adj_list2[v].begin(); u != this->adj_list2[v].end(); u++){
                //print(" barycenter");
                //cout << " barycenter" << *u << " : " << (*pi_k)[*u].get_pos() << endl;
                if((*pi_k)[*u].get_pos()){
                    b+=(*pi_k)[*u].get_pos();
                    K++;
                    //print(this->pi_2[6].get_pos());
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

                    aloc = ( p_i && p_j && p_k && p_l) ? 1 : 0 ;  // Verifying if the all position is alocated ie != 0.

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

            //int i;
            //cin >> i;
            //if (i==1){
              //  printBGraph();
            //}
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

                ut = this->map_v1[i].get_vertex();
                vt = this->map_v2[i].get_vertex();

                u = this->map_v1[i].get_vertex() !=-1  ? to_string(this->map_v1[i].get_vertex()) : "" ;
                v = this->map_v2[i].get_vertex() !=-1  ? to_string(this->map_v2[i].get_vertex()) : "" ;
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

        //cout << *i<<"("<<(G.pi_2[*i].get_pos()) << ") : " << endl ;
        U.insert(*i);

        pos = G.pi_2[*i].get_pos();
        G.pi_2[*i].setPos(0);

        G.map_v2[pos] = *(new Vertex());
    }


    for (auto i = U_1.begin(); i != U_1.end(); i++){

        //cout << *i<<"("<<(G.pi_1[*i].get_pos()) << ") : " << endl ; //print exclude

        U.insert(*i);

        pos = G.pi_1[*i].get_pos();
        G.pi_1[*i].setPos(0);

        G.map_v1[pos] = *(new Vertex());
    }

    /*
    for (auto i = U.begin(); i != U.end(); i++){
        cout << *i << " : ";
        G.pi_1.find(*i) == G.pi_1.end() ? cout << (G.pi_2[*i].get_pos()) << " /" : cout << (G.pi_1[*i].get_pos()) << " /";
        cout << endl;
    }
*/
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
    //print("check pos");
    //cout << G.map_v1[3].get_pos() << endl;
    //cout << G.n_cross() << endl;

    //cout << G.pi_2[v].get_pos() << endl;

    U.erase(v);
    V.insert(v);
    while (U.size()){


       // cout<<"SIZE V:" << V.size() << endl;
        //cout<<"SIZE U:" << U.size() << endl;
        v = G.greedy_selection(alpha, U, V);
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
        //cout << "Layer " << k << " Size " << pos_assing->size() << endl;

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
/*
            print("pos_assing _minus");
            cout << bc_v_minus << " : " << ((*pos_assing)[bc_v_minus].get_pos())
                                        << "("<< ((*pos_assing)[bc_v_minus].get_vertex()) << ")"
                                        << endl;

            print("pos_assing _plus");
            cout << bc_v_plus << " : "  << ((*pos_assing)[bc_v_plus].get_pos())
                                        << "(" <<((*pos_assing)[bc_v_plus].get_vertex()) << ")"
                                        << endl;*/

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

        //cout << i->first << "(" << i->second << ") * ";

        vector_list.push_back(i->first);
        pr_list.push_back(i->second);

    }

    cout << endl;

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
            cout //<< "(" << s << ")"
            << v << " / ";
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

    if (U_1.size() * U_2.size()){
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
    //print("")
    //print("############SAMPLE CHECK################");
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
    Edge e[] = {{0, 49}, {0, 41}, {0, 42}, {1, 57}, {1, 44}, {1, 60}, {2, 53}, {2, 42}, {2, 60}, {3, 36}, {3, 42}, {3, 59}, {3, 44}, {4, 35}, {4, 47}, {5, 35}, {5, 46}, {6, 51}, {6, 40}, {7, 37}, {8, 34}, {9, 54}, {9, 43}, {9, 31}, {10, 50}, {11, 54}, {11, 55}, {11, 56}, {12, 48}, {12, 39}, {13, 57}, {13, 44}, {14, 43}, {15, 40}, {16, 48}, {17, 37}, {17, 58}, {17, 43}, {18, 52}, {18, 58}, {19, 59}, {19, 45}, {20, 32}, {21, 38}, {22, 41}, {22, 60}, {23, 33}, {23, 55}, {23, 57}, {23, 42}, {24, 55}, {24, 42}, {25, 35}, {25, 58}, {26, 35}, {26, 59}, {27, 51}, {27, 55}, {28, 40}, {29, 57}, {29, 61}, {30, 54}, {30, 47}};


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
    v1 = new Vertex[31];
    v2 = new Vertex[31];

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
    */
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


    BGraph New(31, 31, 63,
              v1, v2, e,
              0);

    cout << "Crossing : " << New.n_cross();
    print("matrix")
    New.printMatrix_2();

    print("TESTE <SADASDADADASDASD")
    print_map(New.map_v1);
    //print("test_deg");

    //New.printGraph();

    print("test_deg_2");

    New.printGraph2();

    print("deg_test");
    unordered_map<int,int> aux;

/*
    unordered_set<int> U_1({1,2,3,4,5});
    unordered_set<int> U_2({6,7,8,9});

    unordered_set<int> V({3,4,7,8});
    unordered_set<int> V_sub({1,2,5,6,9});
*/
    unordered_set<int> U_1({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30});
    unordered_set<int> U_2({31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61});

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
/*
    print("")
    print("greedy");
    for (int i=0;i<9;i++){
        print(New.greedy_selection(1,V,V_sub));
    }*/





    clock_t t_start, t_end;
    aux = New.degrees();


    random_choice(New, aux, 42);

    t_start = clock();
    construction_phase(New, U_1, U_2, 1);
    //New.printBGraph();
    improvement_phase(New,42);
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


    //cout << "Hello world!" << endl;
    //cout << v.get_pos() << endl;
    //cout << pi(v);

    delete v1;
    delete v2;
    return 0;
}
