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
        int *pos;

        Vertex(){
            this->v = -1;
            this->pos = new int;
            (*this->pos) = 0;
        }

        Vertex(int i, int p){
            this->v = i;
            this->pos = new int;
            (*this->pos) = p;
        }

        //~Vertex(){
        //    delete this->pos;
        //}

        void setVertex(int i, int p){
            this->v = i;
            this->pos = new int;
            (*this->pos) = p;
        }

        void setPos(int p){
            (*this->pos) = p;
        }

        int get_vertex() const {
            return v;
        }

        int get_pos() const {
            return (*pos);
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

/*

class Pi{ // Hash Table

    private:
        int max_items;


        int getHash(Vertex v){
            return v.get_pos() % max_items;
        }


    public(int max_item=1000);



};
 */
typedef unordered_map<int, unordered_map<int,int>> AdjMatrix;

typedef unordered_set<int> Adj;
typedef unordered_map<int, vector<int> > Adj2;

typedef unordered_map<int, Vertex> Dict;
typedef map<int, Vertex> VSorted;


typedef struct Edge{
    int u;
    int v;
} edge;

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
            for (i=1;i<=n_v2_;i++){
                this->map_v2[v2_[i-1].get_pos()] = v2_[i-1];
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

            cout << " Deg size " << deg_dict.size() << endl ;

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
                    cout << max_q.top().first << " ";
                max_idx++;
                max_q.pop();
            }
            cout << endl;
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
                cout << " barycenter" << *u << " : " << (*pi_k)[*u].get_pos() << endl;
                if((*pi_k)[*u].get_pos()){
                    b+=(*pi_k)[*u].get_pos();
                    K++;
                    print(this->pi_2[6].get_pos());
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
                (*pos_assing)[i->second.get_pos()] = i->second;
            }
            this->maped = 0;
        }

        void printBGraph(){
            if (this->maped) { re_map(this->maped); }

            char u;
            char v;
            int n_max = (this->n_v1 > this->n_v2 ? this->n_v1 : this->n_v2);
            for (int i=1;i <= n_max ; i++){
                u = this->map_v1[i].get_vertex() !=-1  ? (char) this->map_v1[i].get_vertex() + '0' : (char) 0 ;
                v = this->map_v2[i].get_vertex() !=-1  ? (char) this->map_v2[i].get_vertex() + '0' : (char) 0 ;
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

//#define pi(v) v.get_pos()
//typedef hash_map<int,int> Pi;V.insert(*i);


//template<typename T>


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

        cout << *i<<"("<<(G.pi_2[*i].get_pos()) << ") : " << endl ;
        U.insert(*i);

        pos = G.pi_2[*i].get_pos();
        G.pi_2[*i].setPos(0);

        G.map_v2[pos] = *(new Vertex());
    }


    for (auto i = U_1.begin(); i != U_1.end(); i++){

        cout << *i<<"("<<(G.pi_1[*i].get_pos()) << ") : " << endl ; //print exclude

        U.insert(*i);

        pos = G.pi_1[*i].get_pos();
        G.pi_1[*i].setPos(0);

        G.map_v1[pos] = *(new Vertex());
    }

    for (auto i = U.begin(); i != U.end(); i++){
        cout << *i << " : ";
        G.pi_1.find(*i) == G.pi_1.end() ? cout << (G.pi_2[*i].get_pos()) << " /" : cout << (G.pi_1[*i].get_pos()) << " /";
        cout << endl;
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
    print("check pos");
    cout << G.map_v1[3].get_pos() << endl;
    cout << G.n_cross() << endl;

    cout << G.pi_2[v].get_pos() << endl;

    U.erase(v);
    V.insert(v);
    while (U.size()){


        cout<<"SIZE V:" << V.size() << endl;
        cout<<"SIZE U:" << U.size() << endl;
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
        cout << "Layer " << k << " Size " << pos_assing->size() << endl;

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

            print("pos_assing _minus");
            cout << bc_v_minus << " : " << ((*pos_assing)[bc_v_minus].get_pos())
                                        << "("<< ((*pos_assing)[bc_v_minus].get_vertex()) << ")"
                                        << endl;

            print("pos_assing _plus");
            cout << bc_v_plus << " : "  << ((*pos_assing)[bc_v_plus].get_pos())
                                        << "(" <<((*pos_assing)[bc_v_plus].get_vertex()) << ")"
                                        << endl;

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

        cout << i->first << "(" << i->second << ") * ";

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
    cout << endl;
    cout << i ;

    return sample;
    //while(Degrees.size()){
    //}

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

    unsigned int n_cross_aux=-1;
    unsigned int n_cross_bc=-1;
    unsigned int n_cross_plus=-1;
    unsigned int n_cross_minus=-1;


    float bc;
    int v_bc;
    int v_plus;
    int v_minus;

    sample = random_choice(G, deg, seed);
    print("")
    print("############SAMPLE CHECK################");
    for (auto i = sample.begin(); i != sample.end(); i++){
        print(*i);

        v = *i;
        k = G.layer[v];
        bc = G.bc(v, k)

        switch(k){
            case 1:
                pi_aux = G.pi_1;
                n_cross_aux = G.n_cross();

                if ( ((int) bc == 1) ){
                    v_bc = bc + 1;
                    v_minus = v_bc - 1;
                    v_plus =  v_bc + 1 > G.n_v2 ? 0 : v_bc + 1 ;
                }else if ( (int) bc == G.n_v2 ) {
                    v_bc = bc;
                    v_minus = v_bc - 1 > 0 ;
                    v_plus =  v_bc + 1 > G.n_v2 ? 0 : v_bc + 1 ;
                }
                int_bc = int_bc > G.n_v2 - 1



            break;

            case 2:
                pi_aux = G.pi_2;
            break;
        }


        G.move_vertex()

    }
}

int main(){

    Edge e[] = {{1,6}, {2,6}, {3,7}, {3,8}, {4,7}, {5,6}, {5,9}, {2,9},{4,6}};


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
    v1 = new Vertex[5];
    v2 = new Vertex[4];

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

    cout << "Crossing : " << New.n_cross();
    print("matrix")
    New.printMatrix_2();


    //print("test_deg");

    //New.printGraph();

    print("test_deg_2");

    New.printGraph2();

    print("deg_test");
    unordered_map<int,int> aux;


    unordered_set<int> U_1({1,2,3,4,5});
    unordered_set<int> U_2({6,7,8,9});

    unordered_set<int> V({3,4,7,8});
    unordered_set<int> V_sub({1,2,5,6,9});

    print("union")
    srand(1);

    print(New.greedy_selection(1, V, V_sub));

    print("################");
    print("MOVE");

    New.printBGraph();
    print("Crossing " << New.n_cross());

    //auto r_test = New.move_vertex(2,1,0);
    New.pi_1 = New.move_vertex(2,5,0);;
    print("Crossing " << New.n_cross());

    print(New.pi_1[1].get_pos());
    print(New.pi_1[2].get_pos());
    print(New.pi_1[3].get_pos());
    print(New.pi_1[4].get_pos());
    print(New.pi_1[5].get_pos());
    New.n_cross();
    //New.printBGraph();
    print("Crossing " << New.n_cross());

    //New.move_vertex(3,4,1);
    New.printBGraph();
    print("Crossing " << New.n_cross());

    New.move_vertex(6,4,1);
    New.printBGraph();
    print("Crossing " << New.n_cross());

    New.move_vertex(3,5,1);
    New.printBGraph();

    New.move_vertex(3,1,1);
    New.printBGraph();

    New.move_vertex(3,2,1);
    New.printBGraph();

    New.move_vertex(3,3,1);
    New.printBGraph();
    print("Crossing " << New.n_cross());


    construction_phase(New, U_1, U_2, 1);
    New.printBGraph();
/*
    print("")
    print("greedy");
    for (int i=0;i<9;i++){
        print(New.greedy_selection(1,V,V_sub));
    }*/





    time_t t_start, t_end;


    aux = New.degrees();

    time(&t_start);
    random_choice(New, aux, 42);

    improvement_phase(New,42);



    time(&t_end);

    print("time")
    double time_taken = double(t_end - t_start);
    cout << "Time taken by program is : " << fixed
         << time_taken << setprecision(10);
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
