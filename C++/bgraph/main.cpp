#include <iostream>
#include <limits>
#include "bgraph.h"

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;

#define getrandom(min, max) \
    ((rand()%(int)(((max) + 1)-(min)))+ (min))

typedef struct{
    int *pi_1;            /* permutation of the first part vertices  */
    int *pi_2;            /* permutation of the second part vertices */
    long n_crossing;      /* solution value                          */
    double total_time;    /* total time, secs                        */
    double time_to_best;  /* time to the best solution, secs         */
    int characts[10];     /* some characteristics:                   */
                        /*   characts[0] - size of the first part  */
                        /*   characts[1] - size of the second part */
                        /*   characts[2] - number of edges         */
 }  Result;

//void grasp(Result *);


#define print(txt) cout << txt << endl;

int main(){
    cout << getrandom(10,11) << endl;

    Vertex v(1, 6);

    Vertex *pi_1[5];

    pi_1[0] = &v;

    pi_1[1] = new Vertex(5, 6);

    Vertex *v1;//v1[5];
    Vertex *v2;//v2[4];
    v1 = new Vertex[5];
    v2 = new Vertex[4];

    print("test_dadasmain");
    for (int i=0;i<5;i++){
        cout << "v -> " << v1[i].get_vertex() << " pos " << v1[i].get_pos() <<endl;
    }
     print("test_main");
    for (int i=0;i<4;i++){
        cout << "v -> " << v2[i].get_vertex() << " pos " << v2[i].get_pos() <<endl;
    }

    v1[0] = Vertex(1, 1);
    v1[1] = Vertex(2, 2);
    v1[2] = Vertex(3, 3);
    v1[3] = Vertex(4, 4);
    v1[4] = Vertex(5, 5);

    print("test_main");
    for (int i=0;i<5;i++){
        cout << "v -> " << v1[i].get_vertex() << " pos " << v1[i].get_pos() <<endl;
    }


    v2[0] = Vertex(6, 1);
    v2[1] = Vertex(7, 2);
    v2[2] = Vertex(8, 3);
    v2[3] = Vertex(9, 4);

    print("test_main");
    for (int i=0;i<4;i++){
        cout << "v -> " << v2[i].get_vertex() << " pos " << v2[i].get_pos() <<endl;
    }


    cout << "Hello world!" << endl;
    cout << v2[0].get_pos() << endl;
    cout << v2[0].get_vertex();


   return 0;
}
