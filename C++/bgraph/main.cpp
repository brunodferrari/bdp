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

int main(){

    Vertex *v1;
    Vertex *v2;

    Edge e2[] = {{0, 185}, {0, 186}, {1, 187}, {1, 188}, {2, 189}, {3, 190}, {4, 189}, {5, 191}, {6, 192}, {7, 193}, {8, 194}, {9, 195}, {10, 196}, {11, 197}, {11, 198}, {12, 192}, {13, 199}, {13, 200}, {14, 185}, {15, 187}, {16, 201}, {17, 198}, {18, 202}, {19, 199}, {20, 203}, {21, 185}, {21, 203}, {22, 191}, {22, 198}, {23, 202}, {23, 204}, {24, 205}, {25, 206}, {26, 186}, {27, 207}, {27, 208}, {28, 196}, {29, 193}, {30, 209}, {31, 210}, {32, 211}, {33, 212}, {34, 213}, {35, 214}, {36, 202}, {37, 208}, {38, 192}, {38, 193}, {39, 215}, {40, 211}, {41, 210}, {42, 216}, {43, 186}, {44, 205}, {45, 217}, {46, 206}, {47, 218}, {47, 192}, {48, 219}, {49, 220}, {50, 221}, {51, 208}, {52, 222}, {52, 208}, {53, 223}, {54, 205}, {55, 219}, {55, 195}, {56, 201}, {56, 200}, {57, 211}, {58, 186}, {59, 218}, {59, 224}, {60, 186}, {60, 225}, {61, 194}, {62, 220}, {63, 203}, {64, 194}, {64, 212}, {65, 204}, {66, 212}, {66, 214}, {67, 197}, {68, 200}, {69, 204}, {70, 200}, {71, 209}, {72, 193}, {72, 211}, {73, 201}, {73, 194}, {74, 226}, {75, 200}, {76, 222}, {77, 190}, {77, 210}, {78, 223}, {79, 224}, {80, 216}, {81, 218}, {82, 225}, {83, 210}, {84, 188}, {85, 221}, {86, 212}, {87, 214}, {88, 226}, {88, 195}, {89, 207}, {90, 195}, {91, 220}, {92, 219}, {93, 187}, {93, 207}, {94, 226}, {95, 202}, {95, 225}, {96, 224}, {97, 215}, {98, 197}, {98, 199}, {99, 223}, {100, 189}, {100, 216}, {101, 222}, {102, 206}, {102, 227}, {103, 198}, {104, 217}, {104, 221}, {105, 209}, {106, 228}, {107, 190}, {108, 225}, {109, 191}, {110, 213}, {111, 227}, {112, 227}, {112, 215}, {113, 187}, {114, 223}, {114, 214}, {115, 185}, {116, 229}, {117, 207}, {118, 220}, {118, 205}, {119, 185}, {120, 225}, {121, 218}, {122, 219}, {122, 228}, {123, 192}, {124, 188}, {125, 209}, {125, 228}, {126, 194}, {127, 202}, {128, 197}, {129, 204}, {129, 230}, {130, 196}, {131, 222}, {132, 221}, {133, 229}, {134, 189}, {134, 210}, {135, 203}, {136, 201}, {137, 224}, {138, 226}, {138, 229}, {139, 224}, {139, 216}, {140, 227}, {141, 199}, {142, 230}, {143, 213}, {144, 214}, {145, 215}, {146, 199}, {147, 230}, {148, 206}, {149, 213}, {150, 217}, {151, 203}, {151, 211}, {152, 188}, {153, 221}, {153, 206}, {154, 193}, {155, 228}, {156, 226}, {157, 191}, {157, 196}, {158, 201}, {159, 229}, {160, 216}, {161, 198}, {162, 213}, {162, 230}, {163, 219}, {164, 188}, {164, 229}, {165, 207}
};
    BGraph *New;
    v1 = new Vertex[166];
    v2 = new Vertex[46];

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


     New = new BGraph(166, 46, 207,
                        v1, v2, e2);


    print(New->n_cross());

    Edge e[] = {{1,6}, {2,6}, {3,7}, {3,8}, {4,7}, {5,6}, {5,9}, {2,9},{4,6}};

    delete[] v1;
    delete[] v2;

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


    /*(5, 4, 7,
              v1, v2, e);*/

    New = new BGraph(5, 4, 7,
              v1, v2, e);

    New->printBGraph();

    print(New->n_cross());
    cout << "aux";


    return 0;
}
