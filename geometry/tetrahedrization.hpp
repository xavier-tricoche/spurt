#ifndef __TETRAHEDRIZATION_HPP__
#define __TETRAHEDRIZATION_HPP__

#include <math/types.hpp>

// blatantly and shamelessly copied from FAnToM's code...

namespace spurt {
typedef small_vector<int, 8>  hexahedron;
typedef small_vector<int, 6>  prism;
typedef small_vector<int, 4>  tetrahedron;
typedef small_vector<int, 5>  pyramid;
typedef small_vector<int, 4>  face;

/*
       7--------6
      /|       /|
     / |      / |
    4--------5  |
    |  |     |  |
    |  3-----|--2
    | /      | /
    |/       |/
    0--------1

   Faces:   0=LEFT, 1=RIGHT,  2=FRONT,
   3=BACK, 4=DOWN, 5=UP

*/
static const int hexFiveTetsTable[2][5][4] = {
    {
        { 0, 1, 2, 5,},
        { 0, 2, 3, 7,},
        { 0, 5, 7, 4,},
        { 5, 2, 7, 6,},
        { 0, 5, 2, 7,}
    },
    {
        { 0, 1, 3, 4 },
        { 1, 2, 3, 6 },
        { 1, 6, 4, 5 },
        { 4, 6, 3, 7 },
        { 1, 3, 4, 6 }
    }
};

static const int hexPrismsVertsTable[6][2][6] = {
    {
        {0, 7, 4, 1, 6, 5},
        {0, 3, 7, 1, 2, 6}
    },
    {
        {0, 3, 4, 1, 2, 5},
        {4, 3, 7, 5, 2, 6}
    },
    {
        {0, 5, 1, 3, 6, 2},
        {0, 4, 5, 3, 7, 6}
    },
    {
        {0, 4, 1, 3, 7, 2},
        {1, 4, 5, 2, 7, 6}
    },
    {
        {0, 1, 2, 4, 5, 6},
        {0, 2, 3, 4, 6, 7}
    },
    {
        {0, 1, 3, 4, 5, 7},
        {1, 2, 3, 5, 6, 7}
    }
};

/*

      5
     /|\
    / | \
   3-----4
   |  |  |
   |  |  |
   |  2  |
   | / \ |
   |/   \|
   0-----1

A 0 in the binary index stands for a diagonal from the lowest vertex-index
of the side

 */
static const int prismTetsTable[8][3] = {
    {  2,  4,  7 }, /*000*/
    {  2,  3, 10 }, /*001*/
    {  1,  6,  7 }, /*010*/
    { -1, -1, -1 }, /*011*/   // The cases "100" and "011"
    { -1, -1, -1 }, /*100*/   // dont't appear
    { 10,  0,  9 }, /*101*/
    {  1,  5, 11 }, /*110*/
    {  0,  8, 11 } /*111*/
};

static const int prismTetsVertsTable[12][4] = {
    { 0, 1, 2, 3 }, /*0*/
    { 0, 1, 2, 4 }, /*1*/
    { 0, 1, 2, 5 }, /*2*/
    { 0, 3, 1, 5 }, /*3*/
    { 0, 4, 1, 5 }, /*4*/
    { 0, 2, 3, 4 }, /*5*/
    { 0, 4, 2, 5 }, /*6*/
    { 0, 4, 5, 3 }, /*7*/
    { 1, 2, 3, 4 }, /*8*/
    { 1, 2, 3, 5 }, /*9*/
    { 1, 3, 4, 5 }, /*10*/
    { 2, 3, 4, 5 }, /*11*/
};
}

namespace spurt {
inline bool config(const face& f)
{
    int min = std::numeric_limits<int>::max();;
    int minid = 0;
    for (int i = 0; i < 4; i++) {
        if (f[i] < min) {
            min = f[i];
            minid = i;
        }
    }
    return (minid % 2);
}

void decompose_prism(std::list<tetrahedron>&, const prism&);

void decompose_hexahedron(std::list<tetrahedron>& tets, const hexahedron& H)
{
    bool configuration[6];
    face f;

    f[0] = H[0];
    f[1] = H[4];
    f[2] = H[7];
    f[3] = H[3];
    configuration[0] = config(f);

    f[0] = H[1];
    f[1] = H[2];
    f[2] = H[6];
    f[3] = H[5];
    configuration[1] = config(f);

    f[0] = H[0];
    f[1] = H[1];
    f[2] = H[5];
    f[3] = H[4];
    configuration[2] = config(f);

    f[0] = H[3];
    f[1] = H[7];
    f[2] = H[6];
    f[3] = H[2];
    configuration[3] = config(f);

    f[0] = H[0];
    f[1] = H[3];
    f[2] = H[2];
    f[3] = H[1];
    configuration[4] = config(f);

    f[0] = H[4];
    f[1] = H[5];
    f[2] = H[6];
    f[3] = H[7];
    configuration[5] = config(f);

    prism inVerts0, inVerts1;
    if (configuration[0] == configuration[1]) {
        if (!configuration[0]) {
            for (int i = 0; i < 6; i++) {
                inVerts0[i] = H[hexPrismsVertsTable[0][0][i]];
            }
            decompose_prism(tets, inVerts0);
            for (int i = 0; i < 6; i++) {
                inVerts1[i] = H[hexPrismsVertsTable[0][1][i]];
            }
            decompose_prism(tets, inVerts1);
        } else {
            for (int i = 0; i < 6; i++) {
                inVerts0[i] = H[hexPrismsVertsTable[1][0][i]];
            }
            decompose_prism(tets, inVerts0);
            for (int i = 0; i < 6; i++) {
                inVerts1[i] = H[hexPrismsVertsTable[1][1][i]];
            }
            decompose_prism(tets, inVerts1);
        }
    } else {
        if (configuration[2] == configuration[3]) {
            if (!configuration[2]) {
                for (int i = 0; i < 6; i++) {
                    inVerts0[i] = H[hexPrismsVertsTable[2][0][i]];
                }
                decompose_prism(tets, inVerts0);
                for (int i = 0; i < 6; i++) {
                    inVerts1[i] = H[hexPrismsVertsTable[2][1][i]];
                }
                decompose_prism(tets, inVerts1);
            } else {
                for (int i = 0; i < 6; i++) {
                    inVerts0[i] = H[hexPrismsVertsTable[3][0][i]];
                }
                decompose_prism(tets, inVerts0);
                for (int i = 0; i < 6; i++) {
                    inVerts1[i] = H[hexPrismsVertsTable[3][1][i]];
                }
                decompose_prism(tets, inVerts1);
            }
        } else {
            if (configuration[4] == configuration[5]) {
                if (!configuration[4]) {
                    for (int i = 0; i < 6; i++) {
                        inVerts0[i] = H[hexPrismsVertsTable[4][0][i]];
                    }
                    decompose_prism(tets, inVerts0);
                    for (int i = 0; i < 6; i++) {
                        inVerts1[i] = H[hexPrismsVertsTable[4][1][i]];
                    }
                    decompose_prism(tets, inVerts1);
                } else {
                    for (int i = 0; i < 6; i++) {
                        inVerts0[i] = H[hexPrismsVertsTable[5][0][i]];
                    }
                    decompose_prism(tets, inVerts0);
                    for (int i = 0; i < 6; i++) {
                        inVerts1[i] = H[hexPrismsVertsTable[5][1][i]];
                    }
                    decompose_prism(tets, inVerts1);
                }
            } else {
                face inVerts;
                for (int i = 0; i < 5; i++) {
                    for (int j = 0; j < 4; j++) {
                        inVerts[j] = H[hexFiveTetsTable[configuration[0]][i][j]];
                    }
                    tets.push_back(inVerts);
                }
            }
        }
    }
}

void decompose_prism(std::list<tetrahedron>& tets, const prism& P)
{
    int configuration = 0;
    face f;
    f[0] = P[0];
    f[1] = P[1];
    f[2] = P[4];
    f[3] = P[3];
    if (config(f)) {
        configuration += 1;
    }
    f[0] = P[1];
    f[1] = P[2];
    f[2] = P[5];
    f[3] = P[4];
    if (config(f)) {
        configuration += 2;
    }
    f[0] = P[0];
    f[1] = P[3];
    f[2] = P[5];
    f[3] = P[2];
    if (config(f)) {
        configuration += 4;
    }
    tetrahedron inVerts0, inVerts1, inVerts2;
    bool equal[3];

    inVerts0[0] =
        P[prismTetsVertsTable[prismTetsTable[configuration][0]][0]];
    inVerts0[1] =
        P[prismTetsVertsTable[prismTetsTable[configuration][0]][1]];
    inVerts0[2] =
        P[prismTetsVertsTable[prismTetsTable[configuration][0]][2]];
    inVerts0[3] =
        P[prismTetsVertsTable[prismTetsTable[configuration][0]][3]];
    tets.push_back(inVerts0);

    inVerts1[0] =
        P[prismTetsVertsTable[prismTetsTable[configuration][1]][0]];
    inVerts1[1] =
        P[prismTetsVertsTable[prismTetsTable[configuration][1]][1]];
    inVerts1[2] =
        P[prismTetsVertsTable[prismTetsTable[configuration][1]][2]];
    inVerts1[3] =
        P[prismTetsVertsTable[prismTetsTable[configuration][1]][3]];
    tets.push_back(inVerts1);

    inVerts2[0] =
        P[prismTetsVertsTable[prismTetsTable[configuration][2]][0]];
    inVerts2[1] =
        P[prismTetsVertsTable[prismTetsTable[configuration][2]][1]];
    inVerts2[2] =
        P[prismTetsVertsTable[prismTetsTable[configuration][2]][2]];
    inVerts2[3] =
        P[prismTetsVertsTable[prismTetsTable[configuration][2]][3]];
    tets.push_back(inVerts2);
}

void decompose_pyramid(std::list<tetrahedron>& tets, const pyramid& P)
{
    face f;
    bool equal0, equal1;
    f[0] = P[0];
    f[1] = P[1];
    f[2] = P[2];
    f[3] = P[3];

    tetrahedron inVerts0, inVerts1;

    if (!config(f)) {
        inVerts0[0]=P[0];
	    inVerts0[1]=P[1];
	    inVerts0[2]=P[2];
	    inVerts0[3]=P[4];
        tets.push_back(inVerts0);

    	inVerts1[0]=P[0];
    	inVerts1[1]=P[2];
    	inVerts1[2]=P[3];
    	inVerts1[3]=P[4];
    	tets.push_back(inVerts1);
    }
    else {
    	inVerts0[0]=P[0];
    	inVerts0[1]=P[1];
    	inVerts0[2]=P[3];
    	inVerts0[3]=P[4];
    	tets.push_back(inVerts0);

    	inVerts1[0]=P[1];
    	inVerts1[1]=P[2];
    	inVerts1[2]=P[3];
    	inVerts1[3]=P[4];
    	tets.push_back(inVerts1);
    }
}
}


#endif
