/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sammMesh::fillSammAddressingTable()
{
    // SAMM trim type 1: 8 models
    static label SammTrim1Rot0[10] = {1, 5, 6, 2, 8, 10, 4, 7, 3, 9};
    static label SammTrim1Rot1[10] = {2, 6, 7, 3, 8, 10, 5, 4, 0, 9};
    static label SammTrim1Rot2[10] = {3, 7, 4, 0, 8, 10, 6, 5, 1, 9};
    static label SammTrim1Rot3[10] = {0, 4, 5, 1, 8, 10, 7, 6, 2, 9};
    static label SammTrim1Rot4[10] = {7, 3, 2, 6, 8, 10, 0, 1, 5, 9};
    static label SammTrim1Rot5[10] = {4, 0, 3, 7, 8, 10, 1, 2, 6, 9};
    static label SammTrim1Rot6[10] = {5, 1, 0, 4, 8, 10, 2, 3, 7, 9};
    static label SammTrim1Rot7[10] = {6, 2, 1, 5, 8, 10, 3, 0, 4, 9};

    sammAddressingTable[1]   = SammTrim1Rot0;
    sammAddressingTable[2]   = SammTrim1Rot1;
    sammAddressingTable[4]   = SammTrim1Rot2;
    sammAddressingTable[8]   = SammTrim1Rot3;
    sammAddressingTable[16]  = SammTrim1Rot4;
    sammAddressingTable[32]  = SammTrim1Rot5;
    sammAddressingTable[64]  = SammTrim1Rot6;
    sammAddressingTable[128] = SammTrim1Rot7;


    // SAMM trim type 2: 12 models
    static label SammTrim2Rot0[10]  = {9, 3, 7, 4, 10, 8, 2, 6, 5, 11};
    static label SammTrim2Rot1[10]  = {9, 1, 5, 6, 10, 8, 0, 4, 7, 11};
    static label SammTrim2Rot2[10]  = {9, 2, 1, 5, 10, 8, 3, 0, 4, 11};
    static label SammTrim2Rot3[10]  = {9, 0, 3, 7, 10, 8, 1, 2, 6, 11};

    static label SammTrim2Rot4[10]  = {9, 4, 5, 1, 10, 8, 7, 6, 2, 11};
    static label SammTrim2Rot5[10]  = {9, 5, 1, 0, 10, 8, 6, 2, 3, 11};
    static label SammTrim2Rot6[10]  = {9, 1, 0, 4, 10, 8, 2, 3, 7, 11};
    static label SammTrim2Rot7[10]  = {9, 0, 4, 5, 10, 8, 3, 7, 6, 11};

    static label SammTrim2Rot8[10]  = {9, 1, 2, 3, 10, 8, 5, 6, 7, 11};
    static label SammTrim2Rot9[10]  = {9, 2, 3, 0, 10, 8, 6, 7, 4, 11};
    static label SammTrim2Rot10[10] = {9, 3, 0, 1, 10, 8, 7, 4, 5, 11};
    static label SammTrim2Rot11[10] = {9, 0, 1, 2, 10, 8, 4, 5, 6, 11};

    sammAddressingTable[3]   = SammTrim2Rot0;
    sammAddressingTable[12]  = SammTrim2Rot1;
    sammAddressingTable[192] = SammTrim2Rot2;
    sammAddressingTable[48]  = SammTrim2Rot3;
    sammAddressingTable[9]   = SammTrim2Rot4;
    sammAddressingTable[144] = SammTrim2Rot5;
    sammAddressingTable[96]  = SammTrim2Rot6;
    sammAddressingTable[6]   = SammTrim2Rot7;
    sammAddressingTable[17]  = SammTrim2Rot8;
    sammAddressingTable[34]  = SammTrim2Rot9;
    sammAddressingTable[68]  = SammTrim2Rot10;
    sammAddressingTable[136] = SammTrim2Rot11;


    // SAMM trim type 3: 24 models
    static label SammTrim3Rot0[10]  = {5, 4, 7, 6, 11, 10, 9, 3, 8, 12};
    static label SammTrim3Rot1[10]  = {6, 5, 4, 7, 11, 10, 9, 0, 8, 12};
    static label SammTrim3Rot2[10]  = {7, 6, 5, 4, 11, 10, 9, 1, 8, 12};
    static label SammTrim3Rot3[10]  = {4, 7, 6, 5, 11, 10, 9, 2, 8, 12};
    static label SammTrim3Rot4[10]  = {1, 2, 3, 0, 11, 10, 9, 7, 8, 12};
    static label SammTrim3Rot5[10]  = {2, 3, 0, 1, 11, 10, 9, 4, 8, 12};
    static label SammTrim3Rot6[10]  = {3, 0, 1, 2, 11, 10, 9, 5, 8, 12};
    static label SammTrim3Rot7[10]  = {0, 1, 2, 3, 11, 10, 9, 6, 8, 12};
    static label SammTrim3Rot8[10]  = {0, 3, 7, 4, 11, 10, 9, 6, 8, 12};
    static label SammTrim3Rot9[10]  = {3, 7, 4, 0, 11, 10, 9, 5, 8, 12};
    static label SammTrim3Rot10[10] = {7, 4, 0, 3, 11, 10, 9, 1, 8, 12};
    static label SammTrim3Rot11[10] = {4, 0, 3, 7, 11, 10, 9, 2, 8, 12};
    static label SammTrim3Rot12[10] = {1, 5, 6, 2, 11, 10, 9, 7, 8, 12};
    static label SammTrim3Rot13[10] = {2, 1, 5, 6, 11, 10, 9, 4, 8, 12};
    static label SammTrim3Rot14[10] = {6, 2, 1, 5, 11, 10, 9, 0, 8, 12};
    static label SammTrim3Rot15[10] = {5, 6, 1, 2, 11, 10, 9, 3, 8, 12};
    static label SammTrim3Rot16[10] = {7, 3, 2, 6, 11, 10, 9, 1, 8, 12};
    static label SammTrim3Rot17[10] = {6, 7, 3, 2, 11, 10, 9, 0, 8, 12};
    static label SammTrim3Rot18[10] = {2, 6, 7, 3, 11, 10, 9, 4, 8, 12};
    static label SammTrim3Rot19[10] = {3, 2, 6, 7, 11, 10, 9, 5, 8, 12};
    static label SammTrim3Rot20[10] = {4, 5, 1, 0, 11, 10, 9, 2, 8, 12};
    static label SammTrim3Rot21[10] = {5, 1, 0, 4, 11, 10, 9, 3, 8, 12};
    static label SammTrim3Rot22[10] = {1, 0, 4, 5, 11, 10, 9, 7, 8, 12};
    static label SammTrim3Rot23[10] = {0, 4, 5, 1, 11, 10, 9, 6, 8, 12};

    sammAddressingTable[7]   = SammTrim3Rot0;
    sammAddressingTable[14]  = SammTrim3Rot1;
    sammAddressingTable[13]  = SammTrim3Rot2;
    sammAddressingTable[11]  = SammTrim3Rot3;
    sammAddressingTable[112] = SammTrim3Rot4;
    sammAddressingTable[224] = SammTrim3Rot5;
    sammAddressingTable[208] = SammTrim3Rot6;
    sammAddressingTable[176] = SammTrim3Rot7;
    sammAddressingTable[38]  = SammTrim3Rot8;
    sammAddressingTable[70]  = SammTrim3Rot9;
    sammAddressingTable[100] = SammTrim3Rot10;
    sammAddressingTable[98]  = SammTrim3Rot11;
    sammAddressingTable[25]  = SammTrim3Rot12;
    sammAddressingTable[137] = SammTrim3Rot13;
    sammAddressingTable[152] = SammTrim3Rot14;
    sammAddressingTable[145] = SammTrim3Rot15;
    sammAddressingTable[49]  = SammTrim3Rot16;
    sammAddressingTable[50]  = SammTrim3Rot17;
    sammAddressingTable[35]  = SammTrim3Rot18;
    sammAddressingTable[19]  = SammTrim3Rot19;
    sammAddressingTable[200] = SammTrim3Rot20;
    sammAddressingTable[196] = SammTrim3Rot21;
    sammAddressingTable[76]  = SammTrim3Rot22;
    sammAddressingTable[140] = SammTrim3Rot23;


    // SAMM trim type 4: 8 models
    static label SammTrim4Rot0[10] = {6, 7, 2, 5, 13, 12 ,11, 10, 9, 8};
    static label SammTrim4Rot1[10] = {7, 4, 3, 6, 13, 12 ,11, 10, 9, 8};
    static label SammTrim4Rot2[10] = {4, 5, 6, 7, 13, 12 ,11, 10, 9, 8};
    static label SammTrim4Rot3[10] = {5, 6, 1, 4, 13, 12 ,11, 10, 9, 8};
    static label SammTrim4Rot4[10] = {2, 1, 6, 3, 13, 12 ,11, 10, 9, 8};
    static label SammTrim4Rot5[10] = {3, 2, 7, 0, 13, 12 ,11, 10, 9, 8};
    static label SammTrim4Rot6[10] = {0, 3, 4, 1, 13, 12 ,11, 10, 9, 8};
    static label SammTrim4Rot7[10] = {1, 0, 5, 2, 13, 12 ,11, 10, 9, 8};

    sammAddressingTable[27]  = SammTrim4Rot0;
    sammAddressingTable[39]  = SammTrim4Rot1;
    sammAddressingTable[78]  = SammTrim4Rot2;
    sammAddressingTable[141] = SammTrim4Rot3;
    sammAddressingTable[177] = SammTrim4Rot4;
    sammAddressingTable[114] = SammTrim4Rot5;
    sammAddressingTable[228] = SammTrim4Rot6;
    sammAddressingTable[216] = SammTrim4Rot7;


    // SAMM trim type 5: 24 models
    static label SammTrim5Rot0[8]  = {12, 0, 1, 2, 8, 11, 10, 9};
    static label SammTrim5Rot1[8]  = {12, 1, 2, 3, 8, 11, 10, 9};
    static label SammTrim5Rot2[8]  = {12, 2, 3, 0, 8, 11, 10, 9};
    static label SammTrim5Rot3[8]  = {12, 3, 0, 1, 8, 11, 10, 9};
    static label SammTrim5Rot4[8]  = {12, 6, 5, 4, 8, 11, 10, 9};
    static label SammTrim5Rot5[8]  = {12, 7, 6, 5, 8, 11, 10, 9};
    static label SammTrim5Rot6[8]  = {12, 4, 7, 6, 8, 11, 10, 9};
    static label SammTrim5Rot7[8]  = {12, 5, 4, 7, 8, 11, 10, 9};
    static label SammTrim5Rot8[8]  = {12, 2, 1, 5, 8, 11, 10, 9};
    static label SammTrim5Rot9[8]  = {12, 6, 2, 1, 8, 11, 10, 9};
    static label SammTrim5Rot10[8] = {12, 5, 6, 2, 8, 11, 10, 9};
    static label SammTrim5Rot11[8] = {12, 1, 5, 6, 8, 11, 10, 9};
    static label SammTrim5Rot12[8] = {12, 4, 0, 3, 8, 11, 10, 9};
    static label SammTrim5Rot13[8] = {12, 0, 3, 7, 8, 11, 10, 9};
    static label SammTrim5Rot14[8] = {12, 3, 7, 4, 8, 11, 10, 9};
    static label SammTrim5Rot15[8] = {12, 7, 4, 0, 8, 11, 10, 9};
    static label SammTrim5Rot16[8] = {12, 0, 4, 5, 8, 11, 10, 9};
    static label SammTrim5Rot17[8] = {12, 4, 5, 1, 8, 11, 10, 9};
    static label SammTrim5Rot18[8] = {12, 5, 1, 0, 8, 11, 10, 9};
    static label SammTrim5Rot19[8] = {12, 1, 0, 4, 8, 11, 10, 9};
    static label SammTrim5Rot20[8] = {12, 6, 7, 3, 8, 11, 10, 9};
    static label SammTrim5Rot21[8] = {12, 2, 6, 7, 8, 11, 10, 9};
    static label SammTrim5Rot22[8] = {12, 3, 2, 6, 8, 11, 10, 9};
    static label SammTrim5Rot23[8] = {12, 7, 3, 2, 8, 11, 10, 9};

    sammAddressingTable[248] = SammTrim5Rot0;
    sammAddressingTable[241] = SammTrim5Rot1;
    sammAddressingTable[242] = SammTrim5Rot2;
    sammAddressingTable[244] = SammTrim5Rot3;
    sammAddressingTable[143] = SammTrim5Rot4;
    sammAddressingTable[31]  = SammTrim5Rot5;
    sammAddressingTable[47]  = SammTrim5Rot6;
    sammAddressingTable[79]  = SammTrim5Rot7;
    sammAddressingTable[217] = SammTrim5Rot8;
    sammAddressingTable[185] = SammTrim5Rot9;
    sammAddressingTable[155] = SammTrim5Rot10;
    sammAddressingTable[157] = SammTrim5Rot11;
    sammAddressingTable[230] = SammTrim5Rot12;
    sammAddressingTable[118] = SammTrim5Rot13;
    sammAddressingTable[103] = SammTrim5Rot14;
    sammAddressingTable[110] = SammTrim5Rot15;
    sammAddressingTable[206] = SammTrim5Rot16;
    sammAddressingTable[205] = SammTrim5Rot17;
    sammAddressingTable[220] = SammTrim5Rot18;
    sammAddressingTable[236] = SammTrim5Rot19;
    sammAddressingTable[55]  = SammTrim5Rot20;
    sammAddressingTable[59]  = SammTrim5Rot21;
    sammAddressingTable[179] = SammTrim5Rot22;
    sammAddressingTable[115] = SammTrim5Rot23;


    // SAMM trim type 8: 1 model
    static label SammTrim8[12] = {8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

    sammAddressingTable[255]  = SammTrim8;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
