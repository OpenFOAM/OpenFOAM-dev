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

Description
    SAMM cell shape lookup table

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sammMesh::fillSammCellShapeTable()
{
    // Fill the list by hand

    // SAMM trim type 1: 8 models
    sammShapeLookup[1]   = sammTrim1Ptr_;
    sammShapeLookup[2]   = sammTrim1Ptr_;
    sammShapeLookup[4]   = sammTrim1Ptr_;
    sammShapeLookup[8]   = sammTrim1Ptr_;
    sammShapeLookup[16]  = sammTrim1Ptr_;
    sammShapeLookup[32]  = sammTrim1Ptr_;
    sammShapeLookup[64]  = sammTrim1Ptr_;
    sammShapeLookup[128] = sammTrim1Ptr_;

    // SAMM trim type 2: 12 models
    sammShapeLookup[3]   = sammTrim2Ptr_;
    sammShapeLookup[12]  = sammTrim2Ptr_;
    sammShapeLookup[192] = sammTrim2Ptr_;
    sammShapeLookup[48]  = sammTrim2Ptr_;
    sammShapeLookup[9]   = sammTrim2Ptr_;
    sammShapeLookup[144] = sammTrim2Ptr_;
    sammShapeLookup[96]  = sammTrim2Ptr_;
    sammShapeLookup[6]   = sammTrim2Ptr_;
    sammShapeLookup[17]  = sammTrim2Ptr_;
    sammShapeLookup[34]  = sammTrim2Ptr_;
    sammShapeLookup[68]  = sammTrim2Ptr_;
    sammShapeLookup[136] = sammTrim2Ptr_;

    // SAMM trim type 3: 24 models
    sammShapeLookup[7]   = sammTrim3Ptr_;
    sammShapeLookup[14]  = sammTrim3Ptr_;
    sammShapeLookup[13]  = sammTrim3Ptr_;
    sammShapeLookup[11]  = sammTrim3Ptr_;
    sammShapeLookup[112] = sammTrim3Ptr_;
    sammShapeLookup[224] = sammTrim3Ptr_;
    sammShapeLookup[208] = sammTrim3Ptr_;
    sammShapeLookup[176] = sammTrim3Ptr_;
    sammShapeLookup[38]  = sammTrim3Ptr_;
    sammShapeLookup[70]  = sammTrim3Ptr_;
    sammShapeLookup[100] = sammTrim3Ptr_;
    sammShapeLookup[98]  = sammTrim3Ptr_;
    sammShapeLookup[25]  = sammTrim3Ptr_;
    sammShapeLookup[137] = sammTrim3Ptr_;
    sammShapeLookup[152] = sammTrim3Ptr_;
    sammShapeLookup[145] = sammTrim3Ptr_;
    sammShapeLookup[49]  = sammTrim3Ptr_;
    sammShapeLookup[50]  = sammTrim3Ptr_;
    sammShapeLookup[35]  = sammTrim3Ptr_;
    sammShapeLookup[19]  = sammTrim3Ptr_;
    sammShapeLookup[200] = sammTrim3Ptr_;
    sammShapeLookup[196] = sammTrim3Ptr_;
    sammShapeLookup[76]  = sammTrim3Ptr_;
    sammShapeLookup[140] = sammTrim3Ptr_;

    // SAMM trim type 4: 8 models
    sammShapeLookup[27]  = sammTrim4Ptr_;
    sammShapeLookup[39]  = sammTrim4Ptr_;
    sammShapeLookup[78]  = sammTrim4Ptr_;
    sammShapeLookup[141] = sammTrim4Ptr_;
    sammShapeLookup[177] = sammTrim4Ptr_;
    sammShapeLookup[114] = sammTrim4Ptr_;
    sammShapeLookup[228] = sammTrim4Ptr_;
    sammShapeLookup[216] = sammTrim4Ptr_;

    // SAMM trim type 5: 24 models
    sammShapeLookup[248] = sammTrim5Ptr_;
    sammShapeLookup[241] = sammTrim5Ptr_;
    sammShapeLookup[242] = sammTrim5Ptr_;
    sammShapeLookup[244] = sammTrim5Ptr_;
    sammShapeLookup[143] = sammTrim5Ptr_;
    sammShapeLookup[31]  = sammTrim5Ptr_;
    sammShapeLookup[47]  = sammTrim5Ptr_;
    sammShapeLookup[79]  = sammTrim5Ptr_;
    sammShapeLookup[217] = sammTrim5Ptr_;
    sammShapeLookup[185] = sammTrim5Ptr_;
    sammShapeLookup[155] = sammTrim5Ptr_;
    sammShapeLookup[157] = sammTrim5Ptr_;
    sammShapeLookup[230] = sammTrim5Ptr_;
    sammShapeLookup[118] = sammTrim5Ptr_;
    sammShapeLookup[103] = sammTrim5Ptr_;
    sammShapeLookup[110] = sammTrim5Ptr_;
    sammShapeLookup[206] = sammTrim5Ptr_;
    sammShapeLookup[205] = sammTrim5Ptr_;
    sammShapeLookup[220] = sammTrim5Ptr_;
    sammShapeLookup[236] = sammTrim5Ptr_;
    sammShapeLookup[55]  = sammTrim5Ptr_;
    sammShapeLookup[59]  = sammTrim5Ptr_;
    sammShapeLookup[179] = sammTrim5Ptr_;
    sammShapeLookup[115] = sammTrim5Ptr_;

    // SAMM hexagonal prism (trim type 8): 1 model
    sammShapeLookup[255] = sammTrim8Ptr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
