/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::faceZoneList

\*---------------------------------------------------------------------------*/

#ifndef faceZoneList_H
#define faceZoneList_H

#include "faceZone.H"
#include "ZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                          Class cellZones Declaration
\*---------------------------------------------------------------------------*/

class faceZoneList
:
    public ZoneList<faceZone, faceZoneList, polyMesh>
{
public:

    //- Runtime type information
    TypeName("faceZoneList");


    // Constructors

        using ZoneList<faceZone, faceZoneList, polyMesh>::ZoneList;


    // Member Functions

        //- Return number of faces in the mesh
        label allSize() const;

        //- Insert given oriented indices
        void insert
        (
            const List<Map<bool>>& zonesOrientedIndices
        );

        //- Insert given unoriented indices
        void insert
        (
            const List<labelHashSet>& zonesIndices
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
