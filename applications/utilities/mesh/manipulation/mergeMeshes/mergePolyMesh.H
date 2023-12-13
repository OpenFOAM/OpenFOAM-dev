/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Class
    Foam::mergePolyMesh

Description
    Add a given mesh to the original mesh to create a single new mesh

SourceFiles
    mergePolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mergePolyMesh_H
#define mergePolyMesh_H

#include "polyMesh.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class mergePolyMesh Declaration
\*---------------------------------------------------------------------------*/

class mergePolyMesh
:
    public polyMesh
{
    // Private Data

        //- Topological change to accumulated all mesh changes
        polyTopoChange meshMod_;

        //- Patch names
        DynamicList<word> patchNames_;

        //- Patch dictionaries
        DynamicList<dictionary> patchDicts_;

        //- Point zone names
        DynamicList<word> pointZoneNames_;

        //- Face zone names
        DynamicList<word> faceZoneNames_;

        //- Cell zone names
        DynamicList<word> cellZoneNames_;


    // Private Member Functions

        //- Return patch index given a name and type
        label patchIndex(const polyPatch&);

        //- Return zone index given a list of active zones and a name
        label zoneIndex(DynamicList<word>&, const word&);


public:

    //- Runtime type information
    TypeName("mergePolyMesh");


    // Constructors

        //- Construct from IOobject
        mergePolyMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        mergePolyMesh(const mergePolyMesh&) = delete;


    //- Destructor
    virtual ~mergePolyMesh()
    {}


    // Member Functions

        //- Add a mesh
        void addMesh(const polyMesh& m);

        //- Merge meshes
        void merge();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mergePolyMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
