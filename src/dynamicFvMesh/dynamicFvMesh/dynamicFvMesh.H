/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::dynamicFvMesh

Description
    Abstract base class for geometry and/or topology changing fvMesh.

SourceFiles
    dynamicFvMesh.C
    dynamicFvMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicFvMesh_H
#define dynamicFvMesh_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicFvMesh
:
    public fvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicFvMesh(const dynamicFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicFvMesh");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicFvMesh,
            IOobject,
            (const IOobject& io),
            (io)
        );


    // Constructors

        //- Construct from objectRegistry, and read/write options
        explicit dynamicFvMesh(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        dynamicFvMesh
        (
            const IOobject& io,
            const Xfer<pointField>& points,
            const Xfer<faceList>& faces,
            const Xfer<labelList>& allOwner,
            const Xfer<labelList>& allNeighbour,
            const bool syncPar = true
        );

        //- Construct without boundary from cells rather than owner/neighbour.
        //  Boundary is added using addPatches() member function
        dynamicFvMesh
        (
            const IOobject& io,
            const Xfer<pointField>& points,
            const Xfer<faceList>& faces,
            const Xfer<cellList>& cells,
            const bool syncPar = true
        );


    // Selectors

        //- Select null constructed
        static autoPtr<dynamicFvMesh> New(const IOobject& io);


    //- Destructor
    virtual ~dynamicFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
