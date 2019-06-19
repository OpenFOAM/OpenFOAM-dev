/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    // Private Data

        //- The dynamicMeshDict
        dictionary dynamicMeshDict_;


    // Private Member Functions

        //- Return the dynamicMeshDict IOobject
        static IOobject dynamicMeshDictIOobject(const IOobject& io);


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


    //- Helper class to update the velocity boundary conditions
    //  following mesh motion
    class velocityMotionCorrection
    {
        // Private Data

            const dynamicFvMesh& mesh_;

            wordList velocityFields_;

    public:

        // Constructors

            velocityMotionCorrection
            (
                const dynamicFvMesh& mesh,
                const dictionary& dict
            );


        // Member Functions

            void update() const;
    };


    // Constructors

        //- Construct from objectRegistry, and read/write options
        explicit dynamicFvMesh(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        dynamicFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );

        //- Construct without boundary from cells rather than owner/neighbour.
        //  Boundary is added using addPatches() member function
        dynamicFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true
        );

        //- Disallow default bitwise copy construction
        dynamicFvMesh(const dynamicFvMesh&) = delete;


    // Selectors

        //- Select, construct and return the dynamicFvMesh
        //  If the constant/dynamicMeshDict does not exist
        //  a staticFvMesh is returned
        static autoPtr<dynamicFvMesh> New(const IOobject& io);


    //- Destructor
    virtual ~dynamicFvMesh();


    // Member Functions

        //- Return the dynamicMeshDict
        const dictionary& dynamicMeshDict() const
        {
            return dynamicMeshDict_;
        }

        //- Is mesh dynamic
        virtual bool dynamic() const
        {
            return true;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
