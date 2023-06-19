/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::fvMeshDistributor

Description
    Abstract base class for fvMesh movers.

    These classes move the mesh points, update the cell volumes and generate
    the corresponding mesh fluxes without any topology change.

SourceFiles
    fvMeshDistributor.C
    fvMeshDistributorNew.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshDistributor_H
#define fvMeshDistributor_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fvMeshDistributor Declaration
\*---------------------------------------------------------------------------*/

class fvMeshDistributor
{
    // Private Data

        //- Non-const fvMesh reference to allow update
        fvMesh& mesh_;

        //- The dynamicMeshDict
        dictionary dynamicMeshDict_;


public:

    //- Runtime type information
    TypeName("fvMeshDistributor");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvMeshDistributor,
            fvMesh,
            (fvMesh& mesh),
            (mesh)
        );


    //- Helper class to update the velocity boundary conditions
    //  following mesh motion
    class velocityMotionCorrection
    {
        // Private Data

            const fvMesh& mesh_;

            wordList velocityFields_;

    public:

        // Constructors

            velocityMotionCorrection
            (
                const fvMesh& mesh,
                const dictionary& dict
            );


        // Member Functions

            void update() const;
    };


    // Constructors

        //- Construct from fvMesh
        explicit fvMeshDistributor(fvMesh&);

        //- Disallow default bitwise copy construction
        fvMeshDistributor(const fvMeshDistributor&) = delete;


    // Selectors

        //- Select, construct and return the fvMeshDistributor
        //  If the constant/dynamicMeshDict does not exist
        //  a staticFvMesh is returned
        static autoPtr<fvMeshDistributor> New(fvMesh&);


    //- Destructor
    virtual ~fvMeshDistributor();


    // Member Functions

        //- Return the fvMesh
        fvMesh& mesh()
        {
            return mesh_;
        }

        //- Return the fvMesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the dynamicMeshDict/distributor sub-dict
        const dictionary& dict() const
        {
            return dynamicMeshDict_.subDict("distributor");
        }

        //- Is mesh dynamic, i.e. might it change?
        //  Defaults to true, set to false in the fvMeshDistributors::none
        virtual bool dynamic() const
        {
            return true;
        }

        //- Update the mesh
        virtual bool update() = 0;

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&) = 0;

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&) = 0;

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&) = 0;

        //- Write the mover state
        virtual bool write(const bool write = true) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshDistributor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
