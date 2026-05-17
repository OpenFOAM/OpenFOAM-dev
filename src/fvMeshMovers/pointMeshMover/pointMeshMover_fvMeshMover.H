/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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
    Foam::fvMeshMovers::pointMeshMover

Description
    Wrapper class so that a pointMeshMover can be instantiated to move an fvMesh

SourceFiles
    pointMeshMover_fvMeshMover.C

\*---------------------------------------------------------------------------*/

#ifndef pointMeshMover_fvMeshMover_H
#define pointMeshMover_fvMeshMover_H

#include "fvMeshMover.H"
#include "pointMeshMover.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{

/*---------------------------------------------------------------------------*\
                        Class pointMeshMover Declaration
\*---------------------------------------------------------------------------*/

class pointMeshMover
:
    public fvMeshMover
{
    // Private Data

        const word pointMeshMoverName_;

        //- Pointer to the run-time selected pointMeshMover
        autoPtr<Foam::pointMeshMover> pointMeshMoverPtr_;


public:

    //- Runtime type information
    TypeName(Foam::pointMeshMover::typeName_());


    // Constructors

        //- Construct from fvMesh
        pointMeshMover(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        pointMeshMover(const pointMeshMover&) = delete;


    //- Destructor
    ~pointMeshMover();


    // Member Functions

        //- Return the pointMeshMover
        const Foam::pointMeshMover& mover() const;

        //- Is this motion solid body? Delegate to the motion solver.
        virtual bool solidBodyMotion() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Write the motion solver state
        virtual bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointMeshMover&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
