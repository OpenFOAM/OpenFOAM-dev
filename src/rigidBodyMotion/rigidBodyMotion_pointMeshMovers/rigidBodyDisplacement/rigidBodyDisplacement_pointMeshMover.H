/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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
    Foam::pointMeshMovers::rigidBodyDisplacement

Description
    Rigid-body mesh motion solver for fvMesh.

SourceFiles
    rigidBodyDisplacement.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyDisplacement_pointMeshMover_H
#define rigidBodyDisplacement_pointMeshMover_H

#include "displacement_pointMeshMover.H"
#include "rigidBodyMotion.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
                  Class rigidBodyDisplacement Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyDisplacement
:
    public pointMeshMover,
    public RBD::rigidBodyMotion
{
    //- Class containing the patches and point motion weighting for each body
    class bodyMesh
    {
        //- Name of the body
        const word name_;

        //- ID of the body in the RBD::rigidBodyMotion
        const label bodyIndex_;

        //- List of mesh patches associated with this body
        const wordReList patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;


    public:

        friend class rigidBodyDisplacement;

        bodyMesh
        (
            const polyMesh& mesh,
            const word& name,
            const label bodyID,
            const dictionary& dict
        );
    };


    // Private Data

        //- List of the bodyMeshes containing the patches and point motion
        //  weighting for each body
        PtrList<bodyMesh> bodyMeshes_;

        //- Switch for test-mode in which only the
        //  gravitational body-force is applied
        Switch test_;

        //- Number of iterations in test-mode
        label nIter_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

        //- Ramp the forces according to the specified function and period
        autoPtr<Function1<scalar>> ramp_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        autoPtr<pointMeshMover> meshSolverPtr_;

        pointMeshMovers::displacement& meshSolver_;

        //- To avoid warning from clang
        using RBD::rigidBodyMotion::write;


public:

    //- Runtime type information
    TypeName("rigidBodyDisplacement");


    // Constructors

        //- Construct from polyMesh and dictionary
        rigidBodyDisplacement(const polyMesh&, const dictionary& dict);

        //- Disallow default bitwise copy construction
        rigidBodyDisplacement
        (
            const rigidBodyDisplacement&
        ) = delete;


    //- Destructor
    ~rigidBodyDisplacement();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> newPoints();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Write motion state information for restart
        virtual bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rigidBodyDisplacement&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
