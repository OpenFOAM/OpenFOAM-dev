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
    Foam::rigidBodyMeshMotion

Description
    Rigid-body mesh motion solver for fvMesh.

    Applies septernion interpolation of movement as function of distance to the
    object surface.

SourceFiles
    rigidBodyMeshMotion.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyMeshMotion_H
#define rigidBodyMeshMotion_H

#include "multiRigidBodyMeshMotion.H"
#include "rigidBodyMotion.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class rigidBodyMeshMotion Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyMeshMotion
:
    public multiRigidBodyMeshMotion,
    public RBD::rigidBodyMotion
{
    // Private Data

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


    // Private Member Functions

        //- To avoid warning from clang
        using RBD::rigidBodyMotion::write;

        //- Return the list of transformation for each body from time 0
        virtual List<septernion> transforms0() const;

        //- Move the bodies
        virtual void moveBodies();


public:

    //- Runtime type information
    TypeName("rigidBodyMotion");


    // Constructors

        //- Construct from polyMesh and dictionary
        rigidBodyMeshMotion
        (
            const word& name,
            const polyMesh&,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rigidBodyMeshMotion
        (
            const rigidBodyMeshMotion&
        ) = delete;


    //- Destructor
    ~rigidBodyMeshMotion();


    // Member Functions

        //- Write motion state information for restart
        virtual bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rigidBodyMeshMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
