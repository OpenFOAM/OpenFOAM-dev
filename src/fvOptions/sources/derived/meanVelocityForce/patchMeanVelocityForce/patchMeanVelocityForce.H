/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::fv::patchMeanVelocityForce

Description
    Calculates and applies the force necessary to maintain the specified mean
    velocity averaged over the specified patch.

    Note: Currently only handles kinematic pressure (incompressible solvers).

Usage
    Example usage:
    \verbatim
    selectionMode   all;                    // Apply force to all cells
    fields          (U);                    // Name of velocity field
    patch           inlet;                  // Name of the patch
    Ubar            (10.0 0 0);             // Desired mean velocity
    relaxation      0.2;                    // Optional relaxation factor
    \endverbatim

SourceFiles
    patchMeanVelocityForce.C

\*---------------------------------------------------------------------------*/

#ifndef patchMeanVelocityForce_H
#define patchMeanVelocityForce_H

#include "meanVelocityForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class patchMeanVelocityForce Declaration
\*---------------------------------------------------------------------------*/

class patchMeanVelocityForce
:
    public meanVelocityForce
{
protected:

    // Protected data

        //- Patch name
        word patch_;

        //- Patch index
        label patchi_;


    // Protected Member Functions

        //- Calculate and return the magnitude of the mean velocity
        //  averaged over the specified patch
        virtual scalar magUbarAve(const volVectorField& U) const;


public:

    //- Runtime type information
    TypeName("patchMeanVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        patchMeanVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        patchMeanVelocityForce(const patchMeanVelocityForce&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchMeanVelocityForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
