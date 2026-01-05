/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::clouds::particle

Description
    Cloud with spherical particles with thermodynamic modelling

SourceFiles
    particle.C

\*---------------------------------------------------------------------------*/

#ifndef particle_H
#define particle_H

#include "thermal.H"
#include "coupledToThermalFluid.H"
#include "massiveCoupledToFluid.H"
#include "sphericalCoupled.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                        Class particle Declaration
\*---------------------------------------------------------------------------*/

class particle
:
    public cloud,
    public spherical,
    public coupledToThermalFluid,
    public thermal,
    public sphericalCoupled,
    public massiveCoupledToFluid
{
protected:

    // Protected Member Functions

        //- Return the acceleration with which to do second-order tracking
        virtual tmp<LagrangianSubVectorField> dUdt
        (
            const LagrangianSubMesh&
        ) const;

        //- Do we need to re-calculate particles that are modified?
        virtual bool reCalculateModified();

        //- Update the cloud properties
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Partition hook
        virtual void partition();


public:

    // Public Static Data

        //- Run-time type information
        TypeName("particle");


    // Constructors

        //- Construct from a mesh and context
        particle
        (
            LagrangianMesh& mesh,
            const contextType context,
            const dictionary& dict
        );


    //- Destructor
    virtual ~particle();


    // Member Functions

        //- Solve the cloud's evolution over the current time-step
        virtual void solve(const bool initial, const bool final);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
