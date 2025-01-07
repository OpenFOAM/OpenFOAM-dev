/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Cloud with spherical, variable density, particles

SourceFiles
    particle.C

\*---------------------------------------------------------------------------*/

#ifndef particle_H
#define particle_H

#include "massive.H"
#include "coupledToFluid.H"
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
    public massive,
    public coupledToFluid,
    public sphericalCoupled
{
protected:

    // Protected Member Functions

        //- Initialisation hook
        virtual void initialise(const bool predict);

        //- Partition hook
        virtual void partition();

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


public:

    // Public Static Data

        //- Run-time type information
        TypeName("particle");


    // Constructors

        //- Construct from a mesh and a name
        particle
        (
            const polyMesh& mesh,
            const word& name,
            const contextType context,
            const IOobject::readOption readOption=IOobject::READ_IF_PRESENT,
            const IOobject::writeOption writeOption=IOobject::AUTO_WRITE
        );


    //- Destructor
    virtual ~particle();


    // Member Functions

        //- Solve the cloud's evolution over the current time-step
        virtual void solve();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
