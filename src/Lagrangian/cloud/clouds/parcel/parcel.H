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
    Foam::clouds::parcel

Description
    Basic cloud with spherical, variable density, particles, grouped into
    parcels.

    This cloud's modelling is the same as for particle, except that an
    additional number field is used to represent multiple physical particles
    per Lagrangian parcel. This adds flexibility with regards to managing the
    expense of simulating large numbers of physical particles, and also
    provides the ability to specify the number, size (distribution) and volume
    or mass flow rate of particles at an injection.

See also
    Foam::clouds::particle

SourceFiles
    parcel.C

\*---------------------------------------------------------------------------*/

#ifndef parcel_H
#define parcel_H

#include "grouped.H"
#include "massive.H"
#include "coupledToFluid.H"
#include "sphericalCoupled.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                           Class parcel Declaration
\*---------------------------------------------------------------------------*/

class parcel
:
    public cloud,
    public grouped,
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
        TypeName("parcel");


    // Constructors

        //- Construct from a mesh and a name
        parcel
        (
            const polyMesh& mesh,
            const word& name,
            const contextType context,
            const IOobject::readOption readOption=IOobject::READ_IF_PRESENT,
            const IOobject::writeOption writeOption=IOobject::AUTO_WRITE
        );


    //- Destructor
    virtual ~parcel();


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
