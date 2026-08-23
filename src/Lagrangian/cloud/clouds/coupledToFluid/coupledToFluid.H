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
    Foam::clouds::coupledToFluid

Description
    Base class for clouds which are coupled to a variable density fluid

SourceFiles
    coupledToFluid.C

\*---------------------------------------------------------------------------*/

#ifndef coupledToFluid_H
#define coupledToFluid_H

#include "coupled.H"
#include "carrierDynamicViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                           Class coupledToFluid Declaration
\*---------------------------------------------------------------------------*/

class coupledToFluid
:
    private autoPtr<carrierDynamicViscosity>,
    public coupled
{
private:

    // Private Data

        //- Reference to the carrier viscosity
        const carrierDynamicViscosity& muc_;

        //- Optional reference to the fluid surface thermo
        const fluidSurfaceThermo& thermocsOrNull_;

        //- Corresponding Eulerian density
        autoPtr<carrierDensity> rhocPhasePtr_;


protected:

    // Protected Member Functions

        //- Update the cached carrier fields
        void updateCarrier();


public:

    // Public Static Data

        //- Run-time type information
        TypeName("coupledToFluid");


    // Public Data

        //- Carrier density
        const CarrierField<scalar>& rhoc;

        //- Corresponding Eulerian density
        const CarrierField<scalar>& rhocPhase;


    // Constructors

        //- Construct from a reference to the cloud
        coupledToFluid
        (
            const cloud& c,
            const carried&,
            const carrierDynamicViscosity& mucOrNull =
                NullObjectRef<carrierDynamicViscosity>(),
            const fluidSurfaceThermo& thermocsOrNull =
                NullObjectRef<fluidSurfaceThermo>(),
            const carrierDensity& rhocPhaseOrNull =
                NullObjectRef<carrierDensity>()
        );


    //- Destructor
    virtual ~coupledToFluid();


    // Member Functions

        //- Return the carrier dynamic viscosity
        tmp<LagrangianSubScalarSubField> muc(const LagrangianSubMesh&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
