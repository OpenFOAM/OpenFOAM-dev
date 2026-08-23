/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::carrierThermo

Description
    Class which provides thermodynamic properties of the carrier phase

SourceFiles
    fluidCarrierThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidCarrierThermo_H
#define fluidCarrierThermo_H

#include "carrierThermo.H"
#include "carrierDynamicViscosity.H"
#include "fluidLagrangianSubThermo.H"
#include "pureCarrierThermo.H"
#include "CarrierThermo.H"
#include "FluidCarrierThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidThermo;

/*---------------------------------------------------------------------------*\
                   Class fluidCarrierThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidCarrierThermo
:
    virtual public carrierThermo,
    virtual public carrierDynamicViscosity,
    virtual public fluidLagrangianSubThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            FluidCarrierThermo
            <
                CarrierThermo
                <
                    MixtureType,
                    composite
                >
            >;


    //- Runtime type information
    TypeName("fluidCarrierThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidCarrierThermo,
        cloud,
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const basicThermo& thermo
        ),
        (c, carriedCloud, thermo)
    );


    // Selectors

        //- Select a carrier thermo
        static autoPtr<fluidCarrierThermo> New
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const word& phaseName
        );


    //- Destructor
    virtual ~fluidCarrierThermo();


    // Member Functions

        // Thermodynamic State

            using fluidLagrangianSubThermo::p;

            //- Pressure [Pa]
            virtual const CarrierField<scalar>& p() const = 0;


        // Transport Properties

            using fluidLagrangianSubThermo::mu;

            //- Dynamic viscosity [kg/m/s]
            virtual const CloudDerivedField<scalar>& mu() const = 0;

            //- Kinematic viscosity [m^2/s]
            virtual const CloudDerivedField<scalar>& nu() const = 0;
};


/*---------------------------------------------------------------------------*\
             Class fluidCarrierThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class fluidCarrierThermo::implementation
:
    virtual public fluidCarrierThermo
{
    // Private Data

        //- Eulerian thermo
        const fluidThermo& thermo_;

        //- Carrier pressure
        const CarrierField<scalar>& p_;

        //- Carrier dynamic viscosity
        const CarrierField<scalar>& mu_;

        //- Carrier kinematic viscosity
        const CloudDerivedField<scalar>& nu_;


public:

    // Constructors

        //- Construct from a reference to the cloud and its base classes
        implementation
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const basicThermo& thermo
        );

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Thermodynamic State

            //- Pressure [Pa]
            virtual const CarrierField<scalar>& p() const;

            //- Pressure for a sub-mesh [Pa]
            virtual tmp<LagrangianSubScalarSubField> p
            (
                const LagrangianSubMesh&
            ) const;


        // Transport Properties

            //- Dynamic viscosity [kg/m/s]
            virtual const CloudDerivedField<scalar>& mu() const;

            //- Dynamic viscosity for a sub-mesh [kg/m/s]
            virtual tmp<LagrangianSubScalarSubField> mu
            (
                const LagrangianSubMesh&
            ) const;

            //- Kinematic viscosity [m^2/s]
            virtual const CloudDerivedField<scalar>& nu() const;

            //- Kinematic viscosity for a sub-mesh [kg/m/s]
            virtual tmp<LagrangianSubScalarField> nu
            (
                const LagrangianSubMesh&
            ) const;
};


/*---------------------------------------------------------------------------*\
                 Class fluidCarrierThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class fluidCarrierThermo::composite
:
    public carrierThermo::implementation,
    public fluidCarrierThermo::implementation,
    public pureCarrierThermo
{
public:

    // Constructors

        //- Construct from a reference to the cloud and its base classes
        composite
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const basicThermo& thermo
        )
        :
            carrierThermo::implementation(c, carriedCloud, thermo),
            fluidCarrierThermo::implementation(c, carriedCloud, thermo)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
