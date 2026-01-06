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
    Foam::fluidLagrangianThermo

Description
    Base-class for fluid Lagrangian thermodynamic models

SourceFiles
    fluidLagrangianThermo.C
    fluidLagrangianThermoTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fluidLagrangianThermo_H
#define fluidLagrangianThermo_H

#include "FluidLagrangianThermo.H"
#include "pureLagrangianThermo.H"
#include "uniformDimensionedFields.H"
#include "uniformGeometricFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fluidLagrangianThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidLagrangianThermo
:
    virtual public basicLagrangianThermo
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
            FluidLagrangianThermo
            <
                BasicLagrangianThermo
                <
                    MixtureType,
                    composite
                >
            >;

        //- The derived name
        static word derivedThermoName()
        {
            return "heRhoThermo";
        }


    //- Runtime type information
    TypeName("fluidLagrangianThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidLagrangianThermo,
        LagrangianMesh,
        (const LagrangianMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Select a fluid thermo
        static autoPtr<fluidLagrangianThermo> New
        (
            const LagrangianMesh& mesh,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~fluidLagrangianThermo();


    // Member Functions

        // Thermodynamic State

            //- Pressure [Pa]
            virtual const LagrangianScalarDynamicField& p() const = 0;

            //- Modify the Pressure [Pa]
            virtual LagrangianScalarDynamicField& p() = 0;

            //- Pressure for a sub-mesh [Pa]
            virtual tmp<LagrangianSubScalarSubField> p
            (
                const LagrangianSubMesh&
            ) const = 0;

            //- Pressure for an injection [Pa]
            virtual tmp<LagrangianSubScalarField> p
            (
                const LagrangianInjection&,
                const LagrangianSubMesh&
            ) const = 0;


        // Thermodynamic Properties

            //- Compressibility [s^2/m^2]
            virtual const LagrangianScalarDynamicField& psi() const = 0;

            //- Compressibility for an injection [m^2/s^2]
            virtual tmp<LagrangianSubScalarField> psi
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const = 0;


        // Transport Properties

            //- Dynamic viscosity [kg/m/s]
            virtual const LagrangianScalarDynamicField& mu() const = 0;

            //- Dynamic viscosity for an injection [W/m/K]
            virtual tmp<LagrangianSubScalarField> mu
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const = 0;
};


/*---------------------------------------------------------------------------*\
             Class fluidLagrangianThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class fluidLagrangianThermo::implementation
:
    virtual public fluidLagrangianThermo
{
protected:

    // Protected Member Data

        //- Pressure [Pa]
        LagrangianScalarDynamicField p_;

        //- Pressure source condition [Pa]
        autoPtr<LagrangianScalarFieldSource> pSourcePtr_;

        //- Compressibility [s^2/m^2]
        LagrangianScalarDynamicField psi_;

        //- Dynamic Viscosity [kg/m/s]
        LagrangianScalarDynamicField mu_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        implementation
        (
            const dictionary&,
            const LagrangianMesh&,
            const word&,
            const LagrangianScalarDynamicField& T
        );


    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- Initialise state
        virtual void initialise();

        //- Update the pressure
        virtual void correctPressure(const LagrangianSubMesh&);


        // Thermodynamic State

            //- Pressure [Pa]
            virtual const LagrangianScalarDynamicField& p() const;

            //- Modify the Pressure [Pa]
            virtual LagrangianScalarDynamicField& p();

            //- Pressure for a sub-mesh [Pa]
            virtual tmp<LagrangianSubScalarSubField> p
            (
                const LagrangianSubMesh&
            ) const;

            //- Pressure for an injection [Pa]
            virtual tmp<LagrangianSubScalarField> p
            (
                const LagrangianInjection&,
                const LagrangianSubMesh&
            ) const;


        // Thermodynamic Properties

            //- Compressibility [s^2/m^2]
            virtual const LagrangianScalarDynamicField& psi() const;


        // Transport Properties

            //- Dynamic viscosity [kg/m/s]
            virtual const LagrangianScalarDynamicField& mu() const;
};


/*---------------------------------------------------------------------------*\
             Class fluidLagrangianThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class fluidLagrangianThermo::composite
:
    public basicLagrangianThermo::implementation,
    public pureLagrangianThermo,
    public fluidLagrangianThermo::implementation
{
public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        composite
        (
            const dictionary& dict,
            const LagrangianMesh& mesh,
            const word& phaseName
        )
        :
            basicLagrangianThermo::implementation(dict, mesh, phaseName),
            fluidLagrangianThermo::implementation(dict, mesh, phaseName, T_)
        {}

        //- Construct from dictionary, mesh and phase name
        template<class MixtureType>
        composite
        (
            const dictionary& dict,
            const MixtureType& mixture,
            const LagrangianMesh& mesh,
            const word& phaseName
        )
        :
            composite(dict, mesh, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
