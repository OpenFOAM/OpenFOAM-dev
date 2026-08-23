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
    Base-class for thermodynamic models of the carrier phase of a Lagrangian
    cloud. Stores interpolations of the properties stored by the Eulerian
    thermodynamic model. Also constructs low-level modelling in the same way as
    the Eulerian model does, so carrier properties can be evaluated at a
    different thermodynamic state.

SourceFiles
    carrierThermo.C

\*---------------------------------------------------------------------------*/

#ifndef carrierThermo_H
#define carrierThermo_H

#include "carrierDensity.H"
#include "basicLagrangianSubThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;

/*---------------------------------------------------------------------------*\
                   Class carrierThermo Declaration
\*---------------------------------------------------------------------------*/

class carrierThermo
:
    virtual public carrierDensity,
    virtual public basicLagrangianSubThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    // Public Typedefs

        //- The derived name
        static word derivedThermoName()
        {
            return "carrierThermo";
        }


    //- Runtime type information
    TypeName("carrierThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        carrierThermo,
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
        static autoPtr<carrierThermo> New
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const word& phaseName
        );

        //- Select a carrier thermo
        template<class CarrierThermo, class Thermo>
        static autoPtr<CarrierThermo> New
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const word& phaseName
        );


    //- Destructor
    virtual ~carrierThermo();


    // Member Functions

        //- Eulerian thermo
        virtual const basicThermo& thermo() const = 0;

        //- Phase name
        virtual const word& phaseName() const = 0;

        //- Correct
        virtual void correct() = 0;


        // Thermodynamic State

            using basicLagrangianSubThermo::T;

            //- Temperature [K]
            virtual const CarrierField<scalar>& T() const = 0;


        // Thermodynamic Properties

            using basicLagrangianSubThermo::rho;

            //- Density [kg/m^3]
            virtual const CarrierField<scalar>& rho() const = 0;

            //- Enthalpy/energy [J/kg]
            virtual const CarrierField<scalar>& he() const = 0;

            //- Heat capacity [J/kg/K]
            virtual const CarrierField<scalar>& Cpv() const = 0;

            using basicLagrangianSubThermo::Cp;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual const CarrierField<scalar>& Cp() const = 0;

            using basicLagrangianSubThermo::Cv;

            //- Heat capacity at constant volume [J/kg/K]
            virtual const CarrierField<scalar>& Cv() const = 0;


        // Transport Properties

            using basicLagrangianSubThermo::kappa;

            //- Thermal conductivity [W/m/K]
            virtual const CarrierField<scalar>& kappa() const = 0;
};


/*---------------------------------------------------------------------------*\
             Class carrierThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class carrierThermo::implementation
:
    virtual public carrierThermo
{
    // Private Data

        //- Eulerian thermo
        const basicThermo& thermo_;

        //- Lagrangian mesh
        const LagrangianMesh& mesh_;

        //- Carrier temperature
        const CarrierField<scalar>& T_;

        //- Cached carrier density vol-field
        tmp<volScalarField> trhoVf_;

        //- Carrier density
        const CarrierField<scalar>& rho_;

        //- Carrier enthalpy/energy
        const CarrierField<scalar>& he_;

        //- Carrier heat capacity
        const CarrierField<scalar>& Cpv_;

        //- Carrier heat capacity at constant pressure
        const CarrierField<scalar>& Cp_;

        //- Carrier heat capacity at constant volume
        const CarrierField<scalar>& Cv_;

        //- Carrier thermal conductivity
        const CarrierField<scalar>& kappa_;


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

        //- Eulerian thermo
        virtual const basicThermo& thermo() const;

        //- Access the properties dictionary
        virtual const IOdictionary& properties() const;

        //- Access the mesh
        virtual const LagrangianMesh& mesh() const;

        //- Phase name
        virtual const word& phaseName() const;

        //- Correct
        virtual void correct();


        // Thermodynamic State

            //- Temperature [K]
            virtual const CarrierField<scalar>& T() const;

            //- Temperature for a sub-mesh [K]
            virtual tmp<LagrangianSubScalarSubField> T
            (
                const LagrangianSubMesh& subMesh
            ) const;


        // Thermodynamic Properties

            //- Density [kg/m^3]
            virtual const CarrierField<scalar>& rho() const;

            //- Density for a sub-mesh [kg/m^3]
            virtual tmp<LagrangianSubScalarSubField> rho
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Enthalpy/energy [J/kg]
            virtual const CarrierField<scalar>& he() const;

            //- Heat capacity [J/kg/K]
            virtual const CarrierField<scalar>& Cpv() const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual const CarrierField<scalar>& Cp() const;

            //- Heat capacity at constant pressure for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarField> Cp
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual const CarrierField<scalar>& Cv() const;

            //- Heat capacity at constant volume for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarSubField> Cv
            (
                const LagrangianSubMesh& subMesh
            ) const;


        // Transport Properties

            //- Thermal conductivity [W/m/K]
            virtual const CarrierField<scalar>& kappa() const;

            //- Thermal conductivity for a sub-mesh [W/m/K]
            virtual tmp<LagrangianSubScalarSubField> kappa
            (
                const LagrangianSubMesh& subMesh
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "carrierThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
