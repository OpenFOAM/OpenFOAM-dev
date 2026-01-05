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
    Foam::BasicLagrangianThermo

Description
    Lagrangian thermodynamic model implementation and storage of energy.
    Provides overloads of the functions defined in the basic Lagrangian
    thermodynamic type that depend on the primitive thermodynamic model.

SourceFiles
    BasicLagrangianThermo.C

\*---------------------------------------------------------------------------*/

#ifndef BasicLagrangianThermo_H
#define BasicLagrangianThermo_H

#include "basicLagrangianThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BasicThermoName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(BasicLagrangianThermo);


/*---------------------------------------------------------------------------*\
                          Class basicLagrangianThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType, class BasicThermoType>
class BasicLagrangianThermo
:
    public BasicLagrangianThermoName,
    public physicalProperties,
    public MixtureType,
    public BasicThermoType
{
protected:

    // Protected Data

        //- Internal energy [J/kg]
        LagrangianScalarDynamicField e_;


    // Protected Member Functions

        //- Return a LagrangianInternalScalarField of the given property
        template<class Mixture, class Method, class ... Args>
        tmp<LagrangianInternalScalarField> LagrangianInternalScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            Mixture mixture,
            Method psiMethod,
            const Args& ... args
        ) const;

        //- Return a LagrangianSubScalarField of the given property
        template<class Mixture, class Method, class ... Args>
        tmp<LagrangianSubScalarField> LagrangianSubScalarFieldProperty
        (
            const LagrangianSubMesh& subMesh,
            const word& psiName,
            const dimensionSet& psiDim,
            Mixture mixture,
            Method psiMethod,
            const Args& ... args
        ) const;

        //- Return a LagrangianSubScalarField of the given property
        template<class Mixture, class Method, class ... Args>
        tmp<LagrangianSubScalarField> LagrangianInjectionProperty
        (
            const LagrangianInjection& injection,
            const LagrangianSubMesh& subMesh,
            const word& psiName,
            const dimensionSet& psiDim,
            Mixture mixture,
            Method psiMethod,
            const Args& ... args
        ) const;


public:

    // Typedefs

        //- Mixture type
        typedef MixtureType mixtureType;

        //- Basic thermo
        typedef BasicThermoType basicThermoType;


    //- Disambiguate debug switch used by derivations
    using BasicLagrangianThermoName::debug;


    // Constructors

        //- Construct from mesh and phase name
        BasicLagrangianThermo(const LagrangianMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        BasicLagrangianThermo
        (
            const BasicLagrangianThermo<MixtureType, BasicThermoType>&
        ) = delete;


    //- Destructor
    virtual ~BasicLagrangianThermo();


    // Member Functions

        //- Access the properties dictionary
        virtual const IOdictionary& properties() const;

        //- Modify the properties dictionary
        virtual IOdictionary& properties();

        //- Return the full name of the thermodynamic model
        virtual word thermoName() const;


        // Molecular properties

            //- Molecular weight [kg/kmol]
            virtual tmp<LagrangianSubScalarField> W
            (
                const LagrangianSubMesh&
            ) const;


        // Thermodynamic State

            //- Internal energy [J/kg]
            virtual const LagrangianScalarDynamicField& e() const;

            //- Modify the internal energy [J/kg]
            virtual LagrangianScalarDynamicField& e();


        // Derived Thermodynamic Properties

            //- Density for an injection [kg/m^3]
            virtual tmp<LagrangianSubScalarField> rho
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const;

            //- Internal energy for an injection [J/kg]
            virtual tmp<LagrangianSubScalarField> e
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const;

            //- Heat capacity at constant volume for an injection [J/kg/K]
            virtual tmp<LagrangianSubScalarField> Cv
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const;


        // Derived Transport Properties

            //- Thermal conductivity for an injection [W/m/K]
            virtual tmp<LagrangianSubScalarField> kappa
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const;


        //- Read thermophysical properties dictionary
        using BasicThermoType::read;

        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BasicLagrangianThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
