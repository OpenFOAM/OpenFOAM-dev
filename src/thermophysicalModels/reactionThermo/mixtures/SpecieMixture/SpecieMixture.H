/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::SpecieMixture

Description
    Foam::SpecieMixture

SourceFiles
    SpecieMixture.C

\*---------------------------------------------------------------------------*/

#ifndef SpecieMixture_H
#define SpecieMixture_H

#include "scalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;
class fvMesh;


/*---------------------------------------------------------------------------*\
                        Class SpecieMixture Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class SpecieMixture
:
    public MixtureType
{

    // Private Member Functions

        //- Return a volScalarField of the given property
        tmp<volScalarField> volScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            scalar (MixtureType::thermoType::*psiMethod)
            (
                const scalar,
                const scalar
            ) const,
            const label speciei,
            const volScalarField& p,
            const volScalarField& T
        ) const;

        //- Return a field of the given property
        tmp<scalarField> fieldProperty
        (
            scalar (MixtureType::thermoType::*psiMethod)
            (
                const scalar,
                const scalar
            ) const,
            const label speciei,
            const scalarField& p,
            const scalarField& T
        ) const;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        SpecieMixture(const dictionary&, const fvMesh&, const word& phaseName);


    //- Destructor
    virtual ~SpecieMixture()
    {}


    // Member Functions

        // Specie properties

            //- Molecular weight of the given specie [kg/kmol]
            virtual scalar Wi(const label speciei) const;

            //- Enthalpy of formation [J/kg]
            virtual scalar Hf(const label speciei) const;

            //- Density [kg/m^3]
            virtual scalar rho
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual scalar Cp
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy/Internal energy [J/kg]
            virtual scalar HE
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> HE
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<volScalarField> HE
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual tmp<scalarField> Hs
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> Hs
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual tmp<scalarField> Ha
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> Ha
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;


        // Per specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Dynamic viscosity [kg/m/s]
            virtual tmp<volScalarField> mu
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Thermal conductivity [W/m/K]
            virtual scalar kappa
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappa
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;


        // Field properties abstract functions provided here to avoid issues
        // with inheritance of virtual functions in heThermo

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<volScalarField> Cp() const = 0;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const = 0;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Thermal diffusivity for temperature of mixture [W/m/K]
            virtual tmp<volScalarField> kappa() const = 0;

            //- Thermal diffusivity for temperature of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappa(const label patchi) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "SpecieMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
