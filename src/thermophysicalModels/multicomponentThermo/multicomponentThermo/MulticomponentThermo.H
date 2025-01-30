/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::MulticomponentThermo

Description
    Multi-component thermo implementation

SourceFiles
    MulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef MulticomponentThermo_H
#define MulticomponentThermo_H

#include "BasicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class MulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class MulticomponentThermo
:
    public BaseThermo
{
protected:

    // Protected Member Functions

        //- Return a volScalarField of the given property
        template<class Method, class ... Args>
        tmp<volScalarField> volScalarFieldPropertyi
        (
            const word& psiName,
            const dimensionSet& psiDim,
            Method psiMethod,
            const label speciei,
            const Args& ... args
        ) const;

        //- Return a volScalarField::Internal of the given property
        template<class Method, class ... Args>
        tmp<volScalarField::Internal> volInternalScalarFieldPropertyi
        (
            const word& psiName,
            const dimensionSet& psiDim,
            Method psiMethod,
            const label speciei,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property
        template<class Method, class Arg, class ... Args>
        tmp<scalarField> scalarFieldPropertyi
        (
            Method psiMethod,
            const label speciei,
            const Arg& arg,
            const Args& ... args
        ) const;


public:

    // Constructors

        //- Construct from mesh and phase name
        MulticomponentThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        MulticomponentThermo(const MulticomponentThermo<BaseThermo>&) = delete;


    //- Destructor
    virtual ~MulticomponentThermo();


    // Member Functions

        // Specie molecular properties

            //- Molecular weight [kg/kmol]
            virtual scalar WiValue(const label speciei) const;

            //- Molecular weight [kg/kmol]
            virtual dimensionedScalar Wi(const label speciei) const;


        // Specie thermodynamic properties

            //- Density [kg/m^3]
            virtual scalar rhoi
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rhoi
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual scalar Cpi
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cpi
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy/Internal energy [J/kg]
            virtual scalar hei
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> hei
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<volScalarField::Internal> hei
            (
                const label speciei,
                const volScalarField::Internal& p,
                const volScalarField::Internal& T
            ) const;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<volScalarField> hei
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual scalar hsi
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual tmp<scalarField> hsi
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField::Internal> hsi
            (
                const label speciei,
                const volScalarField::Internal& p,
                const volScalarField::Internal& T
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> hsi
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual scalar hai
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual tmp<scalarField> hai
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField::Internal> hai
            (
                const label speciei,
                const volScalarField::Internal& p,
                const volScalarField::Internal& T
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> hai
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy of formation [J/kg]
            virtual scalar hfiValue(const label speciei) const;

            //- Enthalpy of formation [J/kg]
            virtual dimensionedScalar hfi(const label speciei) const;


        // Specie transport properties

            //- Thermal conductivity [W/m/K]
            virtual scalar kappai
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappai
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "MulticomponentThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
