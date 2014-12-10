/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

public:

    // Constructors

        //- Construct from dictionary and mesh
        SpecieMixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~SpecieMixture()
    {}


    // Member functions

        // Per specie properties

            //- Number of moles []
            virtual scalar nMoles(const label speciei) const;

            //- Molecular weight [kg/kmol]
            virtual scalar W(const label speciei) const;


        // Per specie thermo properties

            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar Cp
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar Cv
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Chemical enthalpy [J/kg]
            virtual scalar Hc(const label speciei) const;

            //- Entropy [J/(kg K)]
            virtual scalar S
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Sensible internal energy [J/kg]
            virtual scalar Es
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Gibbs free energy [J/kg]
            virtual scalar G
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Helmholtz free energy [J/kg]
            virtual scalar A
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;


        // Per specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Thermal conductivity [W/m/K]
            virtual scalar kappa
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Thermal diffusivity of enthalpy [kg/m/s]
            virtual scalar alphah
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Density [kg/m3]
            virtual scalar rho
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "SpecieMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
