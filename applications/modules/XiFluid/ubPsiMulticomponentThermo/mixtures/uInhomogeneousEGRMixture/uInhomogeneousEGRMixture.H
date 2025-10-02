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
    Foam::uInhomogeneousEGRMixture

Description
    Unburnt gas ombustion mixture for inhomogeneous lean, stoichiometric
    or rich conditions with exhaust gas recirculation

    in which the fuel mass-fraction \c fu and exhaust gas recirculation
    mass-fraction \c egr define the composition.

SourceFiles
    uInhomogeneousEGRMixture.C

\*---------------------------------------------------------------------------*/

#ifndef uInhomogeneousEGRMixture_H
#define uInhomogeneousEGRMixture_H

#include "dimensionedTypes.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class uPsiMulticomponentThermo;

/*---------------------------------------------------------------------------*\
                     Class uInhomogeneousEGRMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class uInhomogeneousEGRMixture
{
public:

    // Public Typedefs

        //- The type of thermodynamics this mixture is instantiated for
        typedef ThermoType thermoType;

        //- Mixing type for thermodynamic properties
        typedef ThermoType thermoMixtureType;

        //- Mixing type for transport properties
        typedef ThermoType transportMixtureType;


private:

    // Private Enumerations

        //- Specie indices
        enum species { FU, EGR };


    // Private Data

        //- Stoichiometric air/fuel ratio
        scalar stoicRatio_;

        //- Fuel thermodynamic model
        thermoType fuel_;

        //- Oxidant thermodynamic model
        thermoType oxidant_;

        //- Product thermodynamic model
        thermoType products_;

        //- Mutable storage for the mixed thermodynamic model
        mutable thermoType mixture_;


    // Private Member Functions

        //- Return the mixture for the given composition
        const thermoType& mixture
        (
            const scalar fu,
            const scalar egr
        ) const;


public:

    // Constructors

        //- Construct from a dictionary
        uInhomogeneousEGRMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        uInhomogeneousEGRMixture
        (
            const uInhomogeneousEGRMixture<ThermoType>&
        ) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "uInhomogeneousEGRMixture<" + ThermoType::typeName() + '>';
        }

        //- Return the specie names
        static wordList specieNames()
        {
            return {"fu", "egr"};
        }

        const ThermoType& specieThermo(const label speciei) const;

        //- Return the fuel-oxidant equivalence ratio
        scalar Phi(const scalarFieldListSlice& Y) const;

        //- Return the mixture for thermodynamic properties
        const thermoMixtureType& thermoMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&,
            const thermoMixtureType&
        ) const;

        void reset
        (
            const volScalarField& b,
            PtrList<volScalarField>& Yu,
            const volScalarField& c,
            const PtrList<volScalarField>& Yb
        ) const;

        //- Read dictionary
        void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uInhomogeneousEGRMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "uInhomogeneousEGRMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
