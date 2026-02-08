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
    Foam::UInhomogeneousEGRMixture

Description
    Templated Unburnt gas combustion mixture for inhomogeneous lean,
    stoichiometric or rich conditions with exhaust gas recirculation

    in which the fuel mass-fraction \c fu and exhaust gas recirculation
    mass-fraction \c egr define the composition.

See also
    Foam::uInhomogeneousEGRMixture

SourceFiles
    UInhomogeneousEGRMixture.C

\*---------------------------------------------------------------------------*/

#ifndef UInhomogeneousEGRMixture_H
#define UInhomogeneousEGRMixture_H

#include "uInhomogeneousEGRMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class UInhomogeneousEGRMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class UInhomogeneousEGRMixture
:
    public uInhomogeneousEGRMixture
{
public:

    // Public Typedefs

        //- The type of thermodynamics this mixture is instantiated for
        typedef ThermoType thermoType;

        //- Mixing type for thermodynamic properties
        typedef ThermoType thermoMixtureType;

        //- Mixing type for transport properties
        typedef ThermoType transportMixtureType;


protected:

    // Protected Data

        //- Fuel thermodynamic model
        thermoType fuel_;

        //- Oxidant thermodynamic model
        thermoType oxidant_;

        //- Product thermodynamic model
        thermoType products_;

        //- Mutable storage for the mixed thermodynamic model
        mutable thermoType mixture_;


    // Protected Member Functions

        //- Return the mixture for the given composition
        const thermoType& mixture
        (
            const scalar fu,
            const scalar egr
        ) const;


public:

    // Constructors

        //- Construct from a dictionary
        UInhomogeneousEGRMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        UInhomogeneousEGRMixture
        (
            const UInhomogeneousEGRMixture<ThermoType>&
        ) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                uInhomogeneousEGRMixture::typeName
              + '<' + ThermoType::typeName() + '>';
        }

        const ThermoType& specieThermo(const label speciei) const;

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

        //- Read dictionary
        void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const UInhomogeneousEGRMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UInhomogeneousEGRMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
