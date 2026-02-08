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
    Foam::UInhomogeneousMixture

Description
    Templated unburnt gas combustion mixture for inhomogeneous lean,
    stoichiometric or rich conditions

    in which the fuel mass-fraction \c fu defines the composition.

See also
    Foam::uInhomogeneousMixture

SourceFiles
    UInhomogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef UInhomogeneousMixture_H
#define UInhomogeneousMixture_H

#include "uInhomogeneousMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class UInhomogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class UInhomogeneousMixture
:
    public uInhomogeneousMixture
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

    // Private Data

        //- Fuel thermodynamic model
        thermoType fuel_;

        //- Oxidant thermodynamic model
        thermoType oxidant_;

        //- List of specie active flags
        boolList active_;

        //- Mutable storage for the mixed thermodynamic model
        mutable thermoType mixture_;


    // Private Member Functions

        //- Return the mixture for the given composition
        const thermoType& mixture(const scalar fu) const;


public:

    // Constructors

        //- Construct from a dictionary
        UInhomogeneousMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        UInhomogeneousMixture
        (
            const UInhomogeneousMixture<ThermoType>&
        ) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                uInhomogeneousMixture::typeName
              + '<' + ThermoType::typeName() + '>';
        }

        inline const ThermoType& specieThermo(const label speciei) const
        {
            return fuel_;
        }

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
        void operator=(const UInhomogeneousMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UInhomogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
