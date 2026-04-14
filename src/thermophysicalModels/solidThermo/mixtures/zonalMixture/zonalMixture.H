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
    Foam::zonalMixture

Description
    Foam::zonalMixture

SourceFiles
    zonalMixture.C

\*---------------------------------------------------------------------------*/

#ifndef zonalMixture_H
#define zonalMixture_H

#include "PtrList.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class zonalMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class zonalMixture
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

        //- Table of zone names
        hashedWordList zones_;

        //- List of zone thermos
        PtrList<ThermoType> zoneThermos_;


public:

    // Constructors

        //- Construct from a dictionary
        zonalMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        zonalMixture
        (
            const zonalMixture<ThermoType>&
        ) = delete;


    // Member Functions

        //- Return the instantiated type name
        inline static word typeName()
        {
            return "zonalMixture<" + ThermoType::typeName() + '>';
        }

        //- Return the table of zones
        inline const hashedWordList& zones() const
        {
            return zones_;
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return the raw zone thermodynamic data
        inline const PtrList<ThermoType>& zoneThermos() const
        {
            return zoneThermos_;
        }

        //- Return thermo based on index
        inline const ThermoType& zoneThermo(const label zonei) const
        {
            return zoneThermos_[zonei];
        }

        //- Return the mixture for thermodynamic properties
        inline const thermoMixtureType& thermoMixture
        (
            const label zonei
        ) const
        {
            return zoneThermos_[zonei];
        }

        //- Return the mixture for transport properties
        inline const transportMixtureType& transportMixture
        (
            const label zonei
        ) const
        {
            return zoneThermos_[zonei];
        }

        //- Return the mixture for transport properties
        inline const transportMixtureType& transportMixture
        (
            const label zonei,
            const thermoMixtureType&
        ) const
        {
            return zoneThermos_[zonei];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zonalMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

 // ************************************************************************* //
