/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::ThermoRefPair

Description
    Class containing a pair of thermo references. Handles down-casting to more
    specific thermo types by constructing one pair from another (e.g.,
    constructing a multicomponentThermo reference pair from a basicThermo
    pair). Tracks validity of the references.

\*---------------------------------------------------------------------------*/

#ifndef ThermoRefPair_H
#define ThermoRefPair_H

#include "physicalProperties.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ThermoRefPair Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class ThermoRefPair
{
    // Private Member Data

        //- Validity flags
        const Pair<bool> valid_;

        //- The first thermo
        const ThermoType& first_;

        //- The second thermo
        const ThermoType& second_;


public:

    // Constructors

        //- Construct from a database and phase names
        ThermoRefPair
        (
            const objectRegistry& db,
            const Pair<word>& phaseNames
        )
        :
            valid_(true, true),
            first_
            (
                db.lookupObject<ThermoType>
                (
                    IOobject::groupName
                    (
                        physicalProperties::typeName,
                        phaseNames.first()
                    )
                )
            ),
            second_
            (
                db.lookupObject<ThermoType>
                (
                    IOobject::groupName
                    (
                        physicalProperties::typeName,
                        phaseNames.second()
                    )
                )
            )
        {}

        //- Construct by casting a more primitive thermo type
        template<class BasicThermoType>
        ThermoRefPair(const ThermoRefPair<BasicThermoType> basicThermos)
        :
            valid_
            (
                isA<ThermoType>(basicThermos.first()),
                isA<ThermoType>(basicThermos.second())
            ),
            first_
            (
                valid_.first()
              ? refCast<const ThermoType>(basicThermos.first())
              : NullObjectRef<ThermoType>()
            ),
            second_
            (
                valid_.second()
              ? refCast<const ThermoType>(basicThermos.second())
              : NullObjectRef<ThermoType>()
            )
        {}


    // Member Functions

        //- Access the validity flags
        const Pair<bool>& valid() const
        {
            return valid_;
        }

        //- Access the first thermo
        const ThermoType& first() const
        {
            return first_;
        }

        //- Access the second thermo
        const ThermoType& second() const
        {
            return second_;
        }


    // Member Operators

        //- Access a thermo by index
        const ThermoType& operator[](const label i) const
        {
            return
                i == 0 ? first()
              : i == 1 ? second()
              : NullObjectRef<ThermoType>();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
