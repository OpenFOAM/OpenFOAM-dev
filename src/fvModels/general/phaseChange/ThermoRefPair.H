/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "fvModel.H"
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


    // Private Constructors

        //- Construct from components
        inline ThermoRefPair
        (
            Pair<bool> valid,
            const ThermoType& first,
            const ThermoType& second
        );


public:

    //- Declare friendship with other thermo reference pairs
    template<class OtherThermoType>
    friend class ThermoRefPair;


    // Constructors

        //- Construct from a database and phase names
        inline ThermoRefPair
        (
            const objectRegistry& db,
            const Pair<word>& phaseNames
        );


    // Member Functions

        //- Access the validity flags
        inline const Pair<bool>& valid() const;

        //- Return if either validity flag is set
        inline bool either() const;

        //- Return if both validity flags are set
        inline bool both() const;

        //- Access the first thermo
        inline const ThermoType& first() const;

        //- Access the second thermo
        inline const ThermoType& second() const;

        //- Cast to a different thermo type, with error checking
        template<class OtherThermoType>
        inline ThermoRefPair<OtherThermoType> thermos() const;

        //- Cast to a different thermo type, with error checking
        template<class OtherThermoType>
        inline ThermoRefPair<OtherThermoType> thermos
        (
            const Pair<bool> required,
            const fvModel& model,
            const word& thermoType
        ) const;


    // Member Operators

        //- Access a thermo by index
        inline const ThermoType& operator[](const label i) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoRefPairI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
