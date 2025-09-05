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

\*---------------------------------------------------------------------------*/

#include "coupled.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(coupled, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#define ACCESS_CARRIER_EQNS(Type, nullArg)                                     \
namespace Foam                                                                 \
{                                                                              \
    namespace clouds                                                           \
    {                                                                          \
        template<>                                                             \
        PtrDictionary<CarrierEqn<Type>>& coupled::carrierEqns() const          \
        {                                                                      \
            return CAT3(carrier, CAPITALIZE(Type), Eqns_);                     \
        }                                                                      \
    }                                                                          \
}
FOR_ALL_FIELD_TYPES(ACCESS_CARRIER_EQNS)
#undef ACCESS_CARRIER_EQNS


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::clouds::coupled::clearCarrierEqns()
{
    #define CLEAR_TYPE_CARRIER_EQNS(Type, nullArg)                             \
        forAllIter                                                             \
        (                                                                      \
            PtrDictionary<CarrierEqn<Type>>,                                   \
            carrierEqns<Type>(),                                               \
            iter                                                               \
        )                                                                      \
        {                                                                      \
            iter().clear();                                                    \
        }
    FOR_ALL_FIELD_TYPES(CLEAR_TYPE_CARRIER_EQNS);
    #undef CLEAR_TYPE_CARRIER_EQNS
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupled::coupled(const cloud& c)
:
    carried(c),
    nuc(c.derivedField<scalar>(*this, &coupled::calcNuc))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupled::~coupled()
{}


// ************************************************************************* //
