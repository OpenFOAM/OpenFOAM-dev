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

#include "tr.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(tr, 0);
    addToRunTimeSelectionTable(functionObject, tr, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<Foam::direction Rank>
template<class GeoField>
bool Foam::functionObjects::tr::CalcTr<Rank>::calc(tr& fo)
{
    if (!fo.foundObject<GeoField>(fo.fieldName_)) return false;

    Warning
        << "    functionObjects::" << fo.type() << " " << fo.name()
        << " required object " << fo.fieldName_ << " is a field of a "
        << "non-tensor type" << endl;

    return false;
}


template<class GeoField>
bool Foam::functionObjects::tr::CalcTr<Foam::direction(2)>::calc(tr& fo)
{
    if (!fo.foundObject<GeoField>(fo.fieldName_)) return false;

    const GeoField& geoField =
        fo.lookupObject<GeoField>(fo.fieldName_);

    fo.store(fo.resultName_, Foam::tr(geoField));

    return true;
}


bool Foam::functionObjects::tr::calc()
{
    bool processed = false;

    #define calcTrType(Type, GeoField)                                         \
        processed =                                                            \
            processed                                                          \
         || CalcTr<pTraits<Type>::rank>::calc<GeoField<Type>>(*this);
    FOR_ALL_FIELD_TYPES(calcTrType, VolField);
    FOR_ALL_FIELD_TYPES(calcTrType, VolInternalField);
    FOR_ALL_FIELD_TYPES(calcTrType, SurfaceField);
    #undef calcTrType

    if (!processed)
    {
        cannotFindObject(fieldName_);
    }

    return processed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::tr::tr
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::tr::~tr()
{}


// ************************************************************************* //
