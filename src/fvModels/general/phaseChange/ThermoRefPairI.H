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

#include "ThermoRefPair.H"
#include "physicalProperties.H"

// * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::ThermoRefPair<ThermoType>::ThermoRefPair
(
    Pair<bool> valid,
    const ThermoType& first,
    const ThermoType& second
)
:
    valid_(valid),
    first_(first),
    second_(second)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::ThermoRefPair<ThermoType>::ThermoRefPair
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


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
inline const Foam::Pair<bool>& Foam::ThermoRefPair<ThermoType>::valid() const
{
    return valid_;
}


template<class ThermoType>
inline bool Foam::ThermoRefPair<ThermoType>::either() const
{
    return valid_.first() || valid_.second();
}


template<class ThermoType>
inline bool Foam::ThermoRefPair<ThermoType>::both() const
{
    return valid_.first() && valid_.second();
}


template<class ThermoType>
inline const ThermoType& Foam::ThermoRefPair<ThermoType>::first() const
{
    return first_;
}


template<class ThermoType>
inline const ThermoType& Foam::ThermoRefPair<ThermoType>::second() const
{
    return second_;
}


template<class ThermoType>
template<class OtherThermoType>
inline Foam::ThermoRefPair<OtherThermoType>
Foam::ThermoRefPair<ThermoType>::thermos() const
{
    Pair<bool> valid
    (
        this->valid().first() && isA<OtherThermoType>(first()),
        this->valid().second() && isA<OtherThermoType>(second())
    );

    return
        ThermoRefPair<OtherThermoType>
        (
            valid,
            valid.first()
          ? refCast<const OtherThermoType>(first())
          : NullObjectRef<OtherThermoType>(),
            valid.second()
          ? refCast<const OtherThermoType>(second())
          : NullObjectRef<OtherThermoType>()
        );
}


template<class ThermoType>
template<class OtherThermoType>
inline Foam::ThermoRefPair<OtherThermoType>
Foam::ThermoRefPair<ThermoType>::thermos
(
    const Pair<bool> required,
    const fvModel& model,
    const word& thermoType
) const
{
    ThermoRefPair<OtherThermoType> result = thermos<OtherThermoType>();

    forAll(result.valid(), i)
    {
        if (!result.valid()[i] && required[i])
        {
            FatalErrorInFunction
                << "Model " << model.name() << " of type "
                << model.type() << " requires a " << thermoType
                << " thermo for phase " << operator[](i).phaseName()
                << exit(FatalError);
        }
    }

    return result;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class ThermoType>
inline const ThermoType&
Foam::ThermoRefPair<ThermoType>::operator[](const label i) const
{
    return
        i == 0 ? first()
      : i == 1 ? second()
      : NullObjectRef<ThermoType>();
}


// ************************************************************************* //
