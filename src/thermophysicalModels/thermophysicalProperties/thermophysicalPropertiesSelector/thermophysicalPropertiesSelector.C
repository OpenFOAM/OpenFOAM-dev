/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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

#include "thermophysicalPropertiesSelector.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermophysicalProperties>
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::
thermophysicalPropertiesSelector
(
    const word& name
)
:
    name_(name),
    propertiesPtr_(ThermophysicalProperties::New(name))
{}


template<class ThermophysicalProperties>
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::
thermophysicalPropertiesSelector
(
    const word& name,
    const dictionary& dict
)
:
    name_(name)
{
    const word type(dict.first()->keyword());

    if (dict.isDict(type))
    {
        propertiesPtr_ = ThermophysicalProperties::New(dict.subDict(type));
    }
    else
    {
        propertiesPtr_ = ThermophysicalProperties::New(type);
    }
}


template<class ThermophysicalProperties>
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::
thermophysicalPropertiesSelector
(
    const word& name,
    const thermophysicalPropertiesSelector& tps
)
:
    name_(tps.name),
    propertiesPtr_(tps.propertiesPtr_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermophysicalProperties>
void Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::write
(
    Ostream& os
) const
{
    propertiesPtr_->write(os);
}


// ************************************************************************* //
