/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "thermophysicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermophysicalProperties, 0);
    defineRunTimeSelectionTable(thermophysicalProperties,);
    defineRunTimeSelectionTable(thermophysicalProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermophysicalProperties::thermophysicalProperties(scalar W)
:
    W_(W)
{}


Foam::thermophysicalProperties::thermophysicalProperties(const dictionary& dict)
:
    W_(readScalar(dict.lookup("W")))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermophysicalProperties>
Foam::thermophysicalProperties::New
(
    const word& name
)
{
    if (debug)
    {
        InfoInFunction << "Constructing thermophysicalProperties" << endl;
    }

    ConstructorTable::iterator cstrIter = ConstructorTablePtr_->find(name);

    if (cstrIter == ConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown thermophysicalProperties type "
            << name << nl << nl
            << "Valid thermophysicalProperties types are:" << nl
            << ConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<thermophysicalProperties>(cstrIter()());
}


Foam::autoPtr<Foam::thermophysicalProperties>
Foam::thermophysicalProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing thermophysicalProperties" << endl;
    }

    const word& thermophysicalPropertiesTypeName = dict.dictName();

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(thermophysicalPropertiesTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown thermophysicalProperties type "
            << thermophysicalPropertiesTypeName << nl << nl
            << "Valid thermophysicalProperties types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<thermophysicalProperties>(cstrIter()(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermophysicalProperties::readIfPresent(const dictionary &dict)
{
    dict.readIfPresent("W", W_);
}


void Foam::thermophysicalProperties::writeData(Ostream& os) const
{
    os  << W_;
}


void Foam::thermophysicalProperties::write(Ostream& os) const
{
    dictionary dict("thermophysicalProperties");
    dict.add("W", W_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const thermophysicalProperties& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
