/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "Function1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::unitConversions::unitConversions
(
    std::initializer_list<unitConversion> l
)
:
    x(dimless),
    value(dimless)
{
    auto i = l.begin();
    x.reset(*i);
    value.reset(*(++i));
}


Foam::Function1s::unitConversions::unitConversions(Istream& is)
:
    x(dimless),
    value(dimless)
{
    is >> *this;
}


Foam::autoPtr<Foam::Function1s::unitConversions>
Foam::Function1s::unitConversions::clone() const
{
    return autoPtr<unitConversions>(new unitConversions(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::Function1s::unitConversions::readIfPresent
(
    const word& keyword,
    const dictionary& dict
)
{
    const entry* entryPtr = dict.lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        ITstream& is = entryPtr->stream();

        is.readBegin("unitConversions");
        x.read(keyword, dict, is);
        value.read(keyword, dict, is);
        is.readEnd("unitConversions");

        return true;
    }
    else
    {
        if (dictionary::writeOptionalEntries)
        {
            IOInfoInFunction(dict)
                << "Optional entry '" << keyword << "' is not present,"
                << " the default value '" << token::BEGIN_LIST << x.info()
                << token::SPACE << value.info() << token::END_LIST
                << "' will be used." << endl;
        }

        return false;
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::assertNoConvertUnits
(
    const word& typeName,
    const Function1s::unitConversions& units,
    const dictionary& dict
)
{
    if (!units.x.standard() || !units.value.standard())
    {
        FatalIOErrorInFunction(dict)
            << "Unit conversions are not supported by "
            << typeName << " function1 types" << abort(FatalError);
    }
}


void Foam::writeEntry(Ostream& os, const Function1s::unitConversions& units)
{
    os << units;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    Function1s::unitConversions& units
)
{
    is.readBegin("unitConversions");
    is >> units.x >> units.value;
    is.readEnd("unitConversions");

    is.check("Istream& operator>>(Istream&, unitConversions&)");

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Function1s::unitConversions& units
)
{
    os  << token::BEGIN_LIST
        << units.x << token::SPACE << units.value
        << token::END_LIST;

    os.check("Ostream& operator<<(Ostream&, const unitConversions&)");

    return os;
}


// ************************************************************************* //
