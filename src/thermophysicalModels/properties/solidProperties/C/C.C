/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "C.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C, 0);
    addToRunTimeSelectionTable(solidProperties, C,);
    addToRunTimeSelectionTable(solidProperties, C, Istream);
    addToRunTimeSelectionTable(solidProperties, C, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C::C()
:
    solidProperties(2010, 710, 0.04, 0.0, 1.0)
{
    if (debug)
    {
        WarningIn("C::C()")
            << "Properties of graphite need to be checked!!!"
            << endl;
    }
}


Foam::C::C(const solidProperties& s)
:
    solidProperties(s)
{}


Foam::C::C(Istream& is)
:
    solidProperties(is)
{}


Foam::C::C(const dictionary& dict)
:
    solidProperties(dict)
{}


Foam::C::C(const C& s)
:
    solidProperties(s)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C::writeData(Ostream& os) const
{
    solidProperties::writeData(os);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const C& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
