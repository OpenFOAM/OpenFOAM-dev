/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "IOobject.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::Ostream& Foam::operator<<(Ostream& os, const InfoProxy<IOobject>& ip)
{
    const IOobject& io = ip.t_;

    os  << "IOobject: "
        << io.type() << token::SPACE
        << io.name() << token::SPACE
        << "readOpt:" << token::SPACE << io.readOpt() << token::SPACE
        << "writeOpt:" << token::SPACE << io.writeOpt() << token::SPACE
        << "globalObject:" << token::SPACE << io.globalObject() << token::SPACE
        << io.path() << endl;

    return os;
}


// ************************************************************************* //
