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

\*---------------------------------------------------------------------------*/

#include "caseFileConfiguration.H"
#include "Time.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::caseFileConfiguration::beginDict(OFstream& os)
{
    os << indent << "{" << incrIndent << endl;
}


void Foam::caseFileConfiguration::beginDict(OFstream& os,const word& name)
{
    os << indent << name << endl;
    beginDict(os);
}


void Foam::caseFileConfiguration::endDict(OFstream& os, bool newline)
{
    os << decrIndent << indent << "}" << endl;

    if (newline)
    {
        os << endl;
    }
}


void Foam::caseFileConfiguration::beginList(OFstream& os)
{
    os  << indent << "(" << incrIndent << endl;
}


void Foam::caseFileConfiguration::beginList(OFstream& os, const word& name)
{
    os  << indent << name << endl;
    beginList(os);
}


void Foam::caseFileConfiguration::endList(OFstream& os, bool newline)
{
    os  << decrIndent << indent << ");" << endl;

    if (newline)
    {
        os  << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::caseFileConfiguration::caseFileConfiguration
(
    const fileName& name,
    const fileName& dir,
    const Time& time
)
:
    dict_(IOobject(name, dir, time)),
    os_(dict_.objectPath(true))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::caseFileConfiguration::~caseFileConfiguration()
{}


// ************************************************************************* //
