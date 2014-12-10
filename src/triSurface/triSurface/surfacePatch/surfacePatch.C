/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "surfacePatch.H"
#include "dictionary.H"
#include "word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(surfacePatch, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
Foam::surfacePatch::surfacePatch()
:
    geometricSurfacePatch("", "", -1),
    size_(0),
    start_(0)
{}



// Construct from components
Foam::surfacePatch::surfacePatch
(
    const word& geometricType,
    const word& name,
    const label size,
    const label start,
    const label index
)
:
    geometricSurfacePatch(geometricType, name, index),
    size_(size),
    start_(start)
{}


// Construct from Istream
Foam::surfacePatch::surfacePatch(Istream& is, const label index)
:
    geometricSurfacePatch(is, index),
    size_(0),
    start_(0)
{
    size_ = readLabel(is);
    start_ = readLabel(is);
}

// Construct from dictionary
Foam::surfacePatch::surfacePatch
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    geometricSurfacePatch(name, dict, index),
    size_(readLabel(dict.lookup("nFaces"))),
    start_(readLabel(dict.lookup("startFace")))
{}


// Construct as copy
Foam::surfacePatch::surfacePatch(const Foam::surfacePatch& sp)
:
    geometricSurfacePatch(sp),
    size_(sp.size()),
    start_(sp.start())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfacePatch::write(Ostream& os) const
{
    os  << nl
        << static_cast<const geometricSurfacePatch&>(*this) << endl
        << size() << tab << start();
}

void Foam::surfacePatch::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl;

    geometricSurfacePatch::writeDict(os);

    os  << "    nFaces " << size() << ';' << nl
        << "    startFace " << start() << ';' << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::surfacePatch::operator!=(const surfacePatch& p) const
{
    return !(*this == p);
}


bool Foam::surfacePatch::operator==(const surfacePatch& p) const
{
    return
    (
        (geometricType() == p.geometricType())
     && (size() == p.size())
     && (start() == p.start())
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const surfacePatch& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const surfacePatch& p");
    return os;
}


// ************************************************************************* //
