/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "boundaryPatch.H"
#include "dictionary.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryPatch::boundaryPatch
(
    const word& name,
    const label index,
    const label size,
    const label start,
    const word& physicalType
)
:
    patchIdentifier(name, index, physicalType),
    size_(size),
    start_(start)
{}


Foam::boundaryPatch::boundaryPatch
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    patchIdentifier(name, dict, index),
    size_(readLabel(dict.lookup("nFaces"))),
    start_(readLabel(dict.lookup("startFace")))
{}


Foam::boundaryPatch::boundaryPatch(const boundaryPatch& p, const label index)
:
    patchIdentifier(p.name(), index, p.physicalType()),
    size_(p.size()),
    start_(p.start())
{}


Foam::autoPtr<Foam::boundaryPatch> Foam::boundaryPatch::clone() const
{
    return autoPtr<boundaryPatch>(new boundaryPatch(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryPatch::~boundaryPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryPatch::write(Ostream& os) const
{
    patchIdentifier::write(os);
    writeEntry(os, "nFaces", size_);
    writeEntry(os, "startFace", start_);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundaryPatch& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const boundaryPatch&)");
    return os;
}


// ************************************************************************* //
