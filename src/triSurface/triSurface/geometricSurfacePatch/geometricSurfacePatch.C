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

#include "geometricSurfacePatch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(geometricSurfacePatch, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
geometricSurfacePatch::geometricSurfacePatch()
:
    geometricType_("empty"),
    name_("patch"),
    index_(0)
{}


// Construct from components
geometricSurfacePatch::geometricSurfacePatch
(
    const word& geometricType,
    const word& name,
    const label index
)
:
    geometricType_(geometricType),
    name_(name),
    index_(index)

{
    if (geometricType_.empty())
    {
        geometricType_ = "empty";
    }
}


// Construct from Istream
geometricSurfacePatch::geometricSurfacePatch(Istream& is, const label index)
:
    geometricType_(is),
    name_(is),
    index_(index)
{
    if (geometricType_.empty())
    {
        geometricType_ = "empty";
    }
}


// Construct from dictionary
geometricSurfacePatch::geometricSurfacePatch
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    geometricType_(dict.lookup("geometricType")),
    name_(name),
    index_(index)
{
    if (geometricType_.empty())
    {
        geometricType_ = "empty";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Write
void geometricSurfacePatch::write(Ostream& os) const
{
    os  << nl << name_
        << nl << geometricType_;
}


void geometricSurfacePatch::writeDict(Ostream& os) const
{
    os  << "    geometricType " << geometricType_ << ';' << nl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::geometricSurfacePatch::operator!=(const geometricSurfacePatch& p)
    const
{
    return !(*this == p);
}


bool Foam::geometricSurfacePatch::operator==(const geometricSurfacePatch& p)
    const
{
    return
    (
        (geometricType() == p.geometricType())
     && (name() == p.name())
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Istream& operator>>(Istream& is, geometricSurfacePatch& gp)
{
    is >> gp.name_ >> gp.geometricType_;

    return is;
}


Ostream& operator<<(Ostream& os, const geometricSurfacePatch& gp)
{
    gp.write(os);
    os.check
    (
        "Ostream& operator<<(Ostream& f, const geometricSurfacePatch& gp)"
    );
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
