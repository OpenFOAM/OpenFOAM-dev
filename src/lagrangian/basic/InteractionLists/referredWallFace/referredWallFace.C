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

#include "referredWallFace.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::referredWallFace::referredWallFace()
:
    face(),
    pts_(),
    patchI_()
{}


Foam::referredWallFace::referredWallFace
(
    const face& f,
    const pointField& pts,
    label patchI
)
:
    face(f),
    pts_(pts),
    patchI_(patchI)
{
    if (this->size() != pts_.size())
    {
        FatalErrorIn
        (
            "Foam::referredWallFace::referredWallFace"
            "("
                "const face& f, "
                "const pointField& pts, "
                "label patchI"
            ")"
        )   << "Face and pointField are not the same size. " << nl << (*this)
            << abort(FatalError);
    }
}


Foam::referredWallFace::referredWallFace(const referredWallFace& rWF)
:
    face(rWF),
    pts_(rWF.pts_),
    patchI_(rWF.patchI_)
{
    if (this->size() != pts_.size())
    {
        FatalErrorIn
        (
            "Foam::referredWallFace::referredWallFace"
            "("
                "const referredWallFace& rWF"
            ")"
        )   << "Face and pointField are not the same size. " << nl << (*this)
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::referredWallFace::~referredWallFace()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

bool Foam::referredWallFace::operator==(const referredWallFace& rhs) const
{
    return
    (
        static_cast<const face&>(rhs) == static_cast<face>(*this)
     && rhs.pts_ == pts_
     && rhs.patchI_ == patchI_
    );
}


bool Foam::referredWallFace::operator!=(const referredWallFace& rhs) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, referredWallFace& rWF)
{
    is  >> static_cast<face&>(rWF) >> rWF.pts_ >> rWF.patchI_;

    // Check state of Istream
    is.check
    (
        "Foam::Istream& "
        "Foam::operator>>(Foam::Istream&, Foam::referredWallFace&)"
    );

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const referredWallFace& rWF)
{
    os  << static_cast<const face&>(rWF) << token::SPACE
        << rWF.pts_ << token::SPACE
        << rWF.patchI_;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::referredWallFace&)"
    );

    return os;
}


// ************************************************************************* //
