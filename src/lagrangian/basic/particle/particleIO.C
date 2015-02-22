/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "particle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::particle::propertyList_ = Foam::particle::propertyList();

const std::size_t Foam::particle::sizeofPosition_
(
    offsetof(particle, faceI_) - offsetof(particle, position_)
);

const std::size_t Foam::particle::sizeofFields_
(
    sizeof(particle) - offsetof(particle, position_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particle::particle(const polyMesh& mesh, Istream& is, bool readFields)
:
    mesh_(mesh),
    position_(),
    cellI_(-1),
    faceI_(-1),
    stepFraction_(0.0),
    tetFaceI_(-1),
    tetPtI_(-1),
    origProc_(Pstream::myProcNo()),
    origId_(-1)
{
    if (is.format() == IOstream::ASCII)
    {
        is  >> position_ >> cellI_;

        if (readFields)
        {
            is  >> faceI_
                >> stepFraction_
                >> tetFaceI_
                >> tetPtI_
                >> origProc_
                >> origId_;
        }
    }
    else
    {
        if (readFields)
        {
            is.read(reinterpret_cast<char*>(&position_), sizeofFields_);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&position_), sizeofPosition_);
        }
    }

    // Check state of Istream
    is.check("particle::particle(Istream&, bool)");
}


void Foam::particle::writePosition(Ostream& os) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << position_ << token::SPACE << cellI_;
    }
    else
    {
        os.write(reinterpret_cast<const char*>(&position_), sizeofPosition_);
    }

    // Check state of Ostream
    os.check("particle::writePosition(Ostream& os, bool) const");
}


Foam::Ostream& Foam::operator<<(Ostream& os, const particle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << p.position_
            << token::SPACE << p.cellI_
            << token::SPACE << p.faceI_
            << token::SPACE << p.stepFraction_
            << token::SPACE << p.tetFaceI_
            << token::SPACE << p.tetPtI_
            << token::SPACE << p.origProc_
            << token::SPACE << p.origId_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&p.position_),
            particle::sizeofFields_
        );
    }

    return os;
}


// ************************************************************************* //
