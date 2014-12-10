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

#include "particle.H"
#include "IOstreams.H"
#include "IOPosition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::particle::propertyList_ = Foam::particle::propertyList();


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
    // readFields : read additional data. Should be consistent with writeFields.

    if (is.format() == IOstream::ASCII)
    {
        is  >> position_ >> cellI_;

        if (readFields)
        {
            is  >> tetFaceI_ >> tetPtI_ >> origProc_ >> origId_;
        }
    }
    else
    {
        // In binary read all particle data - needed for parallel transfer
        if (readFields)
        {
            is.read
            (
                reinterpret_cast<char*>(&position_),
                sizeof(position_)
              + sizeof(cellI_)
              + sizeof(faceI_)
              + sizeof(stepFraction_)
              + sizeof(tetFaceI_)
              + sizeof(tetPtI_)
              + sizeof(origProc_)
              + sizeof(origId_)
            );
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&position_),
                sizeof(position_)
              + sizeof(cellI_)
              + sizeof(faceI_)
              + sizeof(stepFraction_)
            );
        }
    }

    // Check state of Istream
    is.check("particle::particle(Istream&, bool)");
}


void Foam::particle::write(Ostream& os, bool writeFields) const
{
    if (os.format() == IOstream::ASCII)
    {
        if (writeFields)
        {
            // Write the additional entries
            os  << position_
                << token::SPACE << cellI_
                << token::SPACE << tetFaceI_
                << token::SPACE << tetPtI_
                << token::SPACE << origProc_
                << token::SPACE << origId_;
        }
        else
        {
            os  << position_
                << token::SPACE << cellI_;
        }
    }
    else
    {
        // In binary write both cellI_ and faceI_, needed for parallel transfer
        if (writeFields)
        {
            os.write
            (
                reinterpret_cast<const char*>(&position_),
                sizeof(position_)
              + sizeof(cellI_)
              + sizeof(faceI_)
              + sizeof(stepFraction_)
              + sizeof(tetFaceI_)
              + sizeof(tetPtI_)
              + sizeof(origProc_)
              + sizeof(origId_)
            );
        }
        else
        {
            os.write
            (
                reinterpret_cast<const char*>(&position_),
                sizeof(position_)
              + sizeof(cellI_)
              + sizeof(faceI_)
              + sizeof(stepFraction_)
            );
        }
    }

    // Check state of Ostream
    os.check("particle::write(Ostream& os, bool) const");
}


Foam::Ostream& Foam::operator<<(Ostream& os, const particle& p)
{
    // Write all data
    p.write(os, true);

    return os;
}


// ************************************************************************* //
