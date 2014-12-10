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

Class
    Foam::passiveParticle

Description
    Copy of base particle

SourceFiles
    passiveParticle.H

\*---------------------------------------------------------------------------*/

#ifndef passiveParticle_H
#define passiveParticle_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class passiveParticle Declaration
\*---------------------------------------------------------------------------*/

class passiveParticle
:
    public particle
{

public:

    // Constructors

        //- Construct from components
        passiveParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        )
        :
            particle(mesh, position, cellI, tetFaceI, tetPtI)
        {}

        //- Construct from components, with searching for tetFace and
        //  tetPt unless disabled by doCellFacePt = false.
        passiveParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            bool doCellFacePt = true
        )
        :
            particle(mesh, position, cellI, doCellFacePt)
        {}

        //- Construct from Istream
        passiveParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        )
        :
            particle(mesh, is, readFields)
        {}

        //- Construct as copy
        passiveParticle(const passiveParticle& p)
        :
            particle(p)
        {}

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new passiveParticle(*this));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
