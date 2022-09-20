/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        )
        :
            particle(mesh, coordinates, celli, tetFacei, tetPti)
        {}


        //- Construct from a position and a cell, searching for the rest of the
        //  required topology
        passiveParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        )
        :
            particle(mesh, position, celli)
        {}

        //- Construct from Istream
        passiveParticle(Istream& is, bool readFields = true)
        :
            particle(is, readFields)
        {}

        //- Copy constructor
        passiveParticle(const passiveParticle& p)
        :
            particle(p)
        {}

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new passiveParticle(*this));
        }

        //- Construct from Istream and return
        static autoPtr<passiveParticle> New(Istream& is)
        {
            return autoPtr<passiveParticle>(new passiveParticle(is));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
