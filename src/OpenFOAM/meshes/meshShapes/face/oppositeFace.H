/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::oppositeFace

Description
    Class containing opposite face for a prismatic cell with addressing
    and a possibility of failure.

SourceFiles
    oppositeFace.C

\*---------------------------------------------------------------------------*/

#ifndef oppositeFace_H
#define oppositeFace_H

#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class oppositeFace Declaration
\*---------------------------------------------------------------------------*/

class oppositeFace
:
    public face
{
    // Private Data

        //- Master face index
        const label masterIndex_;

        //- Opposite face index
        const label oppositeIndex_;


public:

    // Constructors

        //- Construct from components
        oppositeFace
        (
            const face& f,
            const label masterIndex,
            const label oppositeIndex
        )
        :
            face(f),
            masterIndex_(masterIndex),
            oppositeIndex_(oppositeIndex)

        {}


    // Member Functions

        //- Master face index
        inline label masterIndex() const
        {
            return masterIndex_;
        }

        //- Slave face index
        inline label oppositeIndex() const
        {
            return oppositeIndex_;
        }

        //- Does the opposite face exist?
        inline bool found() const
        {
            return oppositeIndex_ >= 0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
