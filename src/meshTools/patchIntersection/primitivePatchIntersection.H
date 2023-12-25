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

Class
    Foam::primitivePatchIntersection

Description
    Intersection between two primitive patches

SourceFiles
    primitivePatchIntersection.C

\*---------------------------------------------------------------------------*/

#ifndef primitivePatchIntersection_H
#define primitivePatchIntersection_H

#include "PatchIntersection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class primitivePatchIntersection Declaration
\*---------------------------------------------------------------------------*/

class primitivePatchIntersection
:
    public PatchIntersection<primitivePatch, primitivePatch>
{
public:

    // Runtime type information

        virtual word type() const
        {
            return "primitive" + patchIntersection::typeName.capitalise();
        }


    // Constructors

        //- Construct from a source and a target patch
        primitivePatchIntersection
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const scalar snapTol
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
