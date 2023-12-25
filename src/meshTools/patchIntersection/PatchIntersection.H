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
    Foam::PatchIntersection

Description
    Base class for patch intersections. Provides storage and access to the
    intersection points and faces and their relationship between the source and
    target patches.

SourceFiles
    PatchIntersectionI.H
    PatchIntersection.C

\*---------------------------------------------------------------------------*/

#ifndef PatchIntersection_H
#define PatchIntersection_H

#include "patchIntersection.H"
#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PatchIntersection Declaration
\*---------------------------------------------------------------------------*/

template<class SrcPatchType, class TgtPatchType>
class PatchIntersection
:
    public patchIntersection
{
protected:

    // Protected Data

        //- Reference to the source patch
        const SrcPatchType& srcPatch_;

        //- Reference to the target patch
        const TgtPatchType& tgtPatch_;


public:

    // Constructors

        //- Construct from a source and a target patch
        PatchIntersection
        (
            const SrcPatchType& srcPatch,
            const TgtPatchType& tgtPatch
        );


    //- Destructor
    virtual ~PatchIntersection();


    // Member Functions

        // Access

            //- The source patch
            inline const SrcPatchType& srcPatch() const;

            //- The target patch
            inline const TgtPatchType& tgtPatch() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PatchIntersectionI.H"

#ifdef NoRepository
    #include "PatchIntersection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
