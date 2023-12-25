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
    Foam::FacePatchIntersection

Description
    Patch intersection based on polygonal faces. This triangulates the supplied
    patches and passes them to TriPatchIntersection. It then re-combines the
    triangles that TriPatchIntersection generates.

SourceFiles
    FacePatchIntersection.C

\*---------------------------------------------------------------------------*/

#ifndef FacePatchIntersection_H
#define FacePatchIntersection_H

#include "PatchIntersection.H"
#include "triFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class FacePatchIntersection Declaration
\*---------------------------------------------------------------------------*/

template<class SrcPatchType, class TgtPatchType>
class FacePatchIntersection
:
    public PatchIntersection<SrcPatchType, TgtPatchType>
{
public:

    // Runtime type information

        virtual word type() const
        {
            return "face" + patchIntersection::typeName.capitalise();
        }


    // Constructors

        //- Construct from a source and a target patch
        FacePatchIntersection
        (
            const SrcPatchType& srcPatch,
            const TgtPatchType& tgtPatch,
            const scalar snapTol
        );

        //- Construct from a source and a target patch, and specified source
        //  point normals
        FacePatchIntersection
        (
            const SrcPatchType& srcPatch,
            const vectorField& srcPointNormals,
            const TgtPatchType& tgtPatch,
            const scalar snapTol
        );


    //- Destructor
    virtual ~FacePatchIntersection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FacePatchIntersection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
