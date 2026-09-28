/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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
    Foam::patchToPatchExtrapolation

Description
    Extrapolation procedure for patch-to-patch interpolations which copies the
    value from the nearest successfully mapped face.

SourceFiles
    patchToPatchExtrapolation.C

\*---------------------------------------------------------------------------*/

#ifndef nearestPatchToPatchExtrapolation_H
#define nearestPatchToPatchExtrapolation_H

#include "patchToPatchExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distributionMap;

namespace patchToPatchExtrapolations
{

/*---------------------------------------------------------------------------*\
                            Class nearest Declaration
\*---------------------------------------------------------------------------*/

class nearest
:
    public patchToPatchExtrapolation
{
    // Private Data

        //- For each uncoupled face, the local face used to extrapolate the
        //  interpolation
        labelList uncoupledFaceLocalFaces_;

        //- Map from faces to local extrapolation faces
        autoPtr<distributionMap> extrapolationMapPtr_;


    // Private Member Functions

        //- Extrapolate the given field if necessary
        template<class Type>
        void extrapolateType(Field<Type>& fld) const;


public:

    //- Run-time type information
    TypeName("nearest");


    // Constructors

        //- Construct null
        nearest();


    //- Destructor
    virtual ~nearest();


    // Member Functions

        //- Compute the extrapolation addressing if necessary
        void update
        (
            const polyPatch& patch,
            const PackedBoolList& faceCoupleds
        );

        //- Extrapolate the given field if necessary
        #define declareExtrapolateType(Type, nullArg) \
            virtual void extrapolate(Field<Type>& fld) const;
        FOR_ALL_FIELD_TYPES(declareExtrapolateType);
        #undef declareExtrapolateType
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatchExtrapolations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
