/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::nonConformalMappedFvPatchBase

Description
    Base class for fv patches which provides non-conformal mapping between
    two potentially non-globally conforming fv patches

SourceFiles
    nonConformalMappedFvPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalMappedFvPatchBase_H
#define nonConformalMappedFvPatchBase_H

#include "mappedFvPatchBaseBase.H"
#include "nonConformalMappedPatchBase.H"
#include "nonConformalFvPatch.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class nonConformalMappedFvPatchBase Declaration
\*---------------------------------------------------------------------------*/

class nonConformalMappedFvPatchBase
:
    public mappedFvPatchBaseBase
{
    // Private Data

        //- Poly patch mapper base
        const nonConformalMappedPatchBase& mapper_;


    // Private Member Functions

        //- Map/interpolate a neighbour patch field to this patch
        template<class Type>
        tmp<Field<Type>> fromNeighbour(const Field<Type>& nbrFld) const;

        //- Map/interpolate a patch field to the neighbour patch
        template<class Type>
        tmp<Field<Type>> toNeighbour(const Field<Type>& fld) const;


public:

    //- Runtime type information
    TypeName("nonConformalMappedFvPatchBase");


    // Constructors

        //- Construct from a patch
        nonConformalMappedFvPatchBase(const fvPatch& patch);


    //- Destructor
    virtual ~nonConformalMappedFvPatchBase();


    // Member Functions

        //- Map/interpolate from one patch to another
        template<class Type>
        static tmp<Field<Type>> map
        (
            const fvsPatchLabelField& srcPolyFacesPf,
            const Field<Type>& srcFld,
            const fvsPatchLabelField& tgtPolyFacesPf
        );

        //- Map/interpolate from one patch to another
        template<class Type>
        static tmp<Field<Type>> map
        (
            const fvsPatchLabelField& srcPolyFacesPf,
            Field<Type>&& srcFld,
            const fvsPatchLabelField& tgtPolyFacesPf
        );

        //- Map/interpolate from one patch to another
        template<class Type>
        static tmp<Field<Type>> map
        (
            const fvsPatchLabelField& srcPolyFacesPf,
            const tmp<Field<Type>>& srcFld,
            const fvsPatchLabelField& tgtPolyFacesPf
        );

        //- Map/interpolate
        DEFINE_MAPPED_FV_PATCH_BASE_FROM_AND_TO_NEIGHBOUR(label, );

        //- Map/interpolate
        FOR_ALL_FIELD_TYPES
        (
            DEFINE_MAPPED_FV_PATCH_BASE_FROM_AND_TO_NEIGHBOUR,
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformalMappedFvPatchBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
