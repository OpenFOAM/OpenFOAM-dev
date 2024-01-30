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
    Foam::mappedFvPatchBaseBase

Description
    Base class for fv patches that provide mapping between two fv patches

SourceFiles
    mappedFvPatchBaseBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFvPatchBaseBase_H
#define mappedFvPatchBaseBase_H

#include "mappedPatchBaseBase.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DEFINE_MAPPED_FV_PATCH_BASE__NEIGHBOUR(Type, FromTo, Modifier)         \
                                                                               \
    virtual tmp<Field<Type>> FromTo##Neighbour                                 \
    (                                                                          \
        const Field<Type>&                                                     \
    ) const Modifier;                                                          \
                                                                               \
    virtual tmp<Field<Type>> FromTo##Neighbour                                 \
    (                                                                          \
        const tmp<Field<Type>>&                                                \
    ) const Modifier;


#define DEFINE_MAPPED_FV_PATCH_BASE_FROM_AND_TO_NEIGHBOUR(Type, Modifier)      \
                                                                               \
    DEFINE_MAPPED_FV_PATCH_BASE__NEIGHBOUR(Type, from, Modifier)               \
                                                                               \
    DEFINE_MAPPED_FV_PATCH_BASE__NEIGHBOUR(Type, to, Modifier)


#define IMPLEMENT_MAPPED_FV_PATCH_BASE__NEIGHBOUR(Type, PatchType, FromTo)     \
                                                                               \
    Foam::tmp<Foam::Field<Foam::Type>> Foam::PatchType::FromTo##Neighbour      \
    (                                                                          \
        const Field<Type>& nbrFldFld                                           \
    ) const                                                                    \
    {                                                                          \
        return FromTo##Neighbour<Type>(nbrFldFld);                             \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::Field<Foam::Type>> Foam::PatchType::FromTo##Neighbour      \
    (                                                                          \
        const tmp<Field<Type>>& nbrFldFld                                      \
    ) const                                                                    \
    {                                                                          \
        tmp<Field<Type>> tResult = FromTo##Neighbour<Type>(nbrFldFld());       \
        if (&tResult() == &nbrFldFld())                                        \
        {                                                                      \
            return nbrFldFld;                                                  \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            nbrFldFld.clear();                                                 \
            return tResult;                                                    \
        }                                                                      \
    }


#define IMPLEMENT_MAPPED_FV_PATCH_BASE_FROM_AND_TO_NEIGHBOUR(Type, PatchType)  \
                                                                               \
    IMPLEMENT_MAPPED_FV_PATCH_BASE__NEIGHBOUR(Type, PatchType, from)           \
                                                                               \
    IMPLEMENT_MAPPED_FV_PATCH_BASE__NEIGHBOUR(Type, PatchType, to)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                     Class mappedFvPatchBaseBase Declaration
\*---------------------------------------------------------------------------*/

class mappedFvPatchBaseBase
{
    // Private Data

        //- Patch to map to
        const fvPatch& patch_;

        //- Poly patch mapper base
        const mappedPatchBaseBase& mapper_;


public:

    //- Runtime type information
    TypeName("mappedFvPatchBaseBase");


    // Constructors

        //- Construct from a patch
        mappedFvPatchBaseBase(const fvPatch& patch);


    //- Destructor
    virtual ~mappedFvPatchBaseBase();


    // Member Functions

        //- Cast the given fvPatch to a mappedFvPatchBaseBase. Handle errors.
        static const mappedFvPatchBaseBase& getMap(const fvPatch& patch);


        // Access

            //- Reference to the fvPatch
            inline const fvPatch& patch() const;

            //- Name of the region to map from
            inline const word& nbrRegionName() const;

            //- Name of the patch to map from
            inline const word& nbrPatchName() const;

            //- The transformation between the patches
            inline const transformer& transform() const;

            //- Is the neighbour available?
            bool haveNbr() const;

            //- Get the mesh for the region to map from
            const fvMesh& nbrMesh() const;

            //- Get the patch to map from
            const fvPatch& nbrFvPatch() const;


        //- Map/interpolate
        FOR_ALL_FIELD_TYPES
        (
            DEFINE_MAPPED_FV_PATCH_BASE_FROM_AND_TO_NEIGHBOUR,
            = 0
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedFvPatchBaseBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
