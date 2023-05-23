/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::patchToPatchFvPatchFieldMapper

Description
    FieldMapper which uses a patchToPatch object to map from another patch. The
    source patch may be differently decomposed and/or geometrically and
    topologically different from the target.

\*---------------------------------------------------------------------------*/

#ifndef patchToPatchFvPatchFieldMapper_H
#define patchToPatchFvPatchFieldMapper_H

#include "fvPatchFieldMapper.H"
#include "patchToPatch.H"
#include "patchToPatchStabilisation.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class patchToPatchFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class patchToPatchFvPatchFieldMapper
:
    public fvPatchFieldMapper
{
protected:

    // Protected Data

        //- Patch-to-patch interpolation engine
        const patchToPatch& pToP_;


public:

    // Constructors

        //- Construct given a patch-to-patch interpolation
        patchToPatchFvPatchFieldMapper(const patchToPatch& pToP)
        :
            fvPatchFieldMapper(),
            pToP_(pToP)
        {}


    //- Destructor
    virtual ~patchToPatchFvPatchFieldMapper()
    {}


    // Member Operators

        //- Map a label field (not implemented)
        DEFINE_FIELD_MAPPER_OPERATOR(label, );

        using fvPatchFieldMapper::operator();

        //- Map a temporary field
        template<class Type>
        void operator()(Field<Type>& f, const tmp<Field<Type>>& tmapF) const;

        //- Map a temporary field
        template<class Type>
        tmp<Field<Type>> operator()(const tmp<Field<Type>>& tmapF) const;
};


/*---------------------------------------------------------------------------*\
            Class patchToPatchLeftOverFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class patchToPatchLeftOverFvPatchFieldMapper
:
    public patchToPatchFvPatchFieldMapper
{
    // Private Member Functions

        //- Map from one field to another
        template<class Type>
        void map(Field<Type>& f, const Field<Type>& mapF) const;

        //- Map a field and return the result as tmp
        template<class Type>
        tmp<Field<Type>> map(const Field<Type>& mapF) const;


public:

    // Constructors

        //- Inherit base class constructor
        using patchToPatchFvPatchFieldMapper::patchToPatchFvPatchFieldMapper;


    //- Destructor
    virtual ~patchToPatchLeftOverFvPatchFieldMapper()
    {}


    // Member Functions

        //- Return whether or not all faces receive a mapped value
        virtual bool hasUnmapped() const
        {
            FatalErrorInFunction
                << "Not a valid query for this mapper, which should only be "
                << "used for modifying an existing, valid, field"
                << exit(FatalError);
            return false;
        }


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, );
};


/*---------------------------------------------------------------------------*\
           Class patchToPatchNormalisedFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class patchToPatchNormalisedFvPatchFieldMapper
:
    public patchToPatchFvPatchFieldMapper
{
    // Private Data

        //- Patch stabilisation engine
        const patchToPatchStabilisation& pS_;


    // Private Member Functions

        //- Map from one field to another
        template<class Type>
        void map(Field<Type>& f, const Field<Type>& mapF) const;

        //- Map a field and return the result as tmp
        template<class Type>
        tmp<Field<Type>> map(const Field<Type>& mapF) const;


public:

    // Constructors

        //- Construct given a patch-to-patch interpolation and stabilisation
        patchToPatchNormalisedFvPatchFieldMapper
        (
            const patchToPatch& pToP,
            const patchToPatchStabilisation& pS
        )
        :
            patchToPatchFvPatchFieldMapper(pToP),
            pS_(pS)
        {}


    //- Destructor
    virtual ~patchToPatchNormalisedFvPatchFieldMapper()
    {}


    // Member Functions

        //- Return whether or not all faces receive a mapped value
        virtual bool hasUnmapped() const
        {
            return false;
        }


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, );
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchToPatchFvPatchFieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
