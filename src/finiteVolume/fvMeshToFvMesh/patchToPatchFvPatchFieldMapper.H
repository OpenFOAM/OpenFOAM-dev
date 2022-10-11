/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    // Private Data

        //- Patch-to-patch mapping engine
        const patchToPatch& pToP_;

        //- Mapping direction
        const bool forward_;

        //- Size of the mapped field
        label size_;

        //- Does the mapping engine leave anything unmapped?
        bool hasUnmapped_;


    // Private Member Functions

        //- Map from one field to another
        template<class Type>
        void map(Field<Type>& f, const Field<Type>& mapF) const;

        //- Map a field and return the result as tmp
        template<class Type>
        tmp<Field<Type>> map(const Field<Type>& mapF) const;


public:

    // Constructors

        //- Construct given a patch-to-patch and a mapping direction
        patchToPatchFvPatchFieldMapper
        (
            const patchToPatch& pToP,
            const bool forward
        )
        :
            fvPatchFieldMapper(),
            pToP_(pToP),
            forward_(forward),
            size_(-1),
            hasUnmapped_(false)
        {
            const List<List<remote>> procFaces =
                forward_ ? pToP_.tgtSrcProcFaces() : pToP_.srcTgtProcFaces();

            size_ = procFaces.size();

            forAll(procFaces, i)
            {
                if (procFaces[i].size() == 0)
                {
                    hasUnmapped_ = true;
                }
            }
        }


    //- Destructor
    virtual ~patchToPatchFvPatchFieldMapper()
    {}


    // Member Functions

        //- Return whether or not all faces receive a mapped value
        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, );

        //- Map a label field
        DEFINE_FIELD_MAPPER_OPERATOR(label, );

        //- Map a temporary field
        template<class Type>
        void operator()(Field<Type>& f, const tmp<Field<Type>>& tmapF) const;

        //- Map a temporary field
        template<class Type>
        tmp<Field<Type>> operator()(const tmp<Field<Type>>& tmapF) const;
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
