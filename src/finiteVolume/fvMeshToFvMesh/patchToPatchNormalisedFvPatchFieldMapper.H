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
    Foam::patchToPatchNormalisedFvPatchFieldMapper

Description
    FieldMapper which uses a patchToPatch object to map from another patch. The
    source patch may be differently decomposed and/or geometrically and
    topologically different from the target.

\*---------------------------------------------------------------------------*/

#ifndef patchToPatchNormalisedFvPatchFieldMapper_H
#define patchToPatchNormalisedFvPatchFieldMapper_H

#include "patchToPatchFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class patchToPatchNormalisedFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class patchToPatchNormalisedFvPatchFieldMapper
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

        //- Construct given a patch-to-patch and a mapping direction
        patchToPatchNormalisedFvPatchFieldMapper
        (
            const patchToPatch& pToP,
            const bool forward
        )
        :
            patchToPatchFvPatchFieldMapper(pToP, forward)
        {}


    //- Destructor
    virtual ~patchToPatchNormalisedFvPatchFieldMapper()
    {}


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchToPatchNormalisedFvPatchFieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
