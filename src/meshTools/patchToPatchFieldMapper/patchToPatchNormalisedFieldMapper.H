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
    Foam::patchToPatchNormalisedFieldMapper

Description
    Patch-to-patch fieldMapper which fills values for non-overlapping parts of
    the target patch by normalisation (for partly non-overlapping faces), or by
    referencing values from the nearest overlapping part (for completely
    non-overlapping faces).

\*---------------------------------------------------------------------------*/

#ifndef patchToPatchNormalisedFieldMapper_H
#define patchToPatchNormalisedFieldMapper_H

#include "patchToPatchFieldMapper.H"
#include "patchToPatchStabilisation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class patchToPatchNormalisedFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class patchToPatchNormalisedFieldMapper
:
    public patchToPatchFieldMapper
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
        patchToPatchNormalisedFieldMapper
        (
            const patchToPatch& pToP,
            const patchToPatchStabilisation& pS
        )
        :
            patchToPatchFieldMapper(pToP),
            pS_(pS)
        {}


    //- Destructor
    virtual ~patchToPatchNormalisedFieldMapper()
    {}


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_MAP_OPERATOR, );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
