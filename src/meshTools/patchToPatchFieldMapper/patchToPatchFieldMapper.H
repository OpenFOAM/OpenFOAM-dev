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
    Foam::patchToPatchFieldMapper

Description
    Base class for fieldMappers which uses a patchToPatch object to map from
    another patch. The source patch may be differently decomposed and/or
    geometrically and topologically different from the target. Derived into
    versions which handle unmapped regions differently.

\*---------------------------------------------------------------------------*/

#ifndef patchToPatchFieldMapper_H
#define patchToPatchFieldMapper_H

#include "fieldMapper.H"
#include "patchToPatch.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class patchToPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class patchToPatchFieldMapper
:
    public fieldMapper
{
protected:

    // Protected Data

        //- Patch-to-patch interpolation engine
        const patchToPatch& pToP_;


public:

    // Constructors

        //- Construct given a patch-to-patch interpolation
        patchToPatchFieldMapper(const patchToPatch& pToP)
        :
            fieldMapper(),
            pToP_(pToP)
        {}


    //- Destructor
    virtual ~patchToPatchFieldMapper()
    {}


    // Member Operators

        //- Map a label field (not implemented)
        DEFINE_FIELD_MAPPER_MAP_OPERATOR(label, );

        //- Inherit the base call operator
        using fieldMapper::operator();

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
    #include "patchToPatchFieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
