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
    Foam::patchToPatchLeftOverFieldMapper

Description
    Patch-to-patch fieldMapper which retains values in the target field for
    parts of the patch that do not overlap the source. This process needs an
    input target field as well as a source field, so it can only map a field
    in-place.

\*---------------------------------------------------------------------------*/

#ifndef patchToPatchLeftOverFieldMapper_H
#define patchToPatchLeftOverFieldMapper_H

#include "patchToPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class patchToPatchLeftOverFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class patchToPatchLeftOverFieldMapper
:
    public patchToPatchFieldMapper
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
        using patchToPatchFieldMapper::patchToPatchFieldMapper;


    //- Destructor
    virtual ~patchToPatchLeftOverFieldMapper()
    {}


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
