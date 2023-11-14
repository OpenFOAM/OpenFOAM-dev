/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::setSizeFieldMapper

Description
    Mapper which sets the field size. It does not actually map values.

\*---------------------------------------------------------------------------*/

#ifndef setSizeFieldMapper_H
#define setSizeFieldMapper_H

#include "fieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class setSizeFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class setSizeFieldMapper
:
    public fieldMapper
{
    // Private Data

        //- The size the mapper will set the fields to have
        const label size_;


    // Private Member Functions

        //- Map a field in-place
        template<class Type>
        void map(Field<Type>& f, const Field<Type>&) const;

        //- Map a field and return the result
        template<class Type>
        tmp<Field<Type>> map(const Field<Type>&) const;


public:

    // Constructors

        //- Construct given size
        setSizeFieldMapper(const label size)
        :
            size_(size)
        {}


    //- Destructor
    virtual ~setSizeFieldMapper()
    {}


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, );

        //- Map a label field
        DEFINE_FIELD_MAPPER_OPERATOR(label, );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
