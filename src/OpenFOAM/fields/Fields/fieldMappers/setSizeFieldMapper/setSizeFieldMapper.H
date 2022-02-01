/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2022 OpenFOAM Foundation
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
    Field sizing mapper which sets the field size and does not map values.

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
    virtual public fieldMapper
{
    // Private member data

        //- The size the mapper will set the fields to have
        label size_;


    // Private Member Functions

        template<class Type>
        void map(Field<Type>& f, const Field<Type>&) const;

        template<class Type>
        tmp<Field<Type>> map(const Field<Type>&) const;


public:

    // Constructors

        //- Construct given size
        setSizeFieldMapper(const label size);


    //- Destructor
    virtual ~setSizeFieldMapper()
    {}


    // Member Functions

        //- All values are unmapped
        virtual bool hasUnmapped() const
        {
            return true;
        }


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_OPERATOR, );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "setSizeFieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
