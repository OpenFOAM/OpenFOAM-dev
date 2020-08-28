/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
        void setSize(Field<Type>& f) const;

        template<class Type>
        tmp<Field<Type>> setSize() const;


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

        virtual void operator()
        (
            Field<scalar>& f,
            const Field<scalar>& mapF
        ) const;

        virtual void operator()
        (
            Field<vector>& f,
            const Field<vector>& mapF
        ) const;

        virtual void operator()
        (
            Field<sphericalTensor>& f,
            const Field<sphericalTensor>& mapF
        ) const;

        virtual void operator()
        (
            Field<symmTensor>& f,
            const Field<symmTensor>& mapF
        ) const;

        virtual void operator()
        (
            Field<tensor>& f,
            const Field<tensor>& mapF
        ) const;


        virtual tmp<Field<scalar>> operator()
        (
            const Field<scalar>& mapF
        ) const;

        virtual tmp<Field<vector>> operator()
        (
            const Field<vector>& mapF
        ) const;

        virtual tmp<Field<sphericalTensor>> operator()
        (
            const Field<sphericalTensor>& mapF
        ) const;

        virtual tmp<Field<symmTensor>> operator()
        (
            const Field<symmTensor>& mapF
        ) const;

        virtual tmp<Field<tensor>> operator()
        (
            const Field<tensor>& mapF
        ) const;
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
