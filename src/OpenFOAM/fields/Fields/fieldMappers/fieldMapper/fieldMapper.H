/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::fieldMapper

Description
    Abstract base class for field mapping

\*---------------------------------------------------------------------------*/

#ifndef fieldMapper_H
#define fieldMapper_H

#include "Field.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fieldMapper Declaration
\*---------------------------------------------------------------------------*/

class fieldMapper
{

public:

    // Constructors

        //- Null constructor
        fieldMapper()
        {}


    //- Destructor
    virtual ~fieldMapper()
    {}


    // Member Functions

        //- Are there unmapped values? I.e. do all size() elements get
        //  get value
        virtual bool hasUnmapped() const = 0;


    // Member Operators

        virtual void operator()
        (
            Field<scalar>& f,
            const Field<scalar>& mapF
        ) const = 0;

        virtual void operator()
        (
            Field<vector>& f,
            const Field<vector>& mapF
        ) const = 0;

        virtual void operator()
        (
            Field<sphericalTensor>& f,
            const Field<sphericalTensor>& mapF
        ) const = 0;

        virtual void operator()
        (
            Field<symmTensor>& f,
            const Field<symmTensor>& mapF
        ) const = 0;

        virtual void operator()
        (
            Field<tensor>& f,
            const Field<tensor>& mapF
        ) const = 0;

        template<class Type>
        void operator()(Field<Type>& f, const tmp<Field<Type>>& tmapF) const;


        virtual tmp<Field<scalar>> operator()
        (
            const Field<scalar>& mapF
        ) const = 0;

        virtual tmp<Field<vector>> operator()
        (
            const Field<vector>& mapF
        ) const = 0;

        virtual tmp<Field<sphericalTensor>> operator()
        (
            const Field<sphericalTensor>& mapF
        ) const = 0;

        virtual tmp<Field<symmTensor>> operator()
        (
            const Field<symmTensor>& mapF
        ) const = 0;

        virtual tmp<Field<tensor>> operator()
        (
            const Field<tensor>& mapF
        ) const = 0;

        template<class Type>
        tmp<Field<Type>> operator()(const tmp<Field<Type>>& tmapF) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
