/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::FieldMapper

Description
    Abstract base class to hold the Field mapping addressing and weights.

\*---------------------------------------------------------------------------*/

#ifndef FieldMapper_H
#define FieldMapper_H

#include "mapDistributeBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FieldMapper Declaration
\*---------------------------------------------------------------------------*/

class FieldMapper
{

public:

    // Constructors

        //- Null constructor
        FieldMapper()
        {}


    //- Destructor
    virtual ~FieldMapper()
    {}


    // Member Functions

        virtual label size() const = 0;

        virtual bool direct() const = 0;

        virtual bool distributed() const
        {
            return false;
        }

        virtual const mapDistributeBase& distributeMap() const
        {
            FatalErrorInFunction
                << "attempt to access null distributeMap"
                << abort(FatalError);
            return *(new mapDistributeBase());
        }

        //- Are there unmapped values? I.e. do all size() elements get
        //  get value
        virtual bool hasUnmapped() const = 0;

        virtual const labelUList& directAddressing() const
        {
            FatalErrorInFunction
                << "attempt to access null direct addressing"
                << abort(FatalError);

            return labelUList::null();
        }

        virtual const labelListList& addressing() const
        {
            FatalErrorInFunction
                << "attempt to access null interpolation addressing"
                << abort(FatalError);

            return labelListList::null();
        }

        virtual const scalarListList& weights() const
        {
            FatalErrorInFunction
                << "attempt to access null interpolation weights"
                << abort(FatalError);

            return scalarListList::null();
        }


    // Member Operators

        template<class Type>
        tmp<Field<Type>> operator()(const Field<Type>& f) const
        {
            return tmp<Field<Type>>(new Field<Type>(f, *this));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
