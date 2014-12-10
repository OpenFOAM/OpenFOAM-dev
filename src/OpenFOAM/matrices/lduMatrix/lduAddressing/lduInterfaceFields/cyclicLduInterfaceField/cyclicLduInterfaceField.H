/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::cyclicLduInterfaceField

Description
    Abstract base class for cyclic coupled interfaces.

SourceFiles
    cyclicLduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicLduInterfaceField_H
#define cyclicLduInterfaceField_H

#include "primitiveFieldsFwd.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicLduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicLduInterfaceField
{

public:

    //- Runtime type information
    TypeName("cyclicLduInterfaceField");


    // Constructors

        //- Construct given coupled patch
        cyclicLduInterfaceField()
        {}


    //- Destructor
    virtual ~cyclicLduInterfaceField();


    // Member Functions

        // Access

            //- Is the transform required
            virtual bool doTransform() const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const = 0;

            //- Return rank of component for transform
            virtual int rank() const = 0;


        //- Transform given patch field
        template<class Type>
        void transformCoupleField(Field<Type>& f) const;

        //- Transform given patch component field
        void transformCoupleField
        (
            scalarField& f,
            const direction cmpt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tensorField.H"

template<class Type>
void Foam::cyclicLduInterfaceField::transformCoupleField
(
    Field<Type>& f
) const
{
    if (doTransform())
    {
        if (forwardT().size() == 1)
        {
            transform(f, forwardT()[0], f);
        }
        else
        {
            transform(f, forwardT(), f);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
