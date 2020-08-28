/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::cyclicAMILduInterfaceField

Description
    Abstract base class for cyclic AMI coupled interfaces

SourceFiles
    cyclicAMILduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMILduInterfaceField_H
#define cyclicAMILduInterfaceField_H

#include "primitiveFieldsFwd.H"
#include "transformer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cyclicAMILduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMILduInterfaceField
{

public:

    //- Runtime type information
    TypeName("cyclicAMILduInterfaceField");


    // Constructors

        //- Construct null
        cyclicAMILduInterfaceField()
        {}


    //- Destructor
    virtual ~cyclicAMILduInterfaceField();


    // Member Functions

        // Access

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const = 0;

            //- Return rank of component for transform
            virtual int rank() const = 0;

            //- Is the transform required
            bool transforms() const
            {
                return transform().transforms() && rank() != 0;
            }


        //- Transform given patch field
        template<class Type>
        void transformCoupleField(Field<Type>& f) const;

        //- Transform given patch internal field
        void transformCoupleField
        (
            scalarField& psiInternal,
            const direction cmpt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tensorField.H"

template<class Type>
void Foam::cyclicAMILduInterfaceField::transformCoupleField
(
    Field<Type>& f
) const
{
    if (transforms())
    {
        transform().transform(f, f);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
