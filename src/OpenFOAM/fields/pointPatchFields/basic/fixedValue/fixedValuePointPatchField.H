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
    Foam::fixedValuePointPatchField

Description
    A FixedValue boundary condition for pointField.

SourceFiles
    fixedValuePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValuePointPatchField_H
#define fixedValuePointPatchField_H

#include "valuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fixedValuePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedValuePointPatchField
:
    public valuePointPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("fixedValue");


    // Constructors

        //- Construct from patch and internal field
        fixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        fixedValuePointPatchField
        (
            const fixedValuePointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        fixedValuePointPatchField
        (
            const fixedValuePointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new fixedValuePointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Access

            //- Return true if this patch field fixes a value
            virtual bool fixesValue() const
            {
                return true;
            }


    // Member Operators

        // Disable assignment operators

            virtual void operator=(const Field<Type>&){}
            virtual void operator=(const Type&){}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fixedValuePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
