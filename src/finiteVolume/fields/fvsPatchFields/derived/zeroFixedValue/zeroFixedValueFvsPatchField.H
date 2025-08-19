/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::zeroFixedValueFvsPatchField

Description
    This boundary condition provides a zero fixed value condition.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            zeroFixedValue;
    }
    \endverbatim

See also
    Foam::fixedValueFvsPatchField

SourceFiles
    zeroFixedValueFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef zeroFixedValueFvsPatchField_H
#define zeroFixedValueFvsPatchField_H

#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class zeroFixedValueFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zeroFixedValueFvsPatchField
:
    public fixedValueFvsPatchField<Type>
{
public:

    //- Runtime type information
    TypeName("zeroFixedValue");


    // Constructors

        //- Construct from patch and internal field
        zeroFixedValueFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroFixedValueFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given zeroFixedValueFvsPatchField
        //  onto a new patch
        zeroFixedValueFvsPatchField
        (
            const zeroFixedValueFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        zeroFixedValueFvsPatchField
        (
            const zeroFixedValueFvsPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        zeroFixedValueFvsPatchField
        (
            const zeroFixedValueFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type>>
            (
                new zeroFixedValueFvsPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvsPatchField onto this fvsPatchField
            virtual void map(const fvsPatchField<Type>&, const fieldMapper&);


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvsPatchField<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zeroFixedValueFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
