/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
    Foam::zeroFixedValuePointPatchField

Description
    Enables the specification of a zero fixed value boundary condition.

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            zeroFixedValue;
    }
    \endverbatim

SourceFiles
    zeroFixedValuePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef zeroFixedValuePointPatchField_H
#define zeroFixedValuePointPatchField_H

#include "fixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class zeroFixedValuePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zeroFixedValuePointPatchField
:
    public fixedValuePointPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("zeroFixedValue");


    // Constructors

        //- Construct from patch and internal field
        zeroFixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroFixedValuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        zeroFixedValuePointPatchField
        (
            const zeroFixedValuePointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        zeroFixedValuePointPatchField
        (
            const zeroFixedValuePointPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        zeroFixedValuePointPatchField
        (
            const zeroFixedValuePointPatchField<Type>&,
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
                new zeroFixedValuePointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zeroFixedValuePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
