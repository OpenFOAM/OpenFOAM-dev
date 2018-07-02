/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::nonuniformTransformCyclicFvPatchField

Description
    This boundary condition enforces a cyclic condition between a pair of
    boundaries, incorporating a non-uniform transformation.

SourceFiles
    nonuniformTransformCyclicFvPatchField.C
    nonuniformTransformCyclicFvPatchFields.H
    nonuniformTransformCyclicFvPatchFields.C
    nonuniformTransformCyclicFvPatchFieldsFwd.H

\*---------------------------------------------------------------------------*/

#ifndef nonuniformTransformCyclicFvPatchField_H
#define nonuniformTransformCyclicFvPatchField_H

#include "cyclicFvPatchField.H"
#include "nonuniformTransformCyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class nonuniformTransformCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonuniformTransformCyclicFvPatchField
:
    public cyclicFvPatchField<Type>
{
    // Private data

public:

    //- Runtime type information
    TypeName(nonuniformTransformCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonuniformTransformCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonuniformTransformCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nonuniformTransformCyclicFvPatchField onto a new patch
        nonuniformTransformCyclicFvPatchField
        (
            const nonuniformTransformCyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nonuniformTransformCyclicFvPatchField
        (
            const nonuniformTransformCyclicFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new nonuniformTransformCyclicFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nonuniformTransformCyclicFvPatchField
        (
            const nonuniformTransformCyclicFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new nonuniformTransformCyclicFvPatchField<Type>(*this, iF)
            );
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonuniformTransformCyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
