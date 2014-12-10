/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::nonuniformTransformCyclicFvsPatchField

Description
    Foam::nonuniformTransformCyclicFvsPatchField

SourceFiles
    nonuniformTransformCyclicFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonuniformTransformCyclicFvsPatchField_H
#define nonuniformTransformCyclicFvsPatchField_H

#include "cyclicFvsPatchField.H"
#include "nonuniformTransformCyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class nonuniformTransformCyclicFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonuniformTransformCyclicFvsPatchField
:
    public cyclicFvsPatchField<Type>
{

public:

    //- Runtime type information
    TypeName(nonuniformTransformCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonuniformTransformCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonuniformTransformCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nonuniformTransformCyclicFvsPatchField onto a new patch
        nonuniformTransformCyclicFvsPatchField
        (
            const nonuniformTransformCyclicFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nonuniformTransformCyclicFvsPatchField
        (
            const nonuniformTransformCyclicFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type> > clone() const
        {
            return tmp<fvsPatchField<Type> >
            (
                new nonuniformTransformCyclicFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nonuniformTransformCyclicFvsPatchField
        (
            const nonuniformTransformCyclicFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type> > clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type> >
            (
                new nonuniformTransformCyclicFvsPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "nonuniformTransformCyclicFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
