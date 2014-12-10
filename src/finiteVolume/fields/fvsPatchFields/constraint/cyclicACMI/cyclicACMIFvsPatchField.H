/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::cyclicACMIFvsPatchField

Description
    Foam::cyclicACMIFvsPatchField

SourceFiles
    cyclicACMIFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIFvsPatchField_H
#define cyclicACMIFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "cyclicACMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cyclicACMIFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicACMIFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const cyclicACMIFvPatch& cyclicACMIPatch_;


public:

    //- Runtime type information
    TypeName(cyclicACMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicACMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicACMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cyclicACMIFvsPatchField onto a new patch
        cyclicACMIFvsPatchField
        (
            const cyclicACMIFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cyclicACMIFvsPatchField
        (
            const cyclicACMIFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type> > clone() const
        {
            return tmp<fvsPatchField<Type> >
            (
                new cyclicACMIFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cyclicACMIFvsPatchField
        (
            const cyclicACMIFvsPatchField<Type>&,
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
                new cyclicACMIFvsPatchField<Type>(*this, iF)
            );
        }

    // Member functions

        // Access

            //- Return true if running parallel
            virtual bool coupled() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cyclicACMIFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
