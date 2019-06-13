/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::cyclicFvsPatchField

Description
    Foam::cyclicFvsPatchField

SourceFiles
    cyclicFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicFvsPatchField_H
#define cyclicFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private Data

        //- Local reference cast into the cyclic patch
        const cyclicFvPatch& cyclicPatch_;


public:

    //- Runtime type information
    TypeName(cyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cyclicFvsPatchField onto a new patch
        cyclicFvsPatchField
        (
            const cyclicFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        cyclicFvsPatchField
        (
            const cyclicFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type>> clone() const
        {
            return tmp<fvsPatchField<Type>>
            (
                new cyclicFvsPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        cyclicFvsPatchField
        (
            const cyclicFvsPatchField<Type>&,
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
                new cyclicFvsPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cyclicFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
