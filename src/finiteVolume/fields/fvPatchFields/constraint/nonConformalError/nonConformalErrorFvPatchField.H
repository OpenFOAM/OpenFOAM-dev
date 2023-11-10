/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::nonConformalErrorFvPatchField

Description
    This boundary condition enforces a non-conformal error condition.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nonConformalError;
    }
    \endverbatim

SourceFiles
    nonConformalErrorFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalErrorFvPatchField_H
#define nonConformalErrorFvPatchField_H

#include "zeroGradientFvPatchField.H"
#include "nonConformalErrorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonConformalErrorFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalErrorFvPatchField
:
    public zeroGradientFvPatchField<Type>
{
    // Private Data

        //- Reference to the fvPatch
        const nonConformalErrorFvPatch& nonConformalErrorFvPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalErrorFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonConformalErrorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalErrorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalErrorFvPatchField onto a
        //  new patch
        nonConformalErrorFvPatchField
        (
            const nonConformalErrorFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nonConformalErrorFvPatchField
        (
            const nonConformalErrorFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        nonConformalErrorFvPatchField
        (
            const nonConformalErrorFvPatchField<Type>&,
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
                new nonConformalErrorFvPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformalErrorFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
