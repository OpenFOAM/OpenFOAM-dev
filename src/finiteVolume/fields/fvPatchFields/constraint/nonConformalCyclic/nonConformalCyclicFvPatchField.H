/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::nonConformalCyclicFvPatchField

Description
    This boundary condition enforces a non-conformal cyclic condition between a
    pair of boundaries.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nonConformalCyclic;
    }
    \endverbatim

Note
    The patches must be topologically similar, i.e. if the owner patch is
    transformed to the neighbour patch, the patches should be identical (or
    very similar).

SourceFiles
    nonConformalCyclicFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCyclicFvPatchField_H
#define nonConformalCyclicFvPatchField_H

#include "cyclicFvPatchField.H"
#include "nonConformalCyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonConformalCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalCyclicFvPatchField
:
    public cyclicFvPatchField<Type>
{
    // Private Data

        //- Reference to the fvPatch
        const nonConformalCyclicFvPatch& nonConformalCyclicFvPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonConformalCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalCyclicFvPatchField onto a
        //  new patch
        nonConformalCyclicFvPatchField
        (
            const nonConformalCyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nonConformalCyclicFvPatchField
        (
            const nonConformalCyclicFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        nonConformalCyclicFvPatchField
        (
            const nonConformalCyclicFvPatchField<Type>&,
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
                new nonConformalCyclicFvPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformalCyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
