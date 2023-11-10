/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    Foam::symmetryPlaneFvPatchField

Description
    This boundary condition enforces a symmetryPlane constraint

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            symmetryPlane;
    }
    \endverbatim

SourceFiles
    symmetryPlaneFvPatchField.C
    symmetryPlaneFvPatchFields.C
    symmetryPlaneFvPatchFields.H
    symmetryPlaneFvPatchFieldsFwd.H

\*---------------------------------------------------------------------------*/

#ifndef symmetryPlaneFvPatchField_H
#define symmetryPlaneFvPatchField_H

#include "basicSymmetryFvPatchField.H"
#include "symmetryPlaneFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class symmetryPlaneFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class symmetryPlaneFvPatchField
:
    public basicSymmetryFvPatchField<Type>
{
    // Private Data

        //- Local reference cast into the symmetryPlane patch
        const symmetryPlaneFvPatch& symmetryPlanePatch_;


public:

    //- Runtime type information
    TypeName(symmetryPlaneFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        symmetryPlaneFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        symmetryPlaneFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given symmetryPlaneFvPatchField
        //  onto a new patch
        symmetryPlaneFvPatchField
        (
            const symmetryPlaneFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        symmetryPlaneFvPatchField
        (
            const symmetryPlaneFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        symmetryPlaneFvPatchField
        (
            const symmetryPlaneFvPatchField<Type>&,
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
                new symmetryPlaneFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Return gradient at boundary
        virtual tmp<Field<Type>> snGrad() const;

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Return face-gradient transform diagonal
        virtual tmp<Field<Type>> snGradTransformDiag() const;
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField> symmetryPlaneFvPatchField<scalar>::snGrad() const;

template<>
void symmetryPlaneFvPatchField<scalar>::evaluate
(
    const Pstream::commsTypes commsType
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "symmetryPlaneFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
