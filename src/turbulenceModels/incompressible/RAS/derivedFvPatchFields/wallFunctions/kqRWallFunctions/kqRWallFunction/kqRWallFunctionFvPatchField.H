/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::incompressible::kqRWallFunctionFvPatchField

Group
    grpIcoWallFunctions

Description
    This boundary condition provides a suitable condition for turbulence
    \c k, \c q, and \c R fields for the case of high Reynolds number flow using
    wall functions.

    It is a simple wrapper around the zero-gradient condition.

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            kqRWallFunction;
    }
    \endverbatim

SeeAlso
    Foam::zeroGradientFvPatchField

SourceFiles
    kqRWallFunctionFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef kqRWallFunctionFvPatchField_H
#define kqRWallFunctionFvPatchField_H

#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
                Class kqRWallFunctionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class kqRWallFunctionFvPatchField
:
    public zeroGradientFvPatchField<Type>
{

    // Private Member Functions

        //- Check the type of the patch
        void checkType();


public:

    //- Runtime type information
    TypeName("kqRWallFunction");


    // Constructors

        //- Construct from patch and internal field
        kqRWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        kqRWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // kqRWallFunctionFvPatchField
        //  onto a new patch
        kqRWallFunctionFvPatchField
        (
            const kqRWallFunctionFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        kqRWallFunctionFvPatchField
        (
            const kqRWallFunctionFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new kqRWallFunctionFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        kqRWallFunctionFvPatchField
        (
            const kqRWallFunctionFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new kqRWallFunctionFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::Pstream::blocking
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kqRWallFunctionFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
