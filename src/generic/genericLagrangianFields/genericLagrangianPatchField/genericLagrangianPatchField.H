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
    Foam::genericLagrangianPatchField

Description
    This boundary condition provides a generic version of the \c calculated
    condition, useful as a fallback for handling unknown patch types when
    post-processing or running mesh manipulation utilities. Not generally
    applicable as a user-specified condition.

See also
    Foam::calculatedLagrangianPatchField

SourceFiles
    genericLagrangianPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef genericLagrangianPatchField_H
#define genericLagrangianPatchField_H

#include "genericFieldBase.H"
#include "calculatedLagrangianPatchField.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class genericLagrangianPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class genericLagrangianPatchField
:
    public genericFieldBase,
    public calculatedLagrangianPatchField<Type>
{
    // Private Data

        //- Patch field dictionary
        dictionary dict_;


public:

    //- Runtime type information
    TypeName("generic");


    // Constructors

        //- Construct from patch, internal field and dictionary
        genericLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy
        genericLagrangianPatchField
        (
            const genericLagrangianPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        genericLagrangianPatchField
        (
            const genericLagrangianPatchField<Type>&,
            const regIOobject&
        );

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchField<Type>> clone() const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new calculatedLagrangianPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<LagrangianPatchField<Type>> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new genericLagrangianPatchField<Type>
                (
                    *this,
                    iIo
                )
            );
        }


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "genericLagrangianPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
