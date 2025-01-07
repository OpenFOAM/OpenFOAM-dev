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
    Foam::internalLagrangianPatchField

Description
    An internal boundary condition for Lagrangian

Usage
    Example specification:
    <LagrangianPatchName>
    {
        type            internal;
    }
    \endverbatim

SourceFiles
    internalLagrangianPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef internalLagrangianPatchField_H
#define internalLagrangianPatchField_H

#include "LagrangianPatchField.H"
#include "internalLagrangianPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class internalLagrangianPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class internalLagrangianPatchField
:
    public LagrangianPatchField<Type>
{
private:

    // Private Data

        //- Reference to the internal patch
        const internalLagrangianPatch& internalPatch_;


public:

    //- Runtime type information
    TypeName("internal");


    // Constructors

        //- Construct from patch and internal field
        internalLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct from patch, internal field and dictionary
        internalLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy
        internalLagrangianPatchField
        (
            const internalLagrangianPatchField<Type>&
        );

        //- Construct as copy setting the internal field reference
        internalLagrangianPatchField
        (
            const internalLagrangianPatchField<Type>&,
            const regIOobject&
        );

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchField<Type>> clone() const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new internalLagrangianPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianPatchField<Type>> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new internalLagrangianPatchField<Type>
                (
                    *this,
                    iIo
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "internalLagrangianPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
