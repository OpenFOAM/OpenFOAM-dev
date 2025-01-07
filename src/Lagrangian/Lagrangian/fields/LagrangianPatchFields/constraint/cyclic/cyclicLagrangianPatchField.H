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
    Foam::cyclicLagrangianPatchField

Description
    A cyclic boundary condition for Lagrangian. Properties are transformed by
    the transformation of the associated cyclic patch.

Usage
    Example specification:
    <LagrangianPatchName>
    {
        type            cyclic;
    }
    \endverbatim

SourceFiles
    cyclicLagrangianPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicLagrangianPatchField_H
#define cyclicLagrangianPatchField_H

#include "LagrangianPatchField.H"
#include "cyclicLagrangianPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cyclicLagrangianPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicLagrangianPatchField
:
    public LagrangianPatchField<Type>
{
private:

    // Private Data

        //- Reference to the cyclic patch
        const cyclicLagrangianPatch& cyclicPatch_;


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from patch and internal field
        cyclicLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct from patch, internal field and dictionary
        cyclicLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy
        cyclicLagrangianPatchField
        (
            const cyclicLagrangianPatchField<Type>&
        );

        //- Construct as copy setting the internal field reference
        cyclicLagrangianPatchField
        (
            const cyclicLagrangianPatchField<Type>&,
            const regIOobject&
        );

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchField<Type>> clone() const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new cyclicLagrangianPatchField<Type>
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
                new cyclicLagrangianPatchField<Type>
                (
                    *this,
                    iIo
                )
            );
        }


    // Member Functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            PstreamBuffers&,
            const LagrangianScalarInternalDynamicField& fraction
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cyclicLagrangianPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
