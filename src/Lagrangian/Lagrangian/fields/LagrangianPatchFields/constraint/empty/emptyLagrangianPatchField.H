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
    Foam::emptyLagrangianPatchField

Description
    An empty boundary condition for Lagrangian. This is used for the patches
    that have normals that point outside of the space of the reduced
    dimensionality; e.g., the front and back planes of a two-dimensional slab.

    This operates in a manner identical to a symmetry condition. Properties are
    transformed to represent an identical element entering the domain with an
    opposite normal velocity and other vector and tensor quantities
    transformed.

    Whether or not this is actually used is a matter for the cloud model. If
    the cloud model is also of reduced dimensionality and only moves in the
    two-dimensional plane or the one-dimensional line, then the elements will
    never track into the empty patches and this condition will never be used.
    If, on the other hand, the cloud model is still considered
    three-dimensional even when the Eulerian modelling is of reduced
    dimensionality, then the empty patches will be interacted with and this
    condition will be used.

Usage
    Example specification:
    <LagrangianPatchName>
    {
        type            empty;
    }
    \endverbatim

SourceFiles
    emptyLagrangianPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef emptyLagrangianPatchField_H
#define emptyLagrangianPatchField_H

#include "LagrangianPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class emptyLagrangianPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class emptyLagrangianPatchField
:
    public LagrangianPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("empty");


    // Constructors

        //- Construct from patch and internal field
        emptyLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct from patch, internal field and dictionary
        emptyLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy
        emptyLagrangianPatchField
        (
            const emptyLagrangianPatchField<Type>&
        );

        //- Construct as copy setting the internal field reference
        emptyLagrangianPatchField
        (
            const emptyLagrangianPatchField<Type>&,
            const regIOobject&
        );

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchField<Type>> clone() const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new emptyLagrangianPatchField<Type>
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
                new emptyLagrangianPatchField<Type>
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
    #include "emptyLagrangianPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
