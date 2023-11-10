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
    Foam::nonConformalProcessorCyclicPointPatchField

Description
    Constraint pointPatchField for nonConformalProcessorCyclic couplings.

SourceFiles
    nonConformalProcessorCyclicPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicPointPatchField_H
#define nonConformalProcessorCyclicPointPatchField_H

#include "coupledPointPatchField.H"
#include "nonConformalProcessorCyclicPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nonConformalProcessorCyclicPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalProcessorCyclicPointPatchField
:
    public coupledPointPatchField<Type>
{
    // Private Data

        //- Local reference to processor patch
        const nonConformalProcessorCyclicPointPatch& procPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicPointPatch::typeName_());


    // Constructors

        //- Construct from patch and nonConformalProcessorCyclic field
        nonConformalProcessorCyclicPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, nonConformalProcessorCyclic field and
        //  dictionary
        nonConformalProcessorCyclicPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        nonConformalProcessorCyclicPointPatchField
        (
            const nonConformalProcessorCyclicPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        nonConformalProcessorCyclicPointPatchField
        (
            const nonConformalProcessorCyclicPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new nonConformalProcessorCyclicPointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }

        //- Construct and return a clone onto a new patch
        //  setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const pointPatch& patch,
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new nonConformalProcessorCyclicPointPatchField<Type>
                (
                    patch,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~nonConformalProcessorCyclicPointPatchField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformalProcessorCyclicPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
