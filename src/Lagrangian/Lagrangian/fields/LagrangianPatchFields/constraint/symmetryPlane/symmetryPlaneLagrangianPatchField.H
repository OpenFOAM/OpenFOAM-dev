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
    Foam::symmetryPlaneLagrangianPatchField

Description
    A symmetry-plane boundary condition for Lagrangian. Properties are
    transformed to represent an identical element entering the domain with an
    opposite normal velocity and other vector and tensor quantities
    transformed.

Usage
    Example specification:
    <LagrangianPatchName>
    {
        type            symmetryPlane;
    }
    \endverbatim

SourceFiles
    symmetryPlaneLagrangianPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryPlaneLagrangianPatchField_H
#define symmetryPlaneLagrangianPatchField_H

#include "LagrangianPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class symmetryPlaneLagrangianPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class symmetryPlaneLagrangianPatchField
:
    public LagrangianPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("symmetryPlane");


    // Constructors

        //- Construct from patch and internal field
        symmetryPlaneLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&
        );

        //- Construct from patch, internal field and dictionary
        symmetryPlaneLagrangianPatchField
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy
        symmetryPlaneLagrangianPatchField
        (
            const symmetryPlaneLagrangianPatchField<Type>&
        );

        //- Construct as copy setting the internal field reference
        symmetryPlaneLagrangianPatchField
        (
            const symmetryPlaneLagrangianPatchField<Type>&,
            const regIOobject&
        );

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchField<Type>> clone() const
        {
            return autoPtr<LagrangianPatchField<Type>>
            (
                new symmetryPlaneLagrangianPatchField<Type>
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
                new symmetryPlaneLagrangianPatchField<Type>
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
    #include "symmetryPlaneLagrangianPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
