/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::symmetryPointPatchField

Description
    A Symmetry boundary condition for pointField

SourceFiles
    symmetryPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryPointPatchField_H
#define symmetryPointPatchField_H

#include "basicSymmetryPointPatchField.H"
#include "symmetryPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class symmetryPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class symmetryPointPatchField
:
    public basicSymmetryPointPatchField<Type>
{

public:

    //- Runtime type information
    TypeName(symmetryPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        symmetryPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        symmetryPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        symmetryPointPatchField
        (
            const symmetryPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new symmetryPointPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        symmetryPointPatchField
        (
            const symmetryPointPatchField<Type>&,
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
                new symmetryPointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Return the constraint type this pointPatchField implements
        virtual const word& constraintType() const
        {
            return symmetryPointPatch::typeName;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "symmetryPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
