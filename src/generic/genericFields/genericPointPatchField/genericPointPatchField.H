/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::genericPointPatchField

Description
    This boundary condition provides a generic version of the \c calculated
    condition, useful as a fallback for handling unknown patch types when
    post-processing or running mesh manipulation utilities. Not generally
    applicable as a user-specified condition.

See also
    Foam::calculatedPointPatchField

SourceFiles
    genericPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef genericPointPatchField_H
#define genericPointPatchField_H

#include "genericFieldBase.H"
#include "calculatedPointPatchField.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class genericPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class genericPointPatchField
:
    public genericFieldBase,
    public calculatedPointPatchField<Type>
{
    // Private Data

        //- Patch field dictionary
        dictionary dict_;

        //- Fields
        #define DeclareTypeFields(Type, nullArg) \
            HashPtrTable<Field<Type>> Type##Fields_;
        FOR_ALL_FIELD_TYPES(DeclareTypeFields);
        #undef DeclareTypeFields


public:

    //- Runtime type information
    TypeName("generic");


    // Constructors

        //- Construct from patch, internal field and dictionary
        genericPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given genericPointPatchField onto a new patch
        genericPointPatchField
        (
            const genericPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        genericPointPatchField
        (
            const genericPointPatchField<Type>&,
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
                new genericPointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given pointPatchField onto this pointPatchField
            virtual void map(const pointPatchField<Type>&, const fieldMapper&);

            //- Reset the pointPatchField to the given pointPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const pointPatchField<Type>&);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "genericPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
