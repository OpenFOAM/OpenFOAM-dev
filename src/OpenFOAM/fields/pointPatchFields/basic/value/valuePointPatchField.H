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
    Foam::valuePointPatchField

Description
    Foam::valuePointPatchField

SourceFiles
    valuePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef valuePointPatchField_H
#define valuePointPatchField_H

#include "pointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class valuePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class valuePointPatchField
:
    public pointPatchField<Type>,
    public Field<Type>
{
    // Private Member Functions

        void checkFieldSize() const;


public:

    //- Runtime type information
    TypeName("value");


    // Constructors

        //- Construct from patch and internal field
        valuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        valuePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        valuePointPatchField
        (
            const valuePointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new valuePointPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        valuePointPatchField
        (
            const valuePointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Copy constructor
        valuePointPatchField(const valuePointPatchField<Type>&) = default;

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new valuePointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return Field<Type>::size();
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const pointPatchFieldMapper&);

            //- Reverse map the given pointPatchField onto this pointPatchField
            //  Used to reconstruct fields
            virtual void rmap(const pointPatchField<Type>&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        // Assignment operators

            virtual void operator=
            (
                const valuePointPatchField<Type>&
            );

            virtual void operator=
            (
                const pointPatchField<Type>&
            );

            virtual void operator=(const Field<Type>&);
            virtual void operator=(const Type&);


        // Force an assignment irrespective of form of patch

            virtual void operator==
            (
                const valuePointPatchField<Type>&
            );

            virtual void operator==
            (
                const pointPatchField<Type>&
            );

            virtual void operator==(const Field<Type>&);
            virtual void operator==(const Type&);
};


// This function is added to override the hack in pointPatchField.H
// which enables simple backward compatibility with versions using
// referenceLevel in GeometicField
template<class Type>
tmp<Field<Type>> operator+
(
    const valuePointPatchField<Type>& vsppf,
    const Type& t
)
{
    return static_cast<const Field<Type>&>(vsppf) + t;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "valuePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
