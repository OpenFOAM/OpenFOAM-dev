/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::mappedInternalValueFvPatchField

Description
    This boundary condition maps the values from a internal cells to this
    patch.

Usage
    \table
        Property     | Description                | Required | Default value
        field        | name of field to be mapped | no       | this field name
        setAverage   | set the average value?     | no       | yes if average \\
                                                               is specified, \\
                                                               no otherwise
        average      | average value to apply     | if setAverage is true |
        interplolationScheme | the interpolation scheme to use | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedInternalValue;
        field           T;
        average         300;
        interpolationScheme cellPoint;
        value           uniform 300;
    }
    \endverbatim

    This boundary condition will usually be applied to a patch which is of
    mappedInternalPatchBase type, and which holds all the necessary mapping
    information. It can also create its own mapping data which overrides that
    in the mapped patch, or so that it can be applied to a non-mapped patch.
    This is triggered by the presence of controls relating to
    mappedInternalPatchBase (i.e., neighbourRegion, offset, distance, etc ...).

See also
    Foam::mappedInternalPatchBase
    Foam::interpolation
    Foam::fixedValueFvPatchField

SourceFiles
    mappedInternalValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedInternalValueFvPatchField_H
#define mappedInternalValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "mappedInternalPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mappedInternalValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedInternalValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
protected:

    // Protected Member Data

        //- Name of the field to map
        const word fieldName_;

        //- If true adjust the mapped field to maintain an average value
        const bool setAverage_;

        //- Average value the mapped field is adjusted to
        const Type average_;

        //- Interpolation scheme to use
        const word interpolationScheme_;

        //- Sampling engine
        autoPtr<mappedInternalPatchBase> mapperPtr_;


    // Protected Member Functions

        //- Return the mapping engine
        const mappedInternalPatchBase& mapper() const;


public:

    //- Runtime type information
    TypeName("mappedInternalValue");


    // Constructors

        //- Construct from patch, internal field and dictionary
        mappedInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mappedInternalValueFvPatchField
        //  onto a new patch
        mappedInternalValueFvPatchField
        (
            const mappedInternalValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        mappedInternalValueFvPatchField
        (
            const mappedInternalValueFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        mappedInternalValueFvPatchField
        (
            const mappedInternalValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedInternalValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map
            (
                const fvPatchField<Type>&,
                const fvPatchFieldMapper&
            );

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchField<Type>&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedInternalValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
