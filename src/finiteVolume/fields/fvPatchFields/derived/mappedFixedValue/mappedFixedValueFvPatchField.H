/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::mappedFixedValueFvPatchField

Description
    This boundary condition maps the value at a set of cells or patch faces
    back to *this.

    The sample mode is set by the underlying mapping engine, provided by the
    mappedPatchBase class.

Usage
    \table
        Property     | Description             | Required    | Default value
        fieldName    | name of field to be mapped | no       | this field name
        setAverage   | flag to activate setting of average value | yes |
        average      | average value to apply if \c setAverage = yes | yes |
        interpolationScheme | type of interpolation scheme | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mapped;
        fieldName       T;
        setAverage      no;
        average         0;
        interpolationScheme cell;
        value           uniform 0;
    }
    \endverbatim

    When employing the \c nearestCell sample mode, the user must also specify
    the interpolation scheme using the \c interpolationScheme entry.

    In case of interpolation (where scheme != cell) the limitation is that
    there is only one value per cell.  For example, if you have a cell with two
    boundary faces and both faces sample into the cell, both faces will get the
    same value.

Note
    It is not possible to sample internal faces since volume fields are not
    defined on faces.

See also
    Foam::mappedPatchBase
    Foam::interpolation
    Foam::fixedValueFvPatchField

SourceFiles
    mappedFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFixedValueFvPatchField_H
#define mappedFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
//#include "interpolation.H"
#include "mappedPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mappedFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>,
    public mappedPatchFieldBase<Type>
{

protected:

    // Protected Member Functions

        const mappedPatchBase& mapper
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF
        );


public:

    //- Runtime type information
    TypeName("mapped");


    // Constructors

        //- Construct from patch and internal field
        mappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mappedFixedValueFvPatchField
        //  onto a new patch
        mappedFixedValueFvPatchField
        (
            const mappedFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        mappedFixedValueFvPatchField
        (
            const mappedFixedValueFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        mappedFixedValueFvPatchField
        (
            const mappedFixedValueFvPatchField<Type>&,
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
                new mappedFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

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
    #include "mappedFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
