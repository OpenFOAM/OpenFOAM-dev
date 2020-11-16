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
    Foam::fanPressureFvPatchScalarField

Description
    This boundary condition can be applied to assign either a pressure inlet
    or outlet total pressure condition for a fan.

Usage
    \table
        Property     | Description             | Required    | Default value
        fanCurve     | fan curve file name     | yes         |
        p0           | environmental total pressure | yes    |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            fanPressure;
        fanCurve        table;
        fanCurveCoeffs
        {
            file            "$FOAM_CASE/constant/pressureVsQ";
            format          csv;
            nHeaderLine     1;
            refColumn       0;
            componentColumns 1(1);
            separator       ",";
            mergeSeparators no;
            outOfBounds     clamp;
            interpolationScheme linear;
        }
        direction       in;
        p0              uniform 0;
        value           uniform 0;
    }

    outlet
    {
        type            fanPressure;
        fanCurve        table;
        fanCurveCoeffs
        {
            file            "$FOAM_CASE/constant/pressureVsQ";
            format          csv;
            nHeaderLine     1;
            refColumn       0;
            componentColumns 1(1);
            separator       ",";
            mergeSeparators no;
            outOfBounds     clamp;
            interpolationScheme linear;
        }
        direction       out;
        p0              uniform 0;
        value           uniform 0;
    }
    \endverbatim

    The above example shows the use of a comma separated (CSV) file to specify
    the condition.

See also
    Foam::totalPressureFvPatchScalarField
    Foam::Function1s

SourceFiles
   fanPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fanPressureFvPatchScalarField_H
#define fanPressureFvPatchScalarField_H

#include "totalPressureFvPatchScalarField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fanPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fanPressureFvPatchScalarField
:
    public totalPressureFvPatchScalarField
{
public:

    // Public Enumerations

        //- Fan flow direction
        enum fanFlowDirection
        {
            ffdIn,
            ffdOut
        };

        //- Fan flow directions names
        static const NamedEnum<fanFlowDirection, 2> fanFlowDirectionNames_;


private:

    // Private Data

        //- Fan curve
        const autoPtr<Function1<scalar>> fanCurve_;

        //- Direction of flow through the fan relative to patch
        const fanFlowDirection direction_;


public:

    //- Runtime type information
    TypeName("fanPressure");


    // Constructors

        //- Construct from patch and internal field
        fanPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fanPressureFvPatchScalarField
        //  onto a new patch
        fanPressureFvPatchScalarField
        (
            const fanPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        fanPressureFvPatchScalarField
        (
            const fanPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        fanPressureFvPatchScalarField
        (
            const fanPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fanPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
