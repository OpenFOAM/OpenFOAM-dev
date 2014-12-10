/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::incompressible::nutkWallFunctionFvPatchScalarField

Group
    grpIcoWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions, based on turbulence kinetic energy.
    - replicates OpenFOAM v1.5 (and earlier) behaviour

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            nutkWallFunction;
    }
    \endverbatim

SeeAlso
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutkWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutkWallFunctionFvPatchScalarField_H
#define nutkWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
             Class nutkWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutkWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;


public:

    //- Runtime type information
    TypeName("nutkWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutkWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutkWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutkWallFunctionFvPatchScalarField
        //  onto a new patch
        nutkWallFunctionFvPatchScalarField
        (
            const nutkWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutkWallFunctionFvPatchScalarField
        (
            const nutkWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutkWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutkWallFunctionFvPatchScalarField
        (
            const nutkWallFunctionFvPatchScalarField&,
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
                new nutkWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Calculate and return the yPlus at the boundary
        virtual tmp<scalarField> yPlus() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
