/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::MRFslipFvPatchVectorField

Description
    Rotating wall-velocity condition to be used for a slip-wall rotating with
    the moving frame in an MRF (multi-reference frame) or SRF (single reference
    frame) case.

    SRF cases are simply MRF cases with a single MRF zone which covers the
    entire domain.

Usage
    Example of the boundary condition specification for an SRF case or MRF
    case with a single zone:
    \verbatim
    <patchName>
    {
        type            MRFslip;
    }
    \endverbatim
    or if the case has several MRF zones the particular zone this patch is in
    must be named explicitly, e.g.:
    \verbatim
    <patchName>
    {
        type            MRFslip;
        MRFZoneName     rotor;
    }
    \endverbatim

See also
    Foam::MRFPatchField
    Foam::MRFZone

SourceFiles
    MRFslipFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef MRFslipFvPatchVectorField_H
#define MRFslipFvPatchVectorField_H

#include "MRFPatchField.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class MRFslipFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class MRFslipFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public MRFPatchField
{

public:

    //- Runtime type information
    TypeName("MRFslip");


    // Constructors

        //- Construct from patch, internal field and dictionary
        MRFslipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given MRFslipFvPatchVectorField
        //  onto a new patch
        MRFslipFvPatchVectorField
        (
            const MRFslipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        MRFslipFvPatchVectorField
        (
            const MRFslipFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        MRFslipFvPatchVectorField
        (
            const MRFslipFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new MRFslipFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
