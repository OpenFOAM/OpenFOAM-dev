/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::MRFstaticNoSlipFvPatchVectorField

Description
    Wall-velocity condition to be used on walls that are stationary within the
    moving frame in an MRF (multi-reference frame) or SRF (single reference
    frame) case.

    SRF cases are simply MRF cases with a single MRF zone which covers the
    entire domain.

Usage
    Example of the boundary condition specification for an SRF case or MRF
    case with a single zone:
    \verbatim
    <patchName>
    {
        type            MRFstaticNoSlip;
    }
    \endverbatim
    or if the case has several MRF zones the particular zone this patch is in
    must be named explicitly, e.g.:
    \verbatim
    <patchName>
    {
        type            MRFstaticNoSlip;
        MRFZoneName     rotor;
    }
    \endverbatim

See also
    Foam::MRFnoSlip
    Foam::MRFZone
    Foam::noSlip

SourceFiles
    MRFstaticNoSlipFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef MRFstaticNoSlipFvPatchVectorField_H
#define MRFstaticNoSlipFvPatchVectorField_H

#include "MRFPatchField.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class MRFstaticNoSlipFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class MRFstaticNoSlipFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public MRFPatchField
{

public:

    //- Runtime type information
    TypeName("MRFstaticNoSlip");


    // Constructors

        //- Construct from patch, internal field and dictionary
        MRFstaticNoSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given MRFstaticNoSlipFvPatchVectorField
        //  onto a new patch
        MRFstaticNoSlipFvPatchVectorField
        (
            const MRFstaticNoSlipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        MRFstaticNoSlipFvPatchVectorField
        (
            const MRFstaticNoSlipFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        MRFstaticNoSlipFvPatchVectorField
        (
            const MRFstaticNoSlipFvPatchVectorField&,
            const DimensionedField<vector, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new MRFstaticNoSlipFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchVectorField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
