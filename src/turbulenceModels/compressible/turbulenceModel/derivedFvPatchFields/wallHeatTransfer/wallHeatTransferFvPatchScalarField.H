/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::wallHeatTransferFvPatchScalarField

Group
    grpThermoBoundaryConditions grpWallBoundaryConditions

Description
    This boundary condition provides an enthalpy condition for wall heat
    transfer

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        Tinf         | wall temperature        | yes         |
        alphaWall    | thermal diffusivity     | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            wallHeatTransfer;
        Tinf            uniform 500;
        alphaWall       uniform 1;
    }
    \endverbatim

SourceFiles
    wallHeatTransferFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wallHeatTransferFvPatchScalarField_H
#define wallHeatTransferFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class wallHeatTransferFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wallHeatTransferFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Temperature at the wall
        scalarField Tinf_;

        //- Thermal diffusivity at the wall
        scalarField alphaWall_;


public:

    //- Runtime type information
    TypeName("wallHeatTransfer");


    // Constructors

        //- Construct from patch and internal field
        wallHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wallHeatTransferFvPatchScalarField
        // onto a new patch
        wallHeatTransferFvPatchScalarField
        (
            const wallHeatTransferFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wallHeatTransferFvPatchScalarField
        (
            const wallHeatTransferFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new wallHeatTransferFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wallHeatTransferFvPatchScalarField
        (
            const wallHeatTransferFvPatchScalarField&,
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
                new wallHeatTransferFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return Tinf
            const scalarField& Tinf() const
            {
                return Tinf_;
            }

            //- Return reference to Tinf to allow adjustment
            scalarField& Tinf()
            {
                return Tinf_;
            }

            //- Return alphaWall
            const scalarField& alphaWall() const
            {
                return alphaWall_;
            }

            //- Return reference to alphaWall to allow adjustment
            scalarField& alphaWall()
            {
                return alphaWall_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

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
