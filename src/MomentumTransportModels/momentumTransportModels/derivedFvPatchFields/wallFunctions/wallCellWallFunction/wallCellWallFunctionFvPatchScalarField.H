/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::wallCellWallFunctionFvPatchScalarField

Description
    Base class for wall functions that modify cell values

SourceFiles
    wallCellWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wallCellWallFunctionFvPatchScalarField_H
#define wallCellWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class wallCellWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wallCellWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
private:

    // Private Data

        //- Tolerance used in weighted calculations
        static scalar tol_;

        //- Master patch index
        label masterPatchi_;

        //- Wall cell indices. Master patch only.
        autoPtr<labelList> wallCellsPtr_;

        //- Wall cell fractions; i.e., what proportion of the wall cell value
        //  gets overridden. Master patch only.
        autoPtr<scalarField> wallCellFractionPtr_;


protected:

    // Protected Member Functions

        // Access

            //- Return the master patch index
            inline label masterPatchIndex() const
            {
                return masterPatchi_;
            }

            //- Return the wall cell indices
            inline const labelList& wallCells() const
            {
                return wallCellsPtr_();
            }

            //- Return the wall cell fractions
            inline const scalarField& wallCellFraction() const
            {
                return wallCellFractionPtr_();
            }


        //- Initialise the master cell indices and fractions, and allocate the
        //  turbulence fields
        void initMaster();

        //- Sum values from a patch field into a cell field
        static void patchFieldAddToCellField
        (
            const fvPatch& fvp,
            const scalarField& pf,
            scalarField& vf
        );

        //- Average a set of patch fields into a wall cell field
        tmp<scalarField> patchFieldsToWallCellField
        (
            const PtrList<scalarField>& pfs
        ) const;


public:

    // Constructors

        //- Construct from patch, internal field and dictionary
        wallCellWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired
        );

        //- Construct by mapping given wallCellWallFunctionFvPatchScalarField
        //  onto a new patch
        wallCellWallFunctionFvPatchScalarField
        (
            const wallCellWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&,
            const bool mappingRequired
        );

        //- Disallow copy without setting internal field reference
        wallCellWallFunctionFvPatchScalarField
        (
            const wallCellWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        wallCellWallFunctionFvPatchScalarField
        (
            const wallCellWallFunctionFvPatchScalarField&,
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
                new wallCellWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~wallCellWallFunctionFvPatchScalarField()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
