/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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
    Foam::epsilonmWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulence dissipation wall constraint
    for the mixtureKEpsilon model

SourceFiles
    epsilonmWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonmWallFunctionFvPatchScalarField_H
#define epsilonmWallFunctionFvPatchScalarField_H

#include "wallCellWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class epsilonmWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class epsilonmWallFunctionFvPatchScalarField
:
    public wallCellWallFunctionFvPatchScalarField
{
private:

    // Private Member Functions

        //- Constrain the matrix such that the wall cells take the value
        //  specified by the wall function
        void manipulateMatrixMaster(fvMatrix<scalar>& matrix);


public:

    //- Runtime type information
    TypeName("epsilonmWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        epsilonmWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  epsilonmWallFunctionFvPatchScalarField
        //  onto a new patch
        epsilonmWallFunctionFvPatchScalarField
        (
            const epsilonmWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        epsilonmWallFunctionFvPatchScalarField
        (
            const epsilonmWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        epsilonmWallFunctionFvPatchScalarField
        (
            const epsilonmWallFunctionFvPatchScalarField&,
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
                new epsilonmWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~epsilonmWallFunctionFvPatchScalarField()
    {}


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);


    // Member Operators

        //- Inherit assignment
        using wallCellWallFunctionFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
