/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::epsilonWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulence dissipation wall constraint
    for low- and high-Reynolds number turbulence models.

    The condition can be applied to wall boundaries for which it
    - calculates \c epsilon and \c G
    - specifies the near-wall epsilon value

    where

    \vartable
        epsilon | turbulence dissipation field
        G       | turbulence generation field
    \endvartable

    The model switches between laminar and turbulent functions based on the
    laminar-to-turbulent y+ value derived from the kappa and E specified in the
    corresponding nutWallFunction.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            epsilonWallFunction;
    }
    \endverbatim

See also
    Foam::fixedInternalValueFvPatchField
    Foam::omegaWallFunctionFvPatchScalarField

SourceFiles
    epsilonWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonWallFunctionFvPatchScalarField_H
#define epsilonWallFunctionFvPatchScalarField_H

#include "wallCellWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class momentumTransportModel;

/*---------------------------------------------------------------------------*\
            Class epsilonWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class epsilonWallFunctionFvPatchScalarField
:
    public wallCellWallFunctionFvPatchScalarField
{
private:

    // Private Data

        //- Wall cell turbulence generation values. Master patch only.
        autoPtr<scalarField> wallCellGPtr_;

        //- Wall cell turbulence dissipation values. Master patch only.
        autoPtr<scalarField> wallCellEpsilonPtr_;


    // Private Member Functions

        //- Calculate and return G and epsilon
        Pair<tmp<scalarField>> calculate(const momentumTransportModel&) const;

        //- Calculate the turbulent values in wall cells and apply to the
        //  turbulence model's fields
        void updateCoeffsMaster();

        //- Constrain the matrix such that the wall cells take the value
        //  specified by the wall function
        void manipulateMatrixMaster(fvMatrix<scalar>& matrix);


public:

    //- Runtime type information
    TypeName("epsilonWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        epsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  epsilonWallFunctionFvPatchScalarField
        //  onto a new patch
        epsilonWallFunctionFvPatchScalarField
        (
            const epsilonWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        epsilonWallFunctionFvPatchScalarField
        (
            const epsilonWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        epsilonWallFunctionFvPatchScalarField
        (
            const epsilonWallFunctionFvPatchScalarField&,
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
                new epsilonWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~epsilonWallFunctionFvPatchScalarField()
    {}


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
