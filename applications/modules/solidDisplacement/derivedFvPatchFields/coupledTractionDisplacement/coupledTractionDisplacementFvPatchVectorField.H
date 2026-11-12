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
    Foam::coupledTractionDisplacementFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation in which the traction is caused by
    the pressure of a coupled fluid solver.

    If the coupled fluid region is incompressible and the pressure is kinematic
    the constant uniform \c rhoFluid must be provided to convert the pressure
    to dynamic but for compressible fluids the \c rhoFluid should not be
    specified and defaults to 1;

SourceFiles
    coupledTractionDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledTractionDisplacementFvPatchVectorField_H
#define coupledTractionDisplacementFvPatchVectorField_H

#include "tractionDisplacementFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class coupledTractionDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledTractionDisplacementFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{
    // Private Data

        //- Density of surrounding fluid
        scalar rhoFluid_;


public:

    //- Runtime type information
    TypeName("coupledTractionDisplacement");


    // Constructors

        //- Construct from patch, internal field and dictionary
        coupledTractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  coupledTractionDisplacementFvPatchVectorField onto a new patch
        coupledTractionDisplacementFvPatchVectorField
        (
            const coupledTractionDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        coupledTractionDisplacementFvPatchVectorField
        (
            const coupledTractionDisplacementFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        coupledTractionDisplacementFvPatchVectorField
        (
            const coupledTractionDisplacementFvPatchVectorField&,
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
                new coupledTractionDisplacementFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using tractionDisplacementFvPatchVectorField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
