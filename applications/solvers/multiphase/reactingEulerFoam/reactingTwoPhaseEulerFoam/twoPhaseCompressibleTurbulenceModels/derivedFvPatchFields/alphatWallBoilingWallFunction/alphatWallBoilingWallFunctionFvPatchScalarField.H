/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::compressible::alphatWallBoilingWallFunctionFvPatchScalarField

Group
    grpCmpWallFunctions

Description
    A thermal wall function for simulation of subcooled nucleate wall boiling.

    Implements a version of the well-known RPI wall boiling model
    (Kurul & Podowski, 1991). The model implementation is similar to the model
    described by Peltola & Pättikangas (2012).

    References:
    \verbatim
        "On the modeling of multidimensional effects in boiling channels"
        Kurul, N., Podowski, M.Z.,
        ANS Proceedings, National Heat Transfer Conference,
        Minneapolis, Minnesota, USA, July 28-31, 1991,
        ISBN: 0-89448-162-1, pp. 30-40
    \endverbatim

    \verbatim
        "Development and validation of a boiling model for OpenFOAM
        multiphase solver"
        Peltola, J., Pättikangas, T.J.H.,
        CFD4NRS-4 Conference Proceedings, paper 59,
        Daejeon, Korea, September 10-12 2012
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    alphatWallBoilingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressiblealphatWallBoilingWallFunctionFvPatchScalarField_H
#define compressiblealphatWallBoilingWallFunctionFvPatchScalarField_H

#include "alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
            Class alphatWallBoilingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatWallBoilingWallFunctionFvPatchScalarField
:
    public alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
{
    // Private data

        //- dmdt relaxationFactor
        scalar relax_;

        //- Patch face area by cell volume
        scalarField AbyV_;

        //- Convective turbulent thermal diffusivity
        scalarField alphatConv_;


public:

    //- Runtime type information
    TypeName("compressible::alphatWallBoilingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatWallBoilingWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatWallBoilingWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        alphatWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatWallBoilingWallFunctionFvPatchScalarField&,
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
                new alphatWallBoilingWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
