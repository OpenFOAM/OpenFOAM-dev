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
    Foam::incompressible::
        turbulentMixingLengthDissipationRateInletFvPatchScalarField

Group
    grpIcoRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a turbulence dissipation, \f$\epsilon\f$
    (epsilon) inlet condition based on a specified mixing length.  The patch
    values are calculated using:

        \f[
            \epsilon_p = \frac{C_{\mu}^{0.75} k^{1.5}}{L}
        \f]
   
    where

    \vartable
        \epsilon_p | patch epsilon values
        C_{\mu} | Model coefficient, set to 0.09
        k       | turbulence kinetic energy
        L       | length scale
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        mixingLength | Length scale [m]        | yes         |
        phi          | flux field name         | no          | phi
        k            | turbulence kinetic energy field name | no | k
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            turbulentMixingLengthDissipationRateInlet;
        mixingLength    0.005;
        value           uniform 200;   // placeholder
    }
    \endverbatim

Note
    In the event of reverse flow, a zero-gradient condition is applied

SeeAlso
    Foam::inletOutletFvPatchField

SourceFiles
    turbulentMixingLengthDissipationRateInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleturbulentMixingLengthDissipationRateInlet_H
#define incompressibleturbulentMixingLengthDissipationRateInlet_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
 Class turbulentMixingLengthDissipationRateInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentMixingLengthDissipationRateInletFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- turbulent length scale
        scalar mixingLength_;

        //- Name of the turbulent kinetic energy field
        word kName_;


public:

    //- Runtime type information
    TypeName("turbulentMixingLengthDissipationRateInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentMixingLengthDissipationRateInletFvPatchScalarField
        //  onto a new patch
        turbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const turbulentMixingLengthDissipationRateInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const turbulentMixingLengthDissipationRateInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentMixingLengthDissipationRateInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const turbulentMixingLengthDissipationRateInletFvPatchScalarField&,
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
                new turbulentMixingLengthDissipationRateInletFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
