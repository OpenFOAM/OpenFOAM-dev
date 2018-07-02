/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::activePressureForceBaffleVelocityFvPatchVectorField

Description
    This boundary condition is applied to the flow velocity, to simulate the
    opening or closure of a baffle due to local pressure or force changes,
    by merging the behaviours of wall and cyclic conditions.

    The baffle joins two mesh regions, where the open fraction determines
    the interpolation weights applied to each cyclic- and neighbour-patch
    contribution. This means that this is boundary condition is meant to be
    used in an extra wall beyond an existing cyclic patch pair. See PDRMesh
    for more details.

    Once the threshold is crossed, this condition activated and continues to
    open or close at a fixed rate using

        \f[
            x = x_{old} + s \times \frac{dt}{DT}
        \f]

    where

    \vartable
        x       | baffle open fraction [0-1]
        x_{old} | baffle open fraction on previous evaluation
        s       | sign for orientation: 1 to open or -1 to close
        dt      | simulation time step
        DT      | time taken to open the baffle
    \endvartable

    The open fraction is then applied to scale the patch areas.

Usage
    \table
        Property     | Description             | Required    | Default value
        p            | pressure field name     | no          | p
        cyclicPatch  | cyclic patch name       | yes         |
        orientation  | 1 to open or -1 to close | yes|
        openFraction | current open fraction [0-1] | yes |
        openingTime  | time taken to open or close the baffle | yes |
        maxOpenFractionDelta | max fraction change per timestep | yes |
        minThresholdValue | minimum absolute pressure or
                            force difference for activation | yes |
        forceBased   | force (true) or pressure-based (false) activation | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            activePressureForceBaffleVelocity;
        p               p;
        cyclicPatch     cyclic1;
        orientation     1;
        openFraction    0.2;
        openingTime     5.0;
        maxOpenFractionDelta 0.1;
        minThresholdValue 0.01;
        forceBased      false;
    }
    \endverbatim

SourceFiles
    activePressureForceBaffleVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef activePressureForceBaffleVelocityFvPatchVectorField_H
#define activePressureForceBaffleVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class activePressureForceBaffleVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class activePressureForceBaffleVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of the pressure field used to calculate the force
        //  on the active baffle
        word pName_;

        //- Name of the cyclic patch used when the active baffle is open
        word cyclicPatchName_;

        //- Index of the cyclic patch used when the active baffle is open
        label cyclicPatchLabel_;

        //- Orientation (1 or -1) of the active baffle mode
        //  Used to change the direction of opening or closing the baffle
        label orientation_;

        //- Initial wall patch areas
        vectorField initWallSf_;

        //- Initial cyclic patch areas
        vectorField initCyclicSf_;

        //- Initial neighbour-side cyclic patch areas
        vectorField nbrCyclicSf_;

        //- Current fraction of the active baffle which is open
        scalar openFraction_;

        //- Time taken for the active baffle to open
        scalar openingTime_;

        //- Maximum fractional change to the active baffle openness
        //  per time-step
        scalar maxOpenFractionDelta_;

        label curTimeIndex_;

        //- Minimum value for the active baffle to start opening
        scalar minThresholdValue_;

        //- Force based active baffle
        bool fBased_;

        //- Baffle is activated
        bool baffleActivated_;


public:

    //- Runtime type information
    TypeName("activePressureForceBaffleVelocity");


    // Constructors

        //- Construct from patch and internal field
        activePressureForceBaffleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        activePressureForceBaffleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping
        activePressureForceBaffleVelocityFvPatchVectorField
        (
            const activePressureForceBaffleVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        activePressureForceBaffleVelocityFvPatchVectorField
        (
            const activePressureForceBaffleVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new activePressureForceBaffleVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        activePressureForceBaffleVelocityFvPatchVectorField
        (
            const activePressureForceBaffleVelocityFvPatchVectorField&,
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
                new activePressureForceBaffleVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


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
