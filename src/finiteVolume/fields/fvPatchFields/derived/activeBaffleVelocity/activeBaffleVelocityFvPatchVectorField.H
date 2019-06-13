/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::activeBaffleVelocityFvPatchVectorField

Description
    This velocity boundary condition simulates the opening of a baffle due
    to local flow conditions, by merging the behaviours of wall and cyclic
    conditions.  The baffle joins two mesh regions, where the open fraction
    determines the interpolation weights applied to each cyclic- and
    neighbour-patch contribution.

    We determine whether the baffle is opening or closing from the sign of
    the net force across the baffle, from which the baffle open fraction is
    updated using:

        \f[
            x = x_{old} + sign(F_{net})\frac{dt}{DT}
        \f]

    where

    \vartable
        x       | baffle open fraction [0-1]
        x_{old} | baffle open fraction on previous evaluation
        dt      | simulation time step
        DT      | time taken to open the baffle
        F_{net} | net force across the baffle
    \endvartable

    The open fraction is then applied to scale the patch areas.

Usage
    \table
        Property     | Description             | Required    | Default value
        p            | pressure field name     | no          | p
        cyclicPatch  | cylclic patch name      | yes         |
        orientation  | 1 or -1 used to switch flow direction | yes|
        openFraction | current opatch open fraction [0-1]| yes |
        openingTime  | time taken to open the baffle | yes |
        maxOpenFractionDelta | max open fraction change per timestep | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            activeBaffleVelocity;
        p               p;
        cyclicPatch     cyclic1;
        orientation     1;
        openFraction    0.2;
        openingTime     5.0;
        maxOpenFractionDelta 0.1;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::cyclicFvPatchField

SourceFiles
    activeBaffleVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef activeBaffleVelocityFvPatchVectorField_H
#define activeBaffleVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class activeBaffleVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class activeBaffleVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Name of the pressure field used to calculate the force
        //  on the active baffle
        word pName_;

        //- Name of the cyclic patch used when the active baffle is open
        word cyclicPatchName_;

        //- Index of the cyclic patch used when the active baffle is open
        label cyclicPatchLabel_;

        //- Orientation (1 or -1) of the active baffle patch.
        //  Used to change the direction of opening without the need for
        //  reordering the patch faces
        label orientation_;

        //- Initial wall patch areas
        vectorField initWallSf_;

        //- Initial this-side cyclic patch areas
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


public:

    //- Runtime type information
    TypeName("activeBaffleVelocity");


    // Constructors

        //- Construct from patch and internal field
        activeBaffleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        activeBaffleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given activeBaffleVelocityFvPatchVectorField
        //  onto a new patch
        activeBaffleVelocityFvPatchVectorField
        (
            const activeBaffleVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        activeBaffleVelocityFvPatchVectorField
        (
            const activeBaffleVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new activeBaffleVelocityFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        activeBaffleVelocityFvPatchVectorField
        (
            const activeBaffleVelocityFvPatchVectorField&,
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
                new activeBaffleVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchVectorField&, const labelList&);


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
