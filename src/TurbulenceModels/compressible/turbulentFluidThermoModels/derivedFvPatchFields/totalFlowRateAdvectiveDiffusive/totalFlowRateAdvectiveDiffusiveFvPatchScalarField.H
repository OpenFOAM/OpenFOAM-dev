/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::totalFlowRateAdvectiveDiffusiveFvPatchScalarField

Description
    This BC is used for species inlets. The diffusion and advection fluxes are
    considered to calculate the inlet value for the species
    The massFluxFraction sets the fraction of the flux of each particular
    species.

SourceFiles
    totalFlowRateAdvectiveDiffusiveFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef totalFlowRateAdvectiveDiffusiveFvPatchScalarField_H
#define totalFlowRateAdvectiveDiffusiveFvPatchScalarField_H

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class totalFlowRateAdvectiveDiffusiveFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class totalFlowRateAdvectiveDiffusiveFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
    // Private data

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if necessary
        word rhoName_;

        //- Mass flux fraction
        scalar massFluxFraction_;

public:

    //- Runtime type information
    TypeName("totalFlowRateAdvectiveDiffusive");


    // Constructors

        //- Construct from patch and internal field
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        //  onto a new patch
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const totalFlowRateAdvectiveDiffusiveFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const totalFlowRateAdvectiveDiffusiveFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new
                totalFlowRateAdvectiveDiffusiveFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
        (
            const totalFlowRateAdvectiveDiffusiveFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new
                totalFlowRateAdvectiveDiffusiveFvPatchScalarField
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
