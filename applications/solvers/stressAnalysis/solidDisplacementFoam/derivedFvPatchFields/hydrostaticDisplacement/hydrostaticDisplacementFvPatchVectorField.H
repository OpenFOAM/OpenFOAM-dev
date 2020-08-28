/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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
    Foam::hydrostaticDisplacementFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation in which the traction is caused by
    the hydrostatic pressure of an external liquid reservoir.

SourceFiles
    hydrostaticDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef hydrostaticDisplacementFvPatchVectorField_H
#define hydrostaticDisplacementFvPatchVectorField_H

#include "tractionDisplacementFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class hydrostaticDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class hydrostaticDisplacementFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{
    // Private Data

        //- Density of surrounding liquid
        scalar rhoLiquid_;

        //- Pressure above the liquid
        scalar liquidSurfacePressure_;

        //- Point on liquid surface
        vector liquidSurfacePoint_;


public:

    //- Runtime type information
    TypeName("hydrostaticDisplacement");


    // Constructors

        //- Construct from patch and internal field
        hydrostaticDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        hydrostaticDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  hydrostaticDisplacementFvPatchVectorField onto a new patch
        hydrostaticDisplacementFvPatchVectorField
        (
            const hydrostaticDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        hydrostaticDisplacementFvPatchVectorField
        (
            const hydrostaticDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new hydrostaticDisplacementFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        hydrostaticDisplacementFvPatchVectorField
        (
            const hydrostaticDisplacementFvPatchVectorField&,
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
                new hydrostaticDisplacementFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

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
