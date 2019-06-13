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
    Foam::filmPyrolysisVelocityCoupledFvPatchVectorField

Description
    This boundary condition is designed to be used in conjunction with surface
    film and pyrolysis modelling.

    It provides a velocity boundary condition for patches on the primary region
    based on whether the patch is seen to be 'wet', retrieved from the film
    alpha field.
      - if the patch is wet, the velocity is set using the film velocity
      - otherwise, it is set using pyrolysis out-gassing velocity

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            filmPyrolysisVelocityCoupled;
        phi             phi;      // name of flux field (default = phi)
        rho             rho;      // name of density field (default = rho)
        deltaWet        1e-4;     // threshold height for 'wet' film
        value           uniform   (0 0 0); // initial velocity / [m/s]
    }
    \endverbatim

SourceFiles
    filmPyrolysisVelocityCoupledFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef filmPyrolysisVelocityCoupledFvPatchVectorField_H
#define filmPyrolysisVelocityCoupledFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class filmPyrolysisVelocityCoupledFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class filmPyrolysisVelocityCoupledFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Name of film region
        word filmRegionName_;

        //- Name of pyrolysis region
        word pyrolysisRegionName_;

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;


public:

    //- Runtime type information
    TypeName("filmPyrolysisVelocityCoupled");


    // Constructors

        //- Construct from patch and internal field
        filmPyrolysisVelocityCoupledFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        filmPyrolysisVelocityCoupledFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  filmPyrolysisVelocityCoupledFvPatchVectorField onto a new patch
        filmPyrolysisVelocityCoupledFvPatchVectorField
        (
            const filmPyrolysisVelocityCoupledFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        filmPyrolysisVelocityCoupledFvPatchVectorField
        (
            const filmPyrolysisVelocityCoupledFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new filmPyrolysisVelocityCoupledFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        filmPyrolysisVelocityCoupledFvPatchVectorField
        (
            const filmPyrolysisVelocityCoupledFvPatchVectorField&,
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
                new filmPyrolysisVelocityCoupledFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }

            //- Return the name of rho
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of rho to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }


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
