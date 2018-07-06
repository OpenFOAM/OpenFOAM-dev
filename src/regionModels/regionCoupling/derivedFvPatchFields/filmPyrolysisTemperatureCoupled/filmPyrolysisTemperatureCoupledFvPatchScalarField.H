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
    Foam::filmPyrolysisTemperatureCoupledFvPatchScalarField

Description
    This boundary condition is designed to be used in conjunction with surface
    film and pyrolysis modelling.  It provides a temperature boundary condition
    for patches on the primary region based on whether the patch is seen to
    be 'wet', retrieved from the film alpha field.

    - if the patch is wet, the temperature is set using the film temperature
    - otherwise, it is set using pyrolysis temperature

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            filmPyrolysisTemperatureCoupled;
        phi             phi;      // name of flux field (default = phi)
        rho             rho;      // name of density field (default = rho)
        deltaWet        1e-4;     // threshold height for 'wet' film
        value           uniform   300; // initial temperature / [K]
    }
    \endverbatim

SourceFiles
    filmPyrolysisTemperatureCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef filmPyrolysisTemperatureCoupledFvPatchScalarField_H
#define filmPyrolysisTemperatureCoupledFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class filmPyrolysisTemperatureCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class filmPyrolysisTemperatureCoupledFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of film region
        const word filmRegionName_;

        //- Name of pyrolysis region
        const word pyrolysisRegionName_;

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;


public:

    //- Runtime type information
    TypeName("filmPyrolysisTemperatureCoupled");


    // Constructors

        //- Construct from patch and internal field
        filmPyrolysisTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        filmPyrolysisTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  filmPyrolysisTemperatureCoupledFvPatchScalarField onto a new patch
        filmPyrolysisTemperatureCoupledFvPatchScalarField
        (
            const filmPyrolysisTemperatureCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        filmPyrolysisTemperatureCoupledFvPatchScalarField
        (
            const filmPyrolysisTemperatureCoupledFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new filmPyrolysisTemperatureCoupledFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        filmPyrolysisTemperatureCoupledFvPatchScalarField
        (
            const filmPyrolysisTemperatureCoupledFvPatchScalarField&,
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
                new filmPyrolysisTemperatureCoupledFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

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
