/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::mixtureViscosityModels::HerschelBulkley

Description
    HerschelBulkley viscosity model

    The Herschel-Bulkley model combines the effects of a Bingham plastic
    and power-law behaviour of a fluid. At a low strain rate \f$\gamma\f$
    the fluid adopts a viscosity \f$\mu_0\f$. Beyond a threshold stress
    \f$\tau_0\f$, the viscosity is described as a power law.

    Kinematic viscosity [m^2/s]

    \f[
        \mu = \min(\mu_0, (\tau_0 + k\gamma^{n})/\gamma)
    \f]

Usage
    Example usage:
    \verbatim
        viscosityModel HerschelBulkley;

        rho        1996; // Dispersed phase density

        tau0       10.0;
        k          0.001;
        n          0.2;
    \endverbatim

SourceFiles
    HerschelBulkley.C

\*---------------------------------------------------------------------------*/

#ifndef HerschelBulkley_H
#define HerschelBulkley_H

#include "mixtureViscosityModel.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class HerschelBulkley Declaration
\*---------------------------------------------------------------------------*/

class HerschelBulkley
:
    public mixtureViscosityModel
{
    // Private data

        dimensionedScalar n_;
        dimensionedScalar k_;
        dimensionedScalar tau0_;


public:

    //- Runtime type information
    TypeName("HerschelBulkley");


    // Constructors

        //- Construct from mixture
        HerschelBulkley(const incompressibleDriftFluxMixture& mixture);


    //- Destructor
    virtual ~HerschelBulkley()
    {}


    // Member Functions

        //- Return the mixture viscosity
        //  given the viscosity of the continuous phase
        virtual tmp<volScalarField> mu
        (
            const volScalarField& muc,
            const volVectorField& U
        ) const;

        //- Read phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixtureViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
