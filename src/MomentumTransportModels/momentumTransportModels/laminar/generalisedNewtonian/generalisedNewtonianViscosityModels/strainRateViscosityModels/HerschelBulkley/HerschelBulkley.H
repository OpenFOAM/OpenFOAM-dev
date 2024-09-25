/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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
    Foam::laminarModels::generalisedNewtonianViscosityModels::HerschelBulkley

Description
    Herschel-Bulkley generalised Newtonian viscosity model

    The Herschel-Bulkley model combines the effects of a Bingham plastic
    and power-law behaviour of a fluid. At a low strain rate \f$\gamma\f$
    the fluid adopts a viscosity \f$\nu_0\f$. Beyond a threshold stress
    \f$\tau_0\f$, the viscosity is described as a power law.

    Kinematic viscosity [m^2/s]

    \f[
        \nu = \min(\nu_0, \tau_0/\gamma + k\gamma^{n - 1})
    \f]

    Example specification:
    \verbatim
         viscosityModel HerschelBulkley;

         tau0       0.01;
         k          0.001;
         n          0.2;
    \endverbatim

    Note the viscosity \c nu0 at zero strain rate is a physical property,
    generally specified in the physicalProperties file.

SourceFiles
    HerschelBulkley.C

\*---------------------------------------------------------------------------*/

#ifndef HerschelBulkley_H
#define HerschelBulkley_H

#include "strainRateViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{

/*---------------------------------------------------------------------------*\
                       Class HerschelBulkley Declaration
\*---------------------------------------------------------------------------*/

class HerschelBulkley
:
    public strainRateViscosityModel
{
    // Private Data

        dimensionedScalar n_;
        dimensionedScalar k_;
        dimensionedScalar tau0_;


public:

    //- Runtime type information
    TypeName("HerschelBulkley");


    // Constructors

        //- Construct from components
        HerschelBulkley
        (
            const dictionary& viscosityProperties,
            const Foam::viscosity& viscosity,
            const volVectorField& U
        );


    //- Destructor
    virtual ~HerschelBulkley()
    {}


    // Member Functions

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu
        (
            const volScalarField& nu0,
            const volScalarField& strainRate
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace generalisedNewtonianViscosityModels
} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
