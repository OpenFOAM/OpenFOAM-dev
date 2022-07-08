/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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
    Foam::laminarModels::generalisedNewtonianViscosityModels::BirdCarreau

Description
    Bird-Carreau generalised Newtonian viscosity model

    The Bird-Carreau-Yasuda form is also supported if the optional \c a
    coefficient is specified.  \c a defaults to 2 for the Bird-Carreau model.

    The coefficient applied to strain rate \f$\gamma\f$ can be specified either
    as the constant \c k or the critical stress level at the transition to
    shear thinning \c tauStar if \c tauStar is provided:

    Kinematic viscosity [m^2/s]

    \f[
        \nu = \nu_\infty\,
            + (\nu_0 - \nu_\infty)\,
              \left(1 + (k\gamma)^a \right)^{(n - 1)/a}
    \f]

    or

    \f[
        \nu = \nu_\infty
            + (\nu_0 - \nu_\infty)
              \left(1 + (\frac{\nu_0\gamma}{\tau^*} )^a \right)^{(n - 1)/a}
    \f]

    Example specification for a polymer:
    \verbatim
         viscosityModel BirdCarreau;

         nuInf       0;
         tauStar     90;
         n           0.5;
    \endverbatim

    Note the viscosity \c nu0 at zero strain rate is a physical property,
    generally specified in the physicalProperties file.

SourceFiles
    BirdCarreau.C

\*---------------------------------------------------------------------------*/

#ifndef BirdCarreau_H
#define BirdCarreau_H

#include "strainRateViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{

/*---------------------------------------------------------------------------*\
                         Class BirdCarreau Declaration
\*---------------------------------------------------------------------------*/

class BirdCarreau
:
    public strainRateViscosityModel
{
    // Private Data

        dimensionedScalar nuInf_;
        dimensionedScalar k_;
        dimensionedScalar tauStar_;
        dimensionedScalar n_;
        dimensionedScalar a_;


public:

    //- Runtime type information
    TypeName("BirdCarreau");


    // Constructors

        //- Construct from components
        BirdCarreau
        (
            const dictionary& viscosityProperties,
            const Foam::viscosity& viscosity,
            const volVectorField& U
        );


    //- Destructor
    virtual ~BirdCarreau()
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
