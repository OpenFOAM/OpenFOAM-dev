/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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
    Foam::laminarModels::lambdaThixotropic

Description
    Thixotropic viscosity momentum transport model based on the evolution of
    the structural parameter \f$ \lambda \f$:

        \f[
            \frac{D\lambda}{Dt} = a(1 - \lambda)^b - c \lambda \dot{\gamma}^d
        \f]

    The viscosity is then calculated using the expression

        \f[
            \nu = \frac{\nu_{\infty}}{{1 - K \lambda}^2}
        \f]

    Where the parameter K is given by:

        \f[
            K = 1 - \sqrt{\frac{\nu_{\infty}}{\nu_{0}}}
        \f]

    Here:
    \vartable
        \lambda         | structural parameter
        a               | model coefficient
        b               | model coefficient
        c               | model coefficient
        d               | model coefficient
        \dot{\gamma}    | stress rate [1/s]
        \nu_{0}         | limiting viscosity when \f$ \lambda = 1 \f$
        \nu_{\infty}    | limiting viscosity when \f$ \lambda = 0 \f$
    \endvartable

    Reference:
    \verbatim
        Barnes H A, 1997.  Thixotropy - a review.  J. Non-Newtonian Fluid
        Mech 70, pp 1-33
    \endverbatim

SourceFiles
    lambdaThixotropic.C

\*---------------------------------------------------------------------------*/

#ifndef lambdaThixotropic_H
#define lambdaThixotropic_H

#include "laminarModel.H"
#include "linearViscousStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{

/*---------------------------------------------------------------------------*\
                    Class lambdaThixotropic Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class lambdaThixotropic
:
    public linearViscousStress<laminarModel<BasicMomentumTransportModel>>
{
    // Private data

        //- Model a coefficient
        dimensionedScalar a_;

        //- Model b coefficient
        dimensionedScalar b_;

        //- Model d coefficient
        dimensionedScalar d_;

        //- Model c coefficient (read after d since dims depend on d value)
        dimensionedScalar c_;

        //- Limiting viscosity when lambda = 1
        dimensionedScalar nu0_;

        //- Limiting viscosity when lambda = 0
        dimensionedScalar nuInf_;

        //- Model coefficient
        dimensionedScalar K_;

        //- Switch for optional Bingham plastic handling
        //  Set by the presence of the sigmay entry
        bool BinghamPlastic_;

        //- Optional Bingham plastic yield stress [m^2/s^2]
        dimensionedScalar sigmay_;

        //- Structural parameter
        //  0 = freestream value (most liquid)
        //  1 = fully built (most solid)
        volScalarField lambda_;

        //- The non-Newtonian viscosity field
        volScalarField nu_;


    // Private Member Functions

        //- Calculates and returns the viscosity from the current lambda
        tmp<volScalarField> calcNu
        (
            const volScalarField& strainRate
        ) const;

        //- Returns the current strain rate from the velocity field
        tmp<volScalarField> strainRate() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("lambdaThixotropic");


    // Constructors

        //- Construct from components
        lambdaThixotropic
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    //- Destructor
    virtual ~lambdaThixotropic()
    {}


    // Member Functions

        //- Read momentumTransport dictionary
        virtual bool read();

        //- Return the effective viscosity
        //  i.e. the lambdaThixotropic viscosity
        virtual tmp<volScalarField> nuEff() const;

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const;

        //- The lambdaThixotropic viscosity is not predicted
        virtual void predict()
        {}

        //- Correct the lambdaThixotropic viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lambdaThixotropic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
