/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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
    Foam::surfaceFilmModels::thixotropicViscosity

Description
    Thixotropic viscosity model based on the evolution of the structural
    parameter \f$ \lambda \f$:

        \f[
            \frac{D\lambda}{Dt} = a(1 - \lambda)^b - c \lambda \dot{\gamma}^d
        \f]

    The viscosity is then calculated using the expression

        \f[
            \mu = \frac{\mu_{\infty}}{{1 - K \lambda}^2}
        \f]

    Where the parameter K is given by:

        \f[
            K = 1 - \sqrt{\frac{\mu_{\infty}}{\mu_{0}}}
        \f]

    Here:
    \vartable
        \lambda         | structural parameter
        a               | model coefficient
        b               | model coefficient
        c               | model coefficient
        d               | model coefficient
        \dot{\gamma}    | stress rate [1/s]
        \mu_{0}         | limiting viscosity when \f$ \lambda = 1 \f$
        \mu_{\infty}    | limiting viscosity when \f$ \lambda = 0 \f$
    \endvartable

    Reference:
    \verbatim
        Barnes H A, 1997.  Thixotropy - a review.  J. Non-Newtonian Fluid
        Mech 70, pp 1-33
    \endverbatim

SourceFiles
    thixotropicViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef thixotropicViscosity_H
#define thixotropicViscosity_H

#include "filmViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class thixotropicViscosity Declaration
\*---------------------------------------------------------------------------*/

class thixotropicViscosity
:
    public viscosityModel
{
protected:

    // Protected data

        //- Model `a' coefficient
        dimensionedScalar a_;

        //- Model `b' coefficient
        dimensionedScalar b_;

        //- Model `d' coefficient
        dimensionedScalar d_;

        //- Model `c' coefficient (read after d since dims depend on d value)
        dimensionedScalar c_;

        //- Limiting viscosity when lambda = 1
        dimensionedScalar mu0_;

        //- Limiting viscosity when lambda = 0
        dimensionedScalar muInf_;

        bool BinghamPlastic_;

        dimensionedScalar tauy_;

        //- Model coefficient
        dimensionedScalar K_;

        //- Structural parameter
        //  0 = freestream value (most liquid)
        //  1 = fully built (most solid)
        volScalarField lambda_;


public:

    //- Runtime type information
    TypeName("thixotropic");


    // Constructors

        //- Construct from surface film model
        thixotropicViscosity
        (
            surfaceFilm& film,
            const dictionary& dict,
            volScalarField& mu
        );

        //- Disallow default bitwise copy construction
        thixotropicViscosity(const thixotropicViscosity&) = delete;


    //- Destructor
    virtual ~thixotropicViscosity();


    // Member Functions

        //- Correct
        virtual void correct
        (
            const volScalarField& p,
            const volScalarField& T
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const thixotropicViscosity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
