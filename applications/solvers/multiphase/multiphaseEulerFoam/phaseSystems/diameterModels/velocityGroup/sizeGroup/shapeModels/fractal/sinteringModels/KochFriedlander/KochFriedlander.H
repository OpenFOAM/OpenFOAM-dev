/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2021 OpenFOAM Foundation
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
    Foam::sinteringModels::KochFriedlander

Description
    Sintering model of Koch and Friedlander (1990). The characteristic time for
    sintering is given by

    \f[
        \tau = c_s d_{p_i}^n T^m \exp(T_a/T \cdot [1 - d_{p,min}/d_{p_i}])\;.
    \f]

    Note that the correction factor in the exponential function can be
    eliminated by setting \f$d_{p,min}\f$ to zero which is done by default.

    Reference:
    \verbatim
        Koch, W., & Friedlander, S. K. (1990).
        The effect of particle coalescence on the surface area of a coagulating
        aerosol.
        Journal of Colloid and Interface Science, 140(2), 419-427.
    \endverbatim

Usage
    \table
        Property | Description                | Required    | Default value
        Cs       | Sintering time coefficient | yes         | none
        n        | Particle diameter exponent | yes         | none
        m        | Temperature exponent       | yes         | none
        Ta       | Activation temperature     | yes         | none
        dpMin    | Minimum primary particle diameter | no   | 0
    \endtable

SourceFiles
    KochFriedlander.C

\*---------------------------------------------------------------------------*/

#ifndef KochFriedlander_H
#define KochFriedlander_H

#include "sinteringModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace shapeModels
{
namespace sinteringModels
{

/*---------------------------------------------------------------------------*\
                       Class KochFriedlander Declaration
\*---------------------------------------------------------------------------*/

class KochFriedlander
:
    public sinteringModel
{
    // Private Data

        //- Subdictionary containing model coefficients
        const dictionary& dict_;

        //- Sintering time coefficient
        scalar Cs_;

        //- Sintering time primary particle diameter exponent
        scalar n_;

        //- Sintering time temperature exponent
        scalar m_;

        //- Activation temperature
        scalar Ta_;

        //- Minimum primary particle diameter, defaults to Zero
        scalar dpMin_;


public:

    //- Runtime type information
    TypeName("KochFriedlander");


    // Constructors

        //- Construct from a dictionary
        KochFriedlander
        (
            const dictionary& dict,
            const fractal& fractalShape
        );


    //- Destructor
    virtual ~KochFriedlander();


    // Member Functions

        //- Characteristic time for sintering
        virtual tmp<volScalarField::Internal> tau() const;

        //- Sintering source term
        virtual tmp<fvScalarMatrix> R() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sinteringModels
} // End namespace shapeModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
