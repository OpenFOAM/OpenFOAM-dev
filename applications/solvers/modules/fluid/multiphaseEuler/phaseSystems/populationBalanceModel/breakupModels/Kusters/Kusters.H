/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::diameterModels::breakupModels::Kusters

Description
    Solid particle breakage model of Kusters (1991). The breakage rate is
    calculated by

    \f[
        \sqrt{\frac{4}{15\pi}}\sqrt{\frac{\epsilon}{\nu}}
        \exp\left(-\frac{\epsilon_{crit}}{\epsilon}\right)
    \f]

    \f[
        \epsilon_{crit}=\frac{B}{R_{ci}}
    \f]

    \f[
        R_{ci} = 0.5d_{p}\left(\frac{n_i}{k_c}\right)^{1/D_f}
    \f]

    \f[
        n_i = \left(\frac{d_i}{d_{p}}\right)^{D_f}
    \f]

    where

    \vartable
        B               |   Critical force parameter for breakage [m3/s3]
        d_i             |   Diameter of transported particle [m]
        d_{p}           |   Diameter of primary particles [m]
        D_f             |   Fractal dimension of particle [-]
        k_c             |   Constant relative to packing density [-]
        n_i             |   Number of primary particles in agglomerate [-]
        \nu             |   Kinematic viscosity of continuous phase [m2/s]
        \epsilon        |   Continuous phase turbulent dissipation rate [m2/s3]
        \epsilon_{crit} |   Critical turbulent dissipation rate [m2/s3]
    \endvartable

    Reference:
    \verbatim
        Kusters, K. A. (1991).
        The influence of turbulence on aggregation of small particles in
        agitated vessels.
        PhD Thesis
    \endverbatim

SourceFiles
    Kusters.C

\*---------------------------------------------------------------------------*/

#ifndef Kusters_H
#define Kusters_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                           Class Kusters Declaration
\*---------------------------------------------------------------------------*/

class Kusters
:
    public breakupModel
{
    // Private Data

        //- Model coefficient
        dimensionedScalar  B_;

        //- Primary particle diameter
        dimensionedScalar dP_;

        //- Packing density coefficient
        dimensionedScalar kc_;

        //- Fractal dimension
        dimensionedScalar Df_;


public:

    //- Runtime type information
    TypeName("Kusters");

    // Constructor

        Kusters
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Kusters()
    {}


    // Member Functions

        //- Set total breakupRate
        virtual void setBreakupRate
        (
            volScalarField& breakupRate,
            const label i
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
