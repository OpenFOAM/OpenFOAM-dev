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
    Foam::diameterModels::coalescenceModels::PrinceBlanch

Description
    Model of Prince and Blanch (1990). The coalescence rate is calculated by

    \f[
        \left( \theta_{ij}^{T} + \theta_{ij}^{B} + \theta_{ij}^{LS} \right)
        \lambda_{ij}\;,
    \f]

    with the coalescence efficiency

    \f[
        \lambda_{ij} =
            \mathrm{exp}
            \left(
              - \sqrt{\frac{r_{ij}^3 \rho_c}{16 \sigma}}
                \mathrm{ln} \left(\frac{h_0}{h_f}\right)
                \epsilon_c^{1/3}/r_{ij}^{2/3}
            \right)\;,
    \f]

    the turbulent collision rate

    \f[
        \theta_{ij}^{T} =
            C_1 \pi (d_i + d_j)^{2} \epsilon_c^{1/3}
            \sqrt{d_{i}^{2/3} + d_{j}^{2/3}}\;,
    \f]

    the buoyancy-driven collision rate

    \f[
        \theta_{ij}^{B} = S_{ij} \left| u_{ri} - u_{rj} \right|\;,
    \f]

    and the laminar shear collision rate

    \f[
        \theta_{ij}^{LS} =
            \frac{1}{6} \left(d_i + d_j\right)^{3} \gamma_c\;.
    \f]

    The rise velocity of bubble i is calculated by

    \f[
        u_{ri} = \sqrt{2.14 \sigma / \left(\rho_c d_i \right) + 0.505 g d_i}\;,
    \f]

    the equivalent radius by

    \f[
        r_{ij} = \left( \frac{1}{d_i} + \frac{1}{d_j} \right)^{-1}\;,
    \f]

    the collision cross sectional area by

    \f[
        S_{ij} = \frac{\pi}{4} \left(d_i + d_j\right)^{2}\;,
    \f]

    and the shear strain rate by

    \f[
        \dot{\gamma_{b}} = \mathrm{mag}(\mathrm{symm}(\mathrm{grad}(U_c)))\;.
    \f]

    Note that in equation 2, the bubble radius has been substituted by the
    bubble diameter. Also the expression for the equivalent radius r_ij
    (equation 19 in the paper of Prince and Blanch (1990)) was corrected.

    \vartable
        \theta_{ij}^{T}   |  Turbulent collision rate [m3/s]
        \theta_{ij}^{B}   |  Buoyancy-driven collision rate [m3/s]
        \theta_{ij}^{LS}  |  Laminar shear collision rate [m3/s]
        \lambda_{ij}      |  Coalescence efficiency [-]
        r_{ij}            |  Equivalent radius [m]
        \rho_c            |  Density of continuous phase [kg/m3]
        \sigma            |  Surface tension [N/m]
        h_0               |  Initial film thickness [m]
        h_f               |  Critical film thickness [m]
        \epsilon_c        |  Continuous phase turbulent dissipation rate [m2/s3]
        d_i               |  Diameter of bubble i [m]
        d_j               |  Diameter of bubble j [m]
        u_{ri}            |  Rise velocity of bubble i [m/s]
        S_{ij}            |  Collision cross sectional area [m2]
        g                 |  Gravitational constant [m/s2]
        \gamma_c          |  Continuous phase shear strain rate [1/s]
        U_c               |  Continuous phase velocity field [m/s]
    \endvartable

    Reference:
    \verbatim
        Prince, M. J., & Blanch, H. W. (1990).
        Bubble coalescence and break‐up in air‐sparged bubble columns.
        AIChE journal, 36(10), 1485-1499.
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        C1           | coefficient C1          | no          | 0.089
        h0           | initial film thickness  | no          | 1e-4m
        hf           | critical film thickness | no          | 1e-8m
        turbulence   | Switch for collisions due to turbulence    | yes | none
        buoyancy     | Switch for collisions due to buoyancy      | yes | none
        laminarShear | Switch for collisions due to laminar shear | yes | none
    \endtable

SourceFiles
    PrinceBlanch.C

\*---------------------------------------------------------------------------*/

#ifndef PrinceBlanch_H
#define PrinceBlanch_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                        Class PrinceBlanch Declaration
\*---------------------------------------------------------------------------*/

class PrinceBlanch
:
    public coalescenceModel
{
    // Private Data

        //- Optional coefficient C1, defaults to 0.089
        dimensionedScalar C1_;

        //- Initial film thickness, defaults to 1e-4m
        dimensionedScalar h0_;

        //- Critical film thickness, defaults to 1e-8m
        dimensionedScalar hf_;

        //- Switch for considering turbulent collisions
        Switch turbulence_;

        //- Switch for considering buoyancy-induced collisions
        Switch buoyancy_;

        //- Switch for considering buoyancy-induced collisions
        Switch laminarShear_;

        //- Continuous phase shear strain rate
        autoPtr<volScalarField> shearStrainRate_;


public:

    //- Runtime type information
    TypeName("PrinceBlanch");

    // Constructor

        PrinceBlanch
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~PrinceBlanch()
    {}


    // Member Functions

        //- Correct diameter independent expressions
        virtual void correct();

        //- Add to coalescenceRate
        virtual void addToCoalescenceRate
        (
            volScalarField& coalescenceRate,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
