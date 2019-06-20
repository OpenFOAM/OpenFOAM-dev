/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::diameterModels::coalescenceModels::LehrMilliesMewesCoalescence

Description
    Model of Lehr et al. (2002). The coalescence rate is calculated by

    \f[
        \frac{\pi}{4} (d_i + d_j)^2 \mathrm{min}(u^{\prime}, u_{\mathrm{crit}})
        \mathrm{exp}
        \left[
          - \left(\frac{\alpha_{\mathrm{max}}^{1/3}}{\alpha^{1/3}} - 1\right)^2
        \right]\,,
    \f]

    where

    \f[
        u^{\prime}
      =
        \mathrm{max}
        \left(
            \sqrt{2} \epsilon_c^{1/3} \sqrt{d_i^{2/3} + d_j^{2/3}},
            |\vec{u}_i - \vec{u}_j|
        \right)
    \f]

    is the characteristic velocity for coalescence. Note that a velocity
    difference between bubble i and j is only present if the corresponding
    size groups are assigned to different velocity groups.

    \vartable
        d_i           |  Diameter of bubble i [m]
        d_j           |  Diameter of bubble j [m]
        \epsilon_c    |  Turbulent dissipation rate of continuous phase [m2/s^3]
        \alpha        |  Total void fraction of the bubbles [-]
        \alpha_{max}  |  Maximum packing density of the bubbles [-]
        u_{crit}      |  Critical velocity for coalescence [m/s]
        \vec{u}_i     |  Velocity vector of bubble i [m/s]
        \vec{u}_j     |  Velocity vector of bubble j [m/s]
    \endvartable

    Reference:
    \verbatim
        Lehr, F., Millies, M., & Mewes, D. (2002).
        Bubble‐size distributions and flow fields in bubble columns.
        AIChE Journal, 48(11), 2426-2443.
    \endverbatim

Usage
    \table
        Property     | Description                | Required    | Default value
        uCrit        | Crit. coalescence velocity | no          | 0.08
        alphaMax     | Maximum packing density    | no          | 0.6
    \endtable

SourceFiles
    LehrMilliesMewesCoalescence.C

\*---------------------------------------------------------------------------*/

#ifndef LehrMilliesMewesCoalescence_H
#define LehrMilliesMewesCoalescence_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                 Class LehrMilliesMewesCoalescence Declaration
\*---------------------------------------------------------------------------*/

class LehrMilliesMewesCoalescence
:
    public coalescenceModel
{
    // Private Data

        //- Critical velocity for coalescence, defaults to 0.08 m/s
        dimensionedScalar uCrit_;

        //- Maximum packing density of the bubbles, defaults to 0.6
        dimensionedScalar alphaMax_;


public:

    //- Runtime type information
    TypeName("LehrMilliesMewes");

    // Constructor

        LehrMilliesMewesCoalescence
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LehrMilliesMewesCoalescence()
    {}


    // Member Functions

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
