/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::diameterModels::binaryBreakupModels::LehrMilliesMewes

Description
    Model of Lehr et al. (2002). The breakup rate is calculated by

    \f[
        0.5 d_j^{*^{5/3}} \mathrm{exp}\left(-\frac{\sqrt{2}}{d_j^{*^{3}}}\right)
        \frac{6}{\pi^{3/2}d_i^{*^{3}}}
        \mathrm{exp}
        \left(
          - \frac{9}{4}\left[\mathrm{ln}\left(2^{2/5} d_i^{*}\right)\right]^{2}
        \right)
        \left(
            1
          + \mathrm{erf}
            \left[
                \frac{3}{2}\mathrm{ln} \left(2^{1/15} d_j^{*}\right)
            \right]
        \right)^{-1}
        \frac{1}{L^{3}T}
    \f]

    with the time scale

    \f[
        T = \left(\frac{\sigma}{\rho_c}\right)^{2/5} \frac{1}{\epsilon_c^{3/5}}
    \f]

    the dimensionless diameter

    \f[
        d^{*} = \frac{d}{L}
    \f]

    and the length scale

    \f[
        L = \left(\frac{\sigma}{\rho_c}\right)^{3/5} \frac{1}{\epsilon_c^{2/5}}
    \f]

    \vartable
        \rho_c          |  Density of continuous phase [kg/m^3]
        \sigma          |  Surface tension [N/m]
        \epsilon_c      |  Continuous phase turbulent dissipation rate [m^2/s^3]
        d_i             |  Diameter of daughter bubble i [m]
        d_j             |  Diameter of mother bubble j [m]
    \endvartable

    References:
    \verbatim
        Lehr, F., Millies, M., & Mewes, D. (2002).
        Bubble‐size distributions and flow fields in bubble columns.
        AIChE Journal, 48(11), 2426-2443.
        Eq. 12-16, p. 2429-2430.
    \endverbatim

SourceFiles
    LehrMilliesMewes.C

\*---------------------------------------------------------------------------*/

#ifndef LehrMilliesMewes_H
#define LehrMilliesMewes_H

#include "binaryBreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace binaryBreakupModels
{

/*---------------------------------------------------------------------------*\
                      Class LehrMilliesMewes Declaration
\*---------------------------------------------------------------------------*/

class LehrMilliesMewes
:
    public binaryBreakupModel
{
public:

    //- Runtime type information
    TypeName("LehrMilliesMewes");

    // Constructor

        LehrMilliesMewes
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LehrMilliesMewes()
    {}


    // Member Functions

        //- Add to binary breakupRate
        virtual void addToBinaryBreakupRate
        (
            volScalarField& binaryBreakupRate,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binaryBreakupModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
