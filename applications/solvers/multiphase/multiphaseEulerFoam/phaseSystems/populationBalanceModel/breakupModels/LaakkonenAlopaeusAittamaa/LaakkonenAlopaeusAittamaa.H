/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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
    Foam::diameterModels::breakupModels::LaakkonenAlopaeusAittamaa

Description
    Model of Laakkonen et al. (2007). The total breakup rate is calculated by

    \f[
        C_1 \epsilon_c^{1/3}
        \mathrm{erfc}
        \left(
            \sqrt
            {
                C_2 \frac{\sigma}{\rho_c \epsilon_c^{2/3} (6v_i/\pi)^{5/9}}
              + C_3 \frac{\mu_c}{\sqrt{\rho_c \rho_d}
                \epsilon_c^{1/3} (6v_i/\pi)^{4/9}}
            }
        \right)
    \f]

    where

    \vartable
        \sigma      |  Surface tension [N/m]
        v_i         |  Volume of mother bubble i [m3]
        \epsilon_c  |  Turbulent dissipation rate of continuous phase [m^2/s^3]
        \mu_c       |  Molecular dynamic viscosity of liquid phase [Pa s]
        \rho_c      |  Density of continuous phase [kg/m^3]
        \rho_d      |  Density of disperse phase [kg/m^3]
    \endvartable

    References:
    \verbatim
        Laakkonen, M., Moilanen, P., Alopaeus, V., & Aittamaa, J. (2007).
        Modelling local bubble size distributions in agitated vessels.
        Chemical Engineering Science, 62, 721–740.
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        C1           | coefficient C1          | no          | 2.25
        C2           | coefficient C2          | no          | 0.04
        C3           | coefficient C3          | no          | 0.01
        daughterSizeDistributionModel | inh. from breakupModel | inherited |
    \endtable

SourceFiles
    LaakkonenAlopaeusAittamaa.C

\*---------------------------------------------------------------------------*/

#ifndef LaakkonenAlopaeusAittamaa_H
#define LaakkonenAlopaeusAittamaa_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                  Class LaakkonenAlopaeusAittamaa Declaration
\*---------------------------------------------------------------------------*/

class LaakkonenAlopaeusAittamaa
:
    public breakupModel
{
    // Private Data

        //- Optional coefficient C1, defaults to 2.25
        dimensionedScalar C1_;

        //- Optional coefficient C2, defaults to 0.04
        dimensionedScalar C2_;

        //- Optional coefficient C3, defaults to 0.01
        dimensionedScalar C3_;


public:

    //- Runtime type information
    TypeName("LaakkonenAlopaeusAittamaa");

    // Constructor

        LaakkonenAlopaeusAittamaa
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LaakkonenAlopaeusAittamaa()
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
