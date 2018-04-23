/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Model of Laakkonen et al. (2006). The total breakup rate is calculated by

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
        v_i         |  Volume of mother bubble i [m]
        \epsilon_c  |  Turbulent dissipation rate of continuous phase [m2/s3]
        \mu_c       |  Molecular dynamic viscosity of liquid phase [Pa s]
        \rho_c      |  Density of continous phase [kg/m3]
        \rho_d      |  Density of disperse phase [kg/m3]
    \endvartable

    References:
    \verbatim
        "Validation of bubble breakage, coalescence and mass transfer models for
         gas-liquid dispersion in agitated vessel"
        Laakkonen, M., Alopaeus, V., Aittamaa, J.
        Chemical Engineering Science, Vol. 61, 2006, pp. 218-228
        Eq. 2-3, p. 220
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        C1           | Coefficient C1          | no          | 6.0
        C2           | Coefficient C2          | no          | 0.04
        C3           | Coefficient C3          | no          | 0.01
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
    // Private data

        //- Optional coefficient C1, defaults to 6.0
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
