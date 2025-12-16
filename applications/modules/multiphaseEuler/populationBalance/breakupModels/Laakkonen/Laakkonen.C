/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "Laakkonen.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace breakupModels
{
    defineTypeNameAndDebug(Laakkonen, 0);
    addToRunTimeSelectionTable(breakupModel, Laakkonen, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::breakupModels::Laakkonen::Laakkonen
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    daughterSizeDistribution(popBal, dict),
    C1_("C1", dimensionSet(0, -2.0/3.0, 0, 0, 0), dict, 2.25),
    C2_("C2", dimless, dict, 0.04),
    C3_("C3", dimless, dict, 0.01)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::populationBalance::breakupModels::Laakkonen::rate(const label i) const
{
    const dimensionedScalar& dSphi = popBal_.dSph(i);

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();
    const volScalarField::Internal& rhod = popBal_.phases()[i].rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(i));
    const volScalarField::Internal& sigma = tsigma();

    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volScalarField::Internal& epsilonc = tepsilonc();
    tmp<volScalarField> tmu(popBal_.continuousPhase().fluidThermo().mu());
    const volScalarField::Internal muc = tmu();

    return
        C1_
       *cbrt(epsilonc)
       *erfc
        (
            sqrt
            (
                C2_
               *sigma
               /(
                   rhoc*pow(dSphi, 5.0/3.0)
                  *pow(epsilonc, 2.0/3.0)
                )
              + C3_
               *muc
               /(
                    sqrt(rhoc*rhod)
                   *cbrt(epsilonc)
                   *pow(dSphi, 4.0/3.0)
                )
            )
        );
}


// ************************************************************************* //
