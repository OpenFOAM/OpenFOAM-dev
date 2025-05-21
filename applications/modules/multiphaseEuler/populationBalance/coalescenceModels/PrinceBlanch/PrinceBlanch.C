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

#include "PrinceBlanch.H"
#include "fvcGrad.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(PrinceBlanch, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        PrinceBlanch,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::PrinceBlanch::
PrinceBlanch
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    C1_("C1", dimless, dict, 0.356),
    h0_("h0", dimLength, dict, 1e-4),
    hf_("hf", dimLength, dict, 1e-8),
    turbulence_(dict.lookup("turbulence")),
    buoyancy_(dict.lookup("buoyancy")),
    laminarShear_(dict.lookup("laminarShear"))
{
    if (laminarShear_)
    {
        shearStrainRate_.set
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    "shearStrainRate",
                    popBal_.time().name(),
                    popBal_.mesh()
                ),
                popBal_.mesh(),
                dimensionedScalar
                (
                    "shearStrainRate",
                    dimVelocity/dimLength,
                    Zero
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::coalescenceModels::PrinceBlanch::precompute()
{
    if (laminarShear_)
    {
        shearStrainRate_() =
            sqrt(2.0)*mag(symm(fvc::grad(popBal_.continuousPhase().U())));
    }
}


void Foam::diameterModels::coalescenceModels::PrinceBlanch::
addToCoalescenceRate
(
    volScalarField::Internal& coalescenceRate,
    const label i,
    const label j
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const sizeGroup& fj = popBal_.sizeGroups()[j];

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(fi.phase()));
    const volScalarField::Internal& sigma = tsigma();

    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volScalarField::Internal& epsilonc = tepsilonc();

    const uniformDimensionedVectorField& g =
        popBal_.mesh().lookupObject<uniformDimensionedVectorField>("g");

    const dimensionedScalar rij(1/(1/fi.dSph() + 1/fj.dSph()));

    const volScalarField::Internal collisionEfficiency
    (
        exp
        (
          - sqrt(pow3(rij)*rhoc/(16*sigma))
           *log(h0_/hf_)
           *cbrt(epsilonc)
           /pow(rij, 2.0/3.0)
        )
    );

    if (turbulence_)
    {
        coalescenceRate +=
            (
                C1_
               *pi
               *sqr(fi.dSph() + fj.dSph())
               *cbrt(epsilonc)
               *sqrt(pow(fi.dSph(), 2.0/3.0) + pow(fj.dSph(), 2.0/3.0))
            )
           *collisionEfficiency;
    }

    if (buoyancy_)
    {
        const dimensionedScalar Sij(pi/4*sqr(fi.dSph() + fj.dSph()));

        coalescenceRate +=
            Sij
           *mag
            (
                sqrt(2.14*sigma/(rhoc*fi.dSph()) + 0.505*mag(g)*fi.dSph())
              - sqrt(2.14*sigma/(rhoc*fj.dSph()) + 0.505*mag(g)*fj.dSph())
            )
           *collisionEfficiency;
    }

    if (laminarShear_)
    {
        coalescenceRate +=
            pow3(fi.dSph() + fj.dSph())/6
           *shearStrainRate_()*collisionEfficiency;
    }
}


// ************************************************************************* //
