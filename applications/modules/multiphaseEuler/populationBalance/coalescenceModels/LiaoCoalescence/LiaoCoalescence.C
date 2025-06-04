/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

#include "LiaoCoalescence.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(LiaoCoalescence, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        LiaoCoalescence,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::coalescenceModels::LiaoCoalescence::LiaoCoalescence
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    LiaoBase(popBal, dict),
    PMax_("PMax", dimless, dict, 0.8),
    AH_("AH", dimEnergy, dict, 3.7e-20),
    CEff_("CEff", dimless, dict, 2.5),
    CTurb_("CTurb", dimless, dict, 1),
    CBuoy_("CBuoy", dimless, dict, 1),
    CShear_("CShear", dimless, dict, 1),
    CEddy_("CEddy", dimless, dict, 1),
    CWake_("CWake", dimless, dict, 1),
    turbulence_(dict.lookup("turbulence")),
    buoyancy_(dict.lookup("buoyancy")),
    laminarShear_(dict.lookup("laminarShear")),
    eddyCapture_(dict.lookup("eddyCapture")),
    wakeEntrainment_(dict.lookup("wakeEntrainment")),
    CPack_
    (
        IOobject
        (
            "CPack",
            popBal_.time().name(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar
        (
            "CPack",
            dimless,
            Zero
        )
    ),
    CPackMax_("CPackMax", dimless, dict, 1e5),
    dCrit_
    (
        IOobject
        (
            "dCrit",
            popBal_.time().name(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar
        (
            "dCrit",
            dimLength,
            Zero
        )
    ),
    uRelTurb_
    (
        IOobject
        (
            "uRelTurb",
            popBal_.time().name(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar
        (
            "uRelTurb",
            dimVelocity,
            Zero
        )
    ),
    uRelBuoy_
    (
        IOobject
        (
            "uRelBuoy",
            popBal_.time().name(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar
        (
            "uRelBuoy",
            dimVelocity,
            Zero
        )
    ),
    uRelShear_
    (
        IOobject
        (
            "uRelShear",
            popBal_.time().name(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar
        (
            "uRelShear",
            dimVelocity,
            Zero
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalance::coalescenceModels::LiaoCoalescence::precompute()
{
    LiaoBase::precompute();

    CPack_ = min(PMax_/max(PMax_ - popBal_.alphas(), small), CPackMax_);

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(0));
    const volScalarField::Internal& sigma = tsigma();

    const uniformDimensionedVectorField& g =
        popBal_.mesh().lookupObject<uniformDimensionedVectorField>("g");

    dCrit_ = 4*sqrt(sigma/(mag(g)*(rhoc - popBal_.phases().first().rho()())));
}


void Foam::populationBalance::coalescenceModels::LiaoCoalescence::
addToCoalescenceRate
(
    volScalarField::Internal& coalescenceRate,
    const label i,
    const label j
)
{
    using Foam::constant::mathematical::pi;

    const dimensionedScalar& dSphi = popBal_.dSph(i);
    const dimensionedScalar& dSphj = popBal_.dSph(j);

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(i));
    const volScalarField::Internal& sigma = tsigma();

    tmp<volScalarField> tmuc(popBal_.continuousPhase().fluidThermo().mu());
    const volScalarField::Internal& muc = tmuc();

    dimensionedScalar dEq(2*dSphi*dSphj/(dSphi + dSphj));
    dimensionedScalar Aij(pi*0.25*sqr(dSphi + dSphj));

    if (turbulence_)
    {
        tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
        const volScalarField::Internal& epsilonc = tepsilonc();

        uRelTurb_ =
            CTurb_*sqrt(2.0)
           *sqrt(sqr(cbrt(dSphi)) + sqr(cbrt(dSphj)))
           *cbrt(epsilonc);
    }

    if (buoyancy_)
    {
        uRelBuoy_ = CBuoy_*mag(uTerminal_[i] - uTerminal_[j]);
    }

    if (laminarShear_)
    {
        uRelShear_ = CShear_*0.5/pi*(dSphi + dSphj)*shearStrainRate_;
    }

    const volScalarField::Internal collisionEfficiency
    (
        neg(kolmogorovLengthScale_ - (dSphi + dSphj))
       *exp
        (
          - CEff_
           *sqrt
            (
                rhoc
               *dEq
               /sigma
               *sqr(max(uRelTurb_, max(uRelBuoy_, uRelShear_)))
            )
        )
      + pos0(kolmogorovLengthScale_ - (dSphi + dSphj))
       *exp
        (
          - (3.0/4.0)
           *muc
           *dEq
           *eddyStrainRate_
           /sigma
           *log(cbrt(pi*sigma*sqr(dEq)/(32*AH_)))
        )
    );

    if (turbulence_)
    {
        coalescenceRate +=
            neg(kolmogorovLengthScale_ - (dSphi + dSphj))
           *CPack_
           *Aij
           *uRelTurb_
           *collisionEfficiency;
    }

    if (buoyancy_)
    {
        coalescenceRate += CPack_*0.5*Aij*uRelBuoy_*collisionEfficiency;
    }

    if (laminarShear_)
    {
        coalescenceRate += CPack_*0.5*Aij*uRelShear_*collisionEfficiency;
    }

    if (eddyCapture_)
    {
        const volScalarField::Internal uRelEddy
        (
            CEddy_*0.5/pi*(dSphi + dSphj)*eddyStrainRate_
        );

        coalescenceRate +=
            pos0(kolmogorovLengthScale_ - (dSphi + dSphj))
           *CPack_
           *0.5
           *Aij
           *uRelEddy
           *collisionEfficiency;
    }

    if (wakeEntrainment_)
    {
        const dimensionedScalar uRelWakeI(CWake_*uTerminal_[i]*cbrt(Cd_[i]));

        const dimensionedScalar uRelWakeJ(CWake_*uTerminal_[j]*cbrt(Cd_[j]));

        coalescenceRate +=
            CPack_
           *0.125
           *pi
           *(
                sqr(dSphi)
               *uRelWakeI
               *pos0(dSphi - 0.5*dCrit_)
               *(
                    pow6(dSphi - 0.5*dCrit_)
                   /(pow6(dSphi - 0.5*dCrit_) + pow6(0.5*dCrit_))
                )
              + sqr(dSphj)
               *uRelWakeJ
               *pos0(dSphj - 0.5*dCrit_)
               *(
                    pow6(dSphj - 0.5*dCrit_)
                   /(pow6(dSphj - 0.5*dCrit_) + pow6(0.5*dCrit_))
                )
            );
    }
}


// ************************************************************************* //
