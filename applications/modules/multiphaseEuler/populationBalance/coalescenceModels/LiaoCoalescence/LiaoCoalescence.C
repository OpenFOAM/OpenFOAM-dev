/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
#include "phaseSystem.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "fvcGrad.H"
#include "dragModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
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

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::LiaoCoalescence::
LiaoCoalescence
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

void Foam::diameterModels::coalescenceModels::LiaoCoalescence::precompute()
{
    LiaoBase::precompute();

    CPack_ = min(PMax_/max(PMax_ - popBal_.alphas(), SMALL), CPackMax_);

    const uniformDimensionedVectorField& g =
        popBal_.mesh().lookupObject<uniformDimensionedVectorField>("g");

    dCrit_ =
        4*sqrt
        (
            popBal_.sigmaWithContinuousPhase(popBal_.sizeGroups()[1].phase())()
           /(mag(g)*(popBal_.continuousPhase().rho()
          - popBal_.sizeGroups()[1].phase().rho()))
        );
}


void Foam::diameterModels::coalescenceModels::LiaoCoalescence::
addToCoalescenceRate
(
    volScalarField& coalescenceRate,
    const label i,
    const label j
)
{
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const sizeGroup& fj = popBal_.sizeGroups()[j];

    dimensionedScalar dEq(2*fi.dSph()*fj.dSph()/(fi.dSph() + fj.dSph()));
    dimensionedScalar Aij(pi*0.25*sqr(fi.dSph() + fj.dSph()));

    if (turbulence_)
    {
        uRelTurb_ =
            CTurb_*sqrt(2.0)
           *sqrt(sqr(cbrt(fi.dSph())) + sqr(cbrt(fj.dSph())))
           *cbrt(popBal_.continuousTurbulence().epsilon());
    }

    if (buoyancy_)
    {
        uRelBuoy_ = CBuoy_*mag(uTerminal_[i] - uTerminal_[j]);
    }

    if (laminarShear_)
    {
        uRelShear_ = CShear_*0.5/pi*(fi.dSph() + fj.dSph())*shearStrainRate_;
    }

    const volScalarField collisionEfficiency
    (
        neg(kolmogorovLengthScale_ - (fi.dSph() + fj.dSph()))
       *exp
        (
          - CEff_*sqrt
            (
                continuousPhase.rho()*dEq
               /popBal_.sigmaWithContinuousPhase(fi.phase())
               *sqr(max(uRelTurb_, max(uRelBuoy_, uRelShear_)))
            )
        )
      + pos0(kolmogorovLengthScale_ - (fi.dSph() + fj.dSph()))
       *exp
        (
          - 3*continuousPhase.fluidThermo().mu()*dEq*eddyStrainRate_
           /(4*popBal_.sigmaWithContinuousPhase(fi.phase()))
           *log
            (
                cbrt
                (
                    pi*popBal_.sigmaWithContinuousPhase(fi.phase())*sqr(dEq)
                   /(32*AH_)
                )
            )
        )
    );

    if (turbulence_)
    {
        coalescenceRate +=
            neg(kolmogorovLengthScale_ - (fi.dSph() + fj.dSph()))
           *CPack_*Aij*uRelTurb_*collisionEfficiency;
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
        const volScalarField uRelEddy
        (
            CEddy_*0.5/pi*(fi.dSph() + fj.dSph())*eddyStrainRate_
        );

        coalescenceRate +=
            pos0(kolmogorovLengthScale_ - (fi.dSph() + fj.dSph()))
           *CPack_*0.5*Aij*uRelEddy*collisionEfficiency;
    }

    if (wakeEntrainment_)
    {
        const dimensionedScalar uRelWakeI
        (
            CWake_*uTerminal_[i]*cbrt(Cd_[i])
        );

        const dimensionedScalar uRelWakeJ
        (
            CWake_*uTerminal_[j]*cbrt(Cd_[j])
        );

        coalescenceRate +=
            CPack_*0.125*pi
           *(
                sqr(fi.dSph())*uRelWakeI
               *pos0(fi.dSph() - 0.5*dCrit_)
               *(
                    pow6(fi.dSph() - 0.5*dCrit_)
                   /(pow6(fi.dSph() - 0.5*dCrit_) + pow6(0.5*dCrit_))
                )
              + sqr(fj.dSph())*uRelWakeJ
               *pos0(fj.dSph() - 0.5*dCrit_)
               *(
                    pow6(fj.dSph() - 0.5*dCrit_)
                   /(pow6(fj.dSph() - 0.5*dCrit_) + pow6(0.5*dCrit_))
                )
            );
    }
}


// ************************************************************************* //
