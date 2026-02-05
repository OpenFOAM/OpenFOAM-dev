/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "phaseSurfaceCondensation.H"
#include "multiphaseEuler.H"
#include "fluidMulticomponentThermo.H"
#include "fluidThermophysicalTransportModel.H"
#include "diameterModel.H"
#include "twoResistanceHeatTransfer.H"
#include "generateBlendedInterfacialModels.H"
#include "saturationPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(phaseSurfaceCondensation, 0);
    addToRunTimeSelectionTable(fvModel, phaseSurfaceCondensation, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::phaseSurfaceCondensation::readCoeffs(const dictionary& dict)
{
    reReadSpecie(dict);

    const dictionary& diffusiveMassTransferDict =
        dict.subDict("diffusiveMassTransfer");

    checkBlendedInterfacialModelsDict<blendedSidedDiffusiveMassTransferModel>
    (
        fluid_,
        diffusiveMassTransferDict
    );

    const phaseInterface interface(vapour_, solid_);

    diffusiveMassTransferModel_.reset
    (
        blendedSidedDiffusiveMassTransferModel::New
        (
            diffusiveMassTransferDict,
            interface,
            blendingDict<blendedSidedDiffusiveMassTransferModel>
            (
                fluid_,
                diffusiveMassTransferDict
            )
        ).ptr()
    );

    saturationModelPtr_.reset
    (
        saturationPressureModel::New
        (
            "saturationPressure",
            dict
        ).ptr()
    );

    specieSemiImplicit_ =
        dict.lookupOrDefault<bool>("specieSemiImplicit", false);
}


void Foam::fv::phaseSurfaceCondensation::correctMDot() const
{
    Info<< type() << ": " << name() << endl << incrIndent;

    static const dimensionedScalar rootVSmallH
    (
        heatTransferModel::dimK,
        rootVSmall
    );

    const rhoThermo& solidThermo = solid_.thermo();
    const volScalarField::Internal& solidT = solidThermo.T();

    const fluidMulticomponentThermo& vapourThermo =
        fluidMulticomponentThermos(true, false)[0];
    const volScalarField::Internal& vapourT = vapourThermo.T();

    const label speciei = specieis()[0];

    const Pair<tmp<volScalarField>> Hs =
        solver_.heatTransfer.Hs(vapour_, solid_, scalar(0));
    const volScalarField::Internal& vapourH = Hs.first();
    const volScalarField::Internal& solidH = Hs.second();
    const volScalarField::Internal Tsurface
    (
        (solidH*solidT + vapourH*vapourT + q_)
       /max(solidH + vapourH, rootVSmallH)
    );

    const volScalarField::Internal pSat
    (
        saturationModelPtr_->pSat(Tsurface)
    );

    const volScalarField::Internal L(this->L(Tsurface));

    const fluidThermophysicalTransportModel& ttmVapour =
        mesh().lookupType<fluidThermophysicalTransportModel>
        (
            vapour_.name()
        );

    const volScalarField::Internal freeSurf(vapour_/(1 - solid_));

    const volScalarField::Internal xc
    (
        vapour_.Y(species()[0])/vapourThermo.Wi(speciei)*vapourThermo.W()
    );

    const volScalarField::Internal xw(pSat/vapourThermo.p()());

    // Optional relaxation factor
    const scalar f = mesh().solution().fieldRelaxationFactor(mDot_.member());

    mDot_ =
        (1 - f)*mDot_
      - f*freeSurf*vapourThermo.Wi(speciei)/vapourThermo.W()()
       *diffusiveMassTransferModel_->KinThe(vapour_)()
       *ttmVapour.D(vapour_.Y(species()[0]))()
       *log(max(1 - xc, 0.001)/max(1 - xw, 0.001));
    mDot_.max(0);

    infoField("mDot", mDot_);

    if (specieSemiImplicit_)
    {
        mDotDy_ =
            (1 - f)*mDotDy_
          + f*freeSurf
           *diffusiveMassTransferModel_->KinThe(vapour_)()
           *ttmVapour.D(vapour_.Y(species()[0]))()
           /(1 - min(xc, 0.999))*pos(xc - xw);
    }
    else
    {
        mDotDy_ = Zero;
    }

    // Heat flux
    q_ = mDot_*L;

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::phaseSurfaceCondensation::phaseSurfaceCondensation
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseChange
    (
        name,
        modelType,
        mesh,
        dict,
        readSpecie(coeffs(modelType, dict), false)
    ),
    solver_(mesh().lookupObject<solvers::multiphaseEuler>(solver::typeName)),
    fluid_(solver_.fluid),
    liquid_(fluid_.phases()[phaseNames().second()]),
    vapour_(fluid_.phases()[phaseNames().first()]),
    solid_(fluid_.phases()[dict.lookup("phase")]),
    diffusiveMassTransferModel_(nullptr),
    saturationModelPtr_(nullptr),
    pressureEquationIndex_(-1),
    specieSemiImplicit_(false),
    q_
    (
        IOobject
        (
            name + ":q",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, scalar(0))
    ),
    mDot_
    (
        IOobject
        (
            name + ":mDot",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, scalar(0))
    ),
    mDotDy_
    (
        IOobject
        (
            name + ":mDotDy",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, scalar(0))
    )
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::fv::phaseSurfaceCondensation::addsSupToField(const word& fieldName) const
{
    return
        phaseChange::addsSupToField(fieldName)
     || fieldName == solid_.thermo().he().name();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseSurfaceCondensation::Lfraction() const
{
    // Put all the latent heat into the vapour, additional source term
    // will transfer this to solid phase in addSup
    return
        volScalarField::Internal::New
        (
            name() + ":Lfraction",
            mesh(),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseSurfaceCondensation::mDot() const
{
    return mDot_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseSurfaceCondensation::mDotDy() const
{
    return mDotDy_;
}


void Foam::fv::phaseSurfaceCondensation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    // Pressure equation (i.e., continuity, linearised in the pressure)
    if
    (
        (&alpha == &liquid_ || &alpha == &vapour_)
     && (&rho == &liquid_.rho() || &rho == &vapour_.rho())
     && &eqn.psi() == &solver_.p_rgh
    )
    {
        // Ensure that the source is up-to date if this is the first call in
        // the current phase loop
        if (pressureEquationIndex_ % 2 == 0) correctMDot();
        pressureEquationIndex_ ++;
    }

    // Let the base class add the actual source
    massTransfer::addSup(alpha, rho, eqn);
}


void Foam::fv::phaseSurfaceCondensation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& heOrYi,
    fvMatrix<scalar>& eqn
) const
{
    // Vapour or solid energy equation. Apply the additional explicit latent
    // heat transfers.
    if (&heOrYi == &vapour_.thermo().he() || &heOrYi == &solid_.thermo().he())
    {
        const scalar sign = &heOrYi == &solid_.thermo().he() ? -1 : +1;

        eqn += sign*q_;
    }

    // For solid nothing more is needed
    if (&heOrYi == &solid_.thermo().he())
    {
        return;
    }

    const label i = this->index(phaseNames(), alpha.group());

    // If implicit treatment is not needed or this is liquid, use normal addSup
    if (!specieSemiImplicit_ || i != 0)
    {
        return phaseChange::addSup(alpha, rho, heOrYi, eqn);
    }

    const label s = this->sign(phaseNames(), alpha.group());

    const fluidMulticomponentThermo& thermo =
        fluidMulticomponentThermos(true, false)[0];

    const word specieName = heOrYi.member();

    // Mass fraction equation
    if (thermo.containsSpecie(specieName))
    {
        // A non-transferring specie. Do not add a source.
        if (!species().found(specieName)) return;

        // The transferring specie. Add a linearised source.
        tmp<volScalarField::Internal> tmDot = this->mDot();
        tmp<volScalarField::Internal> tmDotDy = this->mDotDy();

        eqn += s*(tmDot() + correction(fvm::Sp(tmDotDy, eqn.psi())));

        return;
    }

    // Something else. Fall back.
    phaseChange::addSup(alpha, rho, heOrYi, eqn);
}


void Foam::fv::phaseSurfaceCondensation::correct()
{
    // Reset the p_rgh equation solution counter
    pressureEquationIndex_ = 0;

    // Correct the total phase change rate
    correctMDot();
}


bool Foam::fv::phaseSurfaceCondensation::read(const dictionary& dict)
{
    if (phaseChange::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
