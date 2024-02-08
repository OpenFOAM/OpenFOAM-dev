/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "homogeneousLiquidPhaseSeparation.H"
#include "fundamentalConstants.H"
#include "fluidMulticomponentThermo.H"
#include "physicoChemicalConstants.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(homogeneousLiquidPhaseSeparation, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        homogeneousLiquidPhaseSeparation,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::homogeneousLiquidPhaseSeparation::readCoeffs()
{
    solubilityCurve_.reset
    (
        Function1<scalar>::New("solubility", coeffs()).ptr()
    );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousLiquidPhaseSeparation::YSat
(
    const volScalarField::Internal& T
) const
{
    return
        volScalarField::Internal::New
        (
            typedName(IOobject::groupName("YSat", interface_.name())),
            mesh(),
            dimless,
            solubilityCurve_->value(T)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::homogeneousLiquidPhaseSeparation::homogeneousLiquidPhaseSeparation
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    singleComponentPhaseChangeBase
    (
        name,
        modelType,
        mesh,
        dict,
        {true, false},
        {true, false}
    ),
    fluid_
    (
        mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName)
    ),
    interface_
    (
        fluid_.phases()[phaseNames().first()],
        fluid_.phases()[phaseNames().second()]
    ),
    d_
    (
        IOobject
        (
            typedName(IOobject::groupName("d", interface_.name())),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimLength, 0)
    ),
    mDotByAlphaSolution_
    (
        IOobject
        (
            typedName(IOobject::groupName("mDotByAlpha", interface_.name())),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, 0)
    ),
    solubilityCurve_(nullptr)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousLiquidPhaseSeparation::d() const
{
    return d_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousLiquidPhaseSeparation::nDot() const
{
    const volScalarField::Internal& alphaSolution =
        mesh().lookupObject<volScalarField::Internal>(alphaNames().first());

    const multicomponentThermo& thermoSolution = specieThermos().first();

    const volScalarField& p = this->p();
    const volScalarField& T = thermoSolution.T();

    const volScalarField::Internal rhoPrecipitate
    (
        specieThermos().valid().second()
      ? vfToVif(specieThermos().second().rhoi(specieis().second(), p, T))
      : vfToVif(thermos().second().rho())
    );

    const volScalarField::Internal v(constant::mathematical::pi/6*pow3(d_));

    return alphaSolution*mDotByAlphaSolution_/(rhoPrecipitate*v);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousLiquidPhaseSeparation::mDot() const
{
    const volScalarField::Internal& alphaSolution =
        mesh().lookupObject<volScalarField::Internal>(alphaNames().first());

    return alphaSolution*mDotByAlphaSolution_;
}


void Foam::fv::homogeneousLiquidPhaseSeparation::correct()
{
    #define DebugField(field)                                                  \
        DebugInfo                                                              \
            << name() << ": "                                                  \
            << #field << ' ' << field.dimensions() << " min/avg/max = "        \
            << gMin(field) << '/' << gAverage(field) << '/' << gMax(field)     \
            << nl;

    using constant::mathematical::pi;
    using constant::physicoChemical::NA;
    using constant::physicoChemical::k;

    const fluidMulticomponentThermo& thermoSolution =
        refCast<const fluidMulticomponentThermo>(specieThermos().first());

    const volScalarField& p = this->p();
    const volScalarField& T = thermoSolution.T();
    DebugField(p);
    DebugField(T);

    // Phase molecular masses and densities
    const volScalarField::Internal rhoSolution(vfToVif(thermoSolution.rho()));
    const volScalarField::Internal muSolution(vfToVif(thermoSolution.mu()));
    const volScalarField::Internal WPrecipitate
    (
        specieThermos().valid().second()
      ? volScalarField::Internal::New
        (
            "W",
            mesh(),
            specieThermos().second().Wi(specieis().second())
        )
      : vfToVif(thermos().second().W())
    );
    const volScalarField::Internal rhoPrecipitate
    (
        specieThermos().valid().second()
      ? vfToVif(specieThermos().second().rhoi(specieis().second(), p, T))
      : vfToVif(thermos().second().rho())
    );
    DebugField(rhoSolution);
    DebugField(WPrecipitate);
    DebugField(rhoPrecipitate);

    // Surface tension
    const volScalarField::Internal sigma(interface_.sigma());
    DebugField(sigma);

    // Mass fraction of nucleating specie
    const volScalarField::Internal Yi = thermoSolution.Y()[specieis().first()];

    // Saturation mass fraction and concentration
    const volScalarField::Internal solubility
    (
        volScalarField::Internal::New
        (
            "YSat",
            mesh(),
            dimless,
            solubilityCurve_->value(T)
        )
    );
    const volScalarField::Internal YSat(solubility/(solubility + 1));
    const volScalarField::Internal cSat(YSat*rhoSolution/WPrecipitate);
    DebugField(YSat);
    DebugField(cSat);

    // Supersaturation of the nucleating specie
    const volScalarField::Internal S(Yi/YSat);
    DebugField(S);

    // Mass and volume of one molecule in the precipitate
    const volScalarField::Internal mMolc(WPrecipitate/NA);
    const volScalarField::Internal vMolc(mMolc/rhoPrecipitate);
    const volScalarField::Internal dMolc(cbrt(6/pi*vMolc));
    DebugField(mMolc);
    DebugField(vMolc);
    DebugField(dMolc);

    // Diameter of nuclei
    d_ = 4*sigma*vMolc/(k*T()*log(max(S, 1 + small)));
    DebugField(d_);

    // ?
    const volScalarField::Internal deltaPhiStar(pi/3*sigma*sqr(d_));
    DebugField(deltaPhiStar);

    // Ratio of nucleus volume to molecular volume
    const volScalarField::Internal iStar(pi/6*pow3(d_)/vMolc);
    DebugField(iStar);

    // Number-based nucleation rate; i.e., number of nuclei created per second
    // per unit volume
    const volScalarField::Internal intermediate(-deltaPhiStar/(k*T()));
    DebugField(intermediate);
    const volScalarField::Internal J
    (
        cSat*NA*exp(-deltaPhiStar/(k*T()))*k*T()/(3*pi*pow3(dMolc)*muSolution)
    );
    DebugField(J);

    // Mass transfer rate
    mDotByAlphaSolution_ = J*iStar*mMolc;
    DebugField(mDotByAlphaSolution_);

    #undef DebugField
}


void Foam::fv::homogeneousLiquidPhaseSeparation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    const label i = index(alphaNames(), eqn.psi().name());

    if (i != -1)
    {
        if (i == 0)
        {
            eqn -= fvm::Sp(mDotByAlphaSolution_, eqn.psi());
        }
        else
        {
            eqn +=
                mDot()
              - correction(fvm::Sp(mDotByAlphaSolution_, eqn.psi()));
        }
    }
    else
    {
        phaseChangeBase::addSup(alpha, rho, eqn);
    }
}


void Foam::fv::homogeneousLiquidPhaseSeparation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& Yi,
    fvMatrix<scalar>& eqn
) const
{
    const label i = index(alphaNames(), eqn.psi().name());

    if (i == 0 && specieis().first() != -1 && Yi.member() == specie())
    {
        eqn -= fvm::Sp(alpha*mDotByAlphaSolution_, Yi);
    }
    else
    {
        singleComponentPhaseChangeBase::addSup(alpha, rho, Yi, eqn);
    }
}


bool Foam::fv::homogeneousLiquidPhaseSeparation::read(const dictionary& dict)
{
    if (singleComponentPhaseChangeBase::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
