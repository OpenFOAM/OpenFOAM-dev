/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

const Foam::NamedEnum
<
    Foam::fv::homogeneousLiquidPhaseSeparation::nucleateType,
    3
>
Foam::fv::homogeneousLiquidPhaseSeparation::nucleateTypeNames_
{"solid", "liquid", "gas"};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::homogeneousLiquidPhaseSeparation::readCoeffs
(
    const dictionary& dict
)
{
    reReadSpecie(dict);

    solubilityCurve_.reset
    (
        Function1<scalar>::New
        (
            "solubility",
            dimTemperature,
            unitFraction,
            dict
        ).ptr()
    );

    nucleateType_ = nucleateTypeNames_.read(dict.lookup("nucleate"));
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
            name() + ":YSat",
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
    phaseChange(name, modelType, mesh, dict, readSpecie(dict, true)),
    fluid_
    (
        mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName)
    ),
    d_
    (
        IOobject
        (
            name + ":d",
            mesh.time().name(),
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
            name + ":mDotByAlpha",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, 0)
    ),
    solubilityCurve_(nullptr),
    nucleateType_(nucleateType::solid)
{
    readCoeffs(coeffs(dict));
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

    const ThermoRefPair<multicomponentThermo> multicomponentThermos =
        this->multicomponentThermos(true, false);

    const multicomponentThermo& thermoSolution = multicomponentThermos.first();

    const volScalarField& p = this->p();
    const volScalarField& T = thermoSolution.T();

    const volScalarField::Internal rhoPrecipitate
    (
        multicomponentThermos.valid().second()
      ? vfToVif(multicomponentThermos.second().rhoi(specieis().second(), p, T))
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


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousLiquidPhaseSeparation::tau() const
{
    static const dimensionedScalar mDotRootVSmall
    (
        dimDensity/dimTime,
        rootVSmall
    );

    const ThermoRefPair<multicomponentThermo> multicomponentThermos =
        this->multicomponentThermos(true, false);

    const multicomponentThermo& thermoSolution = multicomponentThermos.first();

    // Solution density
    const volScalarField::Internal rhoSolution(vfToVif(thermoSolution.rho()));

    // Mass fraction of nucleating specie
    const volScalarField::Internal Yi = thermoSolution.Y()[specieis().first()];

    return Yi*rhoSolution/max(mDotByAlphaSolution_, mDotRootVSmall);
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
    using constant::physicoChemical::NNA;
    using constant::physicoChemical::k;

    const ThermoRefPair<multicomponentThermo> multicomponentThermos =
        this->multicomponentThermos(true, false);

    const fluidMulticomponentThermo& thermoSolution =
        this->fluidMulticomponentThermos(true, false).first();

    const volScalarField& p = this->p();
    const volScalarField& T = thermoSolution.T();
    DebugField(p);
    DebugField(T);

    // Phase molecular masses and densities
    const volScalarField::Internal rhoSolution(vfToVif(thermoSolution.rho()));
    const volScalarField::Internal muSolution(vfToVif(thermoSolution.mu()));
    const volScalarField::Internal WPrecipitate
    (
        multicomponentThermos.valid().second()
      ? volScalarField::Internal::New
        (
            "W",
            mesh(),
            multicomponentThermos.second().Wi(specieis().second())
        )
      : vfToVif(thermos().second().W())
    );
    const volScalarField::Internal rhoPrecipitate
    (
        multicomponentThermos.valid().second()
      ? vfToVif(multicomponentThermos.second().rhoi(specieis().second(), p, T))
      : vfToVif(thermos().second().rho())
    );
    DebugField(rhoSolution);
    DebugField(WPrecipitate);
    DebugField(rhoPrecipitate);

    // Surface tension
    const volScalarField::Internal sigma
    (
        fluid_.sigma
        (
            phaseInterface
            (
                fluid_.phases()[phaseNames().first()],
                fluid_.phases()[phaseNames().second()]
            )
        )
    );
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
    const volScalarField::Internal mMolc(WPrecipitate/NNA);
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

    // Pre-exponential factor. Depends on the type of nucleates.
    tmp<volScalarField::Internal> talpha;
    switch (nucleateType_)
    {
        case nucleateType::solid:
        case nucleateType::liquid:
            talpha = k*T()/(3*pi*pow3(dMolc)*muSolution);
            break;

        case nucleateType::gas:
            talpha = sqrt(2*sigma/(pi*mMolc));
            break;
    }

    // Number-based nucleation rate; i.e., number of nuclei created per second
    // per unit volume
    const volScalarField::Internal J
    (
        cSat*NNA*talpha*exp(-deltaPhiStar/(k*T()))
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

    // !!! Note at present multiphaseEuler cannot linearise w.r.t alphaA in the
    // continuity equation for alphaB. So we can only create a linearised
    // source for this model in the solution volume-fraction equation.

    if (i == 0)
    {
        eqn -= fvm::Sp(mDotByAlphaSolution_, eqn.psi());
    }
    else
    {
        massTransfer::addSup(alpha, rho, eqn);
    }
}


bool Foam::fv::homogeneousLiquidPhaseSeparation::read(const dictionary& dict)
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
