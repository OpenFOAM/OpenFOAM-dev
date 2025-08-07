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

void Foam::fv::homogeneousLiquidPhaseSeparation::readCoeffs
(
    const dictionary& dict
)
{
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
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::fv::homogeneousLiquidPhaseSeparation::dAndMDotByAlphaSolution() const
{
    #define infoFieldVariable(field, print) infoField(#field, field, print)

    using constant::mathematical::pi;
    using constant::physicoChemical::NNA;
    using constant::physicoChemical::k;

    const ThermoRefPair<multicomponentThermo> multicomponentThermos =
        this->multicomponentThermos(true, false);

    const fluidMulticomponentThermo& thermoSolution =
        this->fluidMulticomponentThermos(true, false).first();

    const volScalarField& p = this->p();
    const volScalarField& T = thermoSolution.T();
    infoFieldVariable(p, debug);
    infoFieldVariable(T, debug);

    // Phase molecular masses and densities
    const volScalarField::Internal rhoSolution(vfToVif(thermoSolution.rho()));
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
    infoFieldVariable(rhoSolution, debug);
    infoFieldVariable(WPrecipitate, debug);
    infoFieldVariable(rhoPrecipitate, debug);

    // Viscosity
    const volScalarField::Internal muSolution(vfToVif(thermoSolution.mu()));

    // Surface tension
    const volScalarField::Internal sigma(this->sigma());
    infoFieldVariable(sigma, debug);

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
    infoFieldVariable(YSat, debug);
    infoFieldVariable(cSat, debug);

    // Supersaturation of the nucleating specie
    const volScalarField::Internal S(Yi/YSat);
    infoFieldVariable(S, true);

    // Mass and diameter of one molecule in the precipitate
    const volScalarField::Internal mMolc(WPrecipitate/NNA);
    const volScalarField::Internal dMolc(cbrt(6/pi*(mMolc/rhoPrecipitate)));
    infoFieldVariable(mMolc, debug);
    infoFieldVariable(dMolc, debug);

    // Diameter and mass of a nucleus
    tmp<volScalarField::Internal> td =
        4*sigma*mMolc/rhoPrecipitate/(k*T()*log(max(S, 1 + small)));
    const volScalarField::Internal& d = td();
    const volScalarField::Internal m(pi/6*pow3(d)*rhoPrecipitate);
    infoFieldVariable(d, true);
    infoField("m/mMolc", m/mMolc, debug);

    // Free energy cost of a nucleus
    const volScalarField::Internal deltaPhiStar(pi/3*sigma*sqr(d));
    infoFieldVariable(deltaPhiStar, debug);

    // Number-based nucleation rate; i.e., number of nuclei created per second
    // per unit volume
    const volScalarField::Internal J
    (
        cSat*NNA*k*T()/(3*pi*pow3(dMolc)*muSolution)*exp(-deltaPhiStar/(k*T()))
    );
    infoFieldVariable(J, debug);

    return Pair<tmp<volScalarField::Internal>>(td, J*m);

    #undef infoFieldVariable
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
    homogeneousNucleation(name, modelType, mesh, dict),
    solubilityCurve_(nullptr)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::homogeneousLiquidPhaseSeparation::read(const dictionary& dict)
{
    if (homogeneousNucleation::read(dict))
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
