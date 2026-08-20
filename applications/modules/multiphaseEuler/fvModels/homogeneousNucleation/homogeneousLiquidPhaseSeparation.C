/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2026 OpenFOAM Foundation
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
            dimensions::temperature,
            units::fraction,
            dict
        ).ptr()
    );
}


Foam::Pair<Foam::tmp<Foam::volInternalScalarField>>
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
    const volInternalScalarField rhoSolution(vfToVif(thermoSolution.rho()));
    const volInternalScalarField WPrecipitate
    (
        multicomponentThermos.valid().second()
      ? volInternalScalarField::New
        (
            "W",
            mesh(),
            multicomponentThermos.second().Wi(specieis().second())
        )
      : vfToVif(thermos().second().W())
    );
    const volInternalScalarField rhoPrecipitate
    (
        multicomponentThermos.valid().second()
      ? vfToVif(multicomponentThermos.second().rhoi(specieis().second(), p, T))
      : vfToVif(thermos().second().rho())
    );
    infoFieldVariable(rhoSolution, debug);
    infoFieldVariable(WPrecipitate, debug);
    infoFieldVariable(rhoPrecipitate, debug);

    // Viscosity
    const volInternalScalarField muSolution(vfToVif(thermoSolution.mu()));

    // Surface tension
    const volInternalScalarField sigma(this->sigma());
    infoFieldVariable(sigma, debug);

    // Mass fraction of nucleating specie
    const volInternalScalarField Yi = thermoSolution.Y()[specieis().first()];

    // Saturation mass fraction and concentration
    const volInternalScalarField solubility
    (
        volInternalScalarField::New
        (
            "YSat",
            mesh(),
            dimless,
            solubilityCurve_->value(T)
        )
    );
    const volInternalScalarField YSat(solubility/(solubility + 1));
    const volInternalScalarField cSat(YSat*rhoSolution/WPrecipitate);
    infoFieldVariable(YSat, debug);
    infoFieldVariable(cSat, debug);

    // Supersaturation of the nucleating specie
    const volInternalScalarField S(Yi/YSat);
    infoFieldVariable(S, true);

    // Mass and diameter of one molecule in the precipitate
    const volInternalScalarField mMolc(WPrecipitate/NNA);
    const volInternalScalarField dMolc(cbrt(6/pi*(mMolc/rhoPrecipitate)));
    infoFieldVariable(mMolc, debug);
    infoFieldVariable(dMolc, debug);

    // Diameter and mass of a nucleus
    tmp<volInternalScalarField> td =
        4*sigma*mMolc/rhoPrecipitate/(k*T()*log(max(S, 1 + small)));
    const volInternalScalarField& d = td();
    const volInternalScalarField m(pi/6*pow3(d)*rhoPrecipitate);
    infoFieldVariable(d, true);
    infoField("m/mMolc", m/mMolc, debug);

    // Free energy cost of a nucleus
    const volInternalScalarField deltaPhiStar(pi/3*sigma*sqr(d));
    infoFieldVariable(deltaPhiStar, debug);

    // Number-based nucleation rate; i.e., number of nuclei created per second
    // per unit volume
    const volInternalScalarField J
    (
        cSat*NNA*k*T()/(3*pi*pow3(dMolc)*muSolution)*exp(-deltaPhiStar/(k*T()))
    );
    infoFieldVariable(J, debug);

    return Pair<tmp<volInternalScalarField>>(td, J*m);

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
