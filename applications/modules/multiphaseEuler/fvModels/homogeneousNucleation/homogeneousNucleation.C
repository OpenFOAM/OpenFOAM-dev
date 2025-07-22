/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "homogeneousNucleation.H"
#include "multicomponentThermo.H"
#include "multiphaseEuler.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(homogeneousNucleation, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::homogeneousNucleation::readCoeffs(const dictionary& dict)
{
    reReadSpecie(dict);
}


void Foam::fv::homogeneousNucleation::correctDAndMDot() const
{
    Info<< type() << ": " << name() << endl << incrIndent;

    Pair<tmp<volScalarField::Internal>> dAndMDotByAlphaSolution =
        this->dAndMDotByAlphaSolution();

    d_ = dAndMDotByAlphaSolution.first();

    mDotByAlphaSolution_ = dAndMDotByAlphaSolution.second();
    infoField("mDotByAlphaSolution", mDotByAlphaSolution_, debug);

    const volScalarField::Internal alphaSolution =
        mesh().lookupObject<volScalarField::Internal>(alphaNames().first());

    mDot_ = alphaSolution*mDotByAlphaSolution_;
    infoField("mDot", mDot_);

    Info<< decrIndent;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousNucleation::sigma() const
{
    return vfToVif(fluid_.sigma(phaseInterface(solution_, nucleate_)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::homogeneousNucleation::homogeneousNucleation
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
    solution_(fluid_.phases()[phaseNames().first()]),
    nucleate_(fluid_.phases()[phaseNames().second()]),
    p_rgh_
    (
        mesh().lookupObject<solvers::multiphaseEuler>(solver::typeName).p_rgh
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
            IOobject::groupName(name + ":mDotByAlpha", solution_.name()),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, 0)
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
        dimensionedScalar(dimDensity/dimTime, 0)
    )
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousNucleation::d() const
{
    return d_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousNucleation::nDot() const
{
    const ThermoRefPair<multicomponentThermo> multicomponentThermos =
        this->multicomponentThermos(true, false);

    const multicomponentThermo& thermoSolution = multicomponentThermos.first();

    const volScalarField& p = this->p();
    const volScalarField& T = thermoSolution.T();

    const volScalarField::Internal rhoNucleate
    (
        multicomponentThermos.valid().second()
      ? vfToVif(multicomponentThermos.second().rhoi(specieis().second(), p, T))
      : vfToVif(thermos().second().rho())
    );

    const volScalarField::Internal v(constant::mathematical::pi/6*pow3(d_));

    return mDot_/(rhoNucleate*v);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousNucleation::mDot() const
{
    return mDot_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousNucleation::tau() const
{
    static const dimensionedScalar mDotRootVSmall
    (
        dimDensity/dimTime,
        rootVSmall
    );

    const ThermoRefPair<multicomponentThermo> multicomponentThermos =
        this->multicomponentThermos(true, false);

    const multicomponentThermo& thermoSolution = multicomponentThermos.first();

    // Solution molecular mass and density
    const volScalarField::Internal WSolution(vfToVif(thermoSolution.W()));
    const volScalarField::Internal rhoSolution(vfToVif(thermoSolution.rho()));

    // Mass and mole fractions of the nucleating specie in the fluid
    const volScalarField::Internal& Yi = thermoSolution.Y()[specieis().first()];
    const volScalarField::Internal Xi
    (
        Yi*WSolution/thermoSolution.Wi(specieis().first())
    );

    return Xi*rhoSolution/max(mDotByAlphaSolution_, mDotRootVSmall);
}


void Foam::fv::homogeneousNucleation::correct()
{
    // Reset the p_rgh equation solution counter
    pressureEquationIndex_ = 0;

    // Correct the total phase change rate
    correctDAndMDot();
}


void Foam::fv::homogeneousNucleation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    // Pressure equation (i.e., continuity, linearised in the pressure)
    if
    (
        (&alpha == &solution_ || &alpha == &nucleate_)
     && (&rho == &solution_.rho() || &rho == &nucleate_.rho())
     && &eqn.psi() == &p_rgh_
    )
    {
        // Ensure that the source is up-to date if this is the first call in
        // the current phase loop
        if (pressureEquationIndex_ % 2 == 0) correctDAndMDot();
        pressureEquationIndex_ ++;
    }

    // Let the base class add the actual source
    massTransfer::addSup(alpha, rho, eqn);
}


bool Foam::fv::homogeneousNucleation::read(const dictionary& dict)
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
