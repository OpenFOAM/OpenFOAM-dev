/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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

#include "solidThermalEquilibrium.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(solidThermalEquilibrium, 0);
    addToRunTimeSelectionTable(fvModel, solidThermalEquilibrium, dictionary);
    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvModel,
        solidThermalEquilibrium,
        dictionary,
        solidEquilibriumEnergySource,
        "solidEquilibriumEnergySource"
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::solidThermalEquilibrium::readCoeffs()
{
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

    solidPhaseName_ = coeffs().lookup<word>("solidPhase");
}


const Foam::volScalarField&
Foam::fv::solidThermalEquilibrium::solidAlpha() const
{
    const word alphaName = IOobject::groupName("alpha", solidPhaseName_);

    if (!mesh().foundObject<volScalarField>(alphaName))
    {
        volScalarField* alphaPtr =
            new volScalarField
            (
                IOobject
                (
                    alphaName,
                    mesh().time().constant(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh()
            );

        alphaPtr->store();
    }

    return mesh().lookupObject<volScalarField>(alphaName);
}


const Foam::solidThermo&
Foam::fv::solidThermalEquilibrium::solidThermo() const
{
    const word thermoName =
        IOobject::groupName(physicalProperties::typeName, solidPhaseName_);

    if (!mesh().foundObject<Foam::solidThermo>(thermoName))
    {
        Foam::solidThermo* thermoPtr =
            solidThermo::New(mesh(), solidPhaseName_).ptr();

        thermoPtr->properties().store();
    }

    return mesh().lookupObject<Foam::solidThermo>(thermoName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::solidThermalEquilibrium::solidThermalEquilibrium
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    phaseName_(word::null),
    solidPhaseName_(word::null)
{
    read(dict);
    solidAlpha();
    solidThermo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::solidThermalEquilibrium::~solidThermalEquilibrium()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList Foam::fv::solidThermalEquilibrium::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        );

    return wordList(1, thermo.he().name());
}


void Foam::fv::solidThermalEquilibrium::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const volScalarField alphahe
    (
        "alphahe",
        solidThermo().kappa()/solidThermo().Cpv()
    );

    const volScalarField& A = solidAlpha();
    const volScalarField B(1 - A);

    eqn -=
        A/B*fvm::ddt(solidThermo().rho(), eqn.psi());
      - 1/B*fvm::laplacian
        (
            A*alphahe,
            eqn.psi(),
            "laplacian(" + alphahe.name() + "," + eqn.psi().name() + ")"
        );
}


void Foam::fv::solidThermalEquilibrium::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const volScalarField alphahe
    (
        "alphahe",
        solidThermo().kappa()/solidThermo().Cpv()
    );

    const volScalarField& A = solidAlpha();
    const volScalarField B(1 - A);

    eqn -=
        A/B*fvm::ddt(alpha, solidThermo().rho(), eqn.psi());
      - 1/B*fvm::laplacian
        (
            A*alphahe,
            eqn.psi(),
            "laplacian(" + alphahe.name() + "," + eqn.psi().name() + ")"
        );
}


bool Foam::fv::solidThermalEquilibrium::movePoints()
{
    return true;
}


void Foam::fv::solidThermalEquilibrium::topoChange
(
    const polyTopoChangeMap&
)
{}


void Foam::fv::solidThermalEquilibrium::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::solidThermalEquilibrium::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fv::solidThermalEquilibrium::read(const dictionary& dict)
{
    if (fvModel::read(dict))
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
