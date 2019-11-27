/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "solidEquilibriumEnergySource.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(solidEquilibriumEnergySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        solidEquilibriumEnergySource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::volScalarField&
Foam::fv::solidEquilibriumEnergySource::alpha() const
{
    const word alphaName = IOobject::groupName("alpha", phaseName_);

    if (!mesh_.foundObject<volScalarField>(alphaName))
    {
        volScalarField* alphaPtr =
            new volScalarField
            (
                IOobject
                (
                    alphaName,
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            );

        alphaPtr->store();
    }

    return mesh_.lookupObject<volScalarField>(alphaName);
}


const Foam::solidThermo& Foam::fv::solidEquilibriumEnergySource::thermo() const
{
    const word thermoName =
        IOobject::groupName(basicThermo::dictName, phaseName_);

    if (!mesh_.foundObject<solidThermo>(thermoName))
    {
        solidThermo* thermoPtr = solidThermo::New(mesh_, phaseName_).ptr();

        thermoPtr->store();
    }

    return mesh_.lookupObject<solidThermo>(thermoName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::solidEquilibriumEnergySource::solidEquilibriumEnergySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    phaseName_(dict.lookup<word>("phase"))
{
    read(dict);
    alpha();
    thermo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::solidEquilibriumEnergySource::~solidEquilibriumEnergySource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fv::solidEquilibriumEnergySource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volScalarField alphahe(thermo().alphahe());

    const volScalarField& A = this->alpha();
    const volScalarField B(1 - A);

    eqn -=
        A/B*fvm::ddt(thermo().rho(), eqn.psi());
      - 1/B*fvm::laplacian
        (
            A*alphahe,
            eqn.psi(),
            "laplacian(" + alphahe.name() + "," + eqn.psi().name() + ")"
        );
}


void Foam::fv::solidEquilibriumEnergySource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volScalarField alphahe(alpha*thermo().alphahe());

    const volScalarField& A = this->alpha();
    const volScalarField B(1 - A);

    eqn -=
        A/B*fvm::ddt(alpha, thermo().rho(), eqn.psi());
      - 1/B*fvm::laplacian
        (
            A*alphahe,
            eqn.psi(),
            "laplacian(" + alphahe.name() + "," + eqn.psi().name() + ")"
        );
}


bool Foam::fv::solidEquilibriumEnergySource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        fieldNames_ = wordList(1, coeffs_.lookup<word>("field"));

        applied_.setSize(fieldNames_.size(), false);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
