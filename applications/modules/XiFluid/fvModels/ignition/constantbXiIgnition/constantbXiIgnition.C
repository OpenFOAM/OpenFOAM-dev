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

#include "constantbXiIgnition.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(constantbXiIgnition, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            constantbXiIgnition,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::constantbXiIgnition::readCoeffs(const dictionary& dict)
{
    start_.read(dict, mesh().time().userUnits());
    duration_.read(dict, mesh().time().userUnits());
    strength_.read(dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::constantbXiIgnition::constantbXiIgnition
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    bXiIgnition(name, modelType, mesh, dict),
    set_(mesh, dict),
    XiCorrModel_(XiCorrModel::New(mesh, dict)),
    start_("start", mesh().time().userUnits(), dict),
    duration_("duration", mesh().time().userUnits(), dict),
    strength_("strength", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::constantbXiIgnition::igniting() const
{
    const dimensionedScalar curTime = mesh().time();
    const dimensionedScalar deltaT = mesh().time().deltaT();

    return
    (
        (curTime > start_ - 0.5*deltaT)
     && (curTime < start_ + max(duration_, deltaT))
    );
}


bool Foam::fv::constantbXiIgnition::ignited() const
{
    const dimensionedScalar curTime = mesh().time();
    const dimensionedScalar deltaT = mesh().time().deltaT();

    return (curTime > start_ - 0.5*deltaT);
}


void Foam::fv::constantbXiIgnition::addSup
(
    const volScalarField& rho,
    const volScalarField& b,
    fvMatrix<scalar>& eqn
) const
{
    if (!igniting()) return;

    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const volScalarField& rhou = mesh().lookupObject<volScalarField>("rhou");

    scalarField& Sp = eqn.diag();
    const scalarField& V = mesh().V();

    const labelUList cells = set_.cells();

    const scalar strength = strength_.value();
    const scalar duration = duration_.value();

    forAll(cells, i)
    {
        const label celli = cells[i];
        const scalar Vc = V[celli];
        Sp[celli] -= Vc*rhou[celli]*strength/(duration*(b[celli] + 0.001));
    }
}


void Foam::fv::constantbXiIgnition::XiCorr
(
    volScalarField& Xi,
    const volScalarField& b,
    const volScalarField& mgb
) const
{
    if (igniting())
    {
        XiCorrModel_->XiCorr(Xi, b, mgb);
    }
}


void Foam::fv::constantbXiIgnition::topoChange
(
    const polyTopoChangeMap& map
)
{
    set_.topoChange(map);
    XiCorrModel_->topoChange(map);
}


void Foam::fv::constantbXiIgnition::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
    XiCorrModel_->mapMesh(map);
}


void Foam::fv::constantbXiIgnition::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
    XiCorrModel_->distribute(map);
}


bool Foam::fv::constantbXiIgnition::movePoints()
{
    set_.movePoints();
    XiCorrModel_->movePoints();
    return true;
}


bool Foam::fv::constantbXiIgnition::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        set_.read(coeffs(dict));
        XiCorrModel_->read(coeffs(dict));
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }

    return false;
}


// ************************************************************************* //
