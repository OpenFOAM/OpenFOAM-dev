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
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::constantbXiIgnition::readCoeffs()
{
    start_ = coeffs().lookup<scalar>("start", mesh().time().userUnits());
    duration_ = coeffs().lookup<scalar>("duration", mesh().time().userUnits());
    strength_ = coeffs().lookup<scalar>("strength", dimless);
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
    set_(mesh, coeffs()),
    XiCorrSet_(mesh, coeffs().subDict("XiCorr"))
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::constantbXiIgnition::igniting() const
{
    scalar curTime = mesh().time().value();
    scalar deltaT = mesh().time().deltaTValue();

    return
    (
        (curTime - deltaT >= start_)
     && (curTime - deltaT < start_ + max(duration_, deltaT) + small)
    );
}


bool Foam::fv::constantbXiIgnition::ignited() const
{
    scalar curTime = mesh().time().value();
    scalar deltaT = mesh().time().deltaTValue();

    return (curTime - deltaT >= start_);
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

    forAll(cells, i)
    {
        const label celli = cells[i];
        const scalar Vc = V[celli];
        Sp[celli] -= Vc*rhou[celli]*strength_/(duration_*(b[celli] + 0.001));
    }
}


void Foam::fv::constantbXiIgnition::XiCorr
(
    volScalarField& Xi,
    const volScalarField& b,
    const volScalarField& mgb
) const
{
    if (!igniting()) return;

    // Calculate volume of ignition kernel
    const dimensionedScalar Vk
    (
        "Vk",
        dimVolume,
        gSum((1.0 - b)*mesh().V().primitiveField())
    );

    if (Vk.value() > small)
    {
        // Calculate kernel area from its volume
        const dimensionedScalar Ak(this->Ak(Vk));

        // Calculate kernel area from b field
        const dimensionedScalar AkEst
        (
            gSum(mgb*mesh().V().primitiveField())
        );

        const scalar XiCorr = max(min((Ak/AkEst).value(), 10.0), 1.0);
        Info<< "XiCorr = " << XiCorr << endl;

        Xi *= XiCorr;
    }
}


void Foam::fv::constantbXiIgnition::topoChange
(
    const polyTopoChangeMap& map
)
{
    set_.topoChange(map);
}


void Foam::fv::constantbXiIgnition::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::constantbXiIgnition::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
}


bool Foam::fv::constantbXiIgnition::movePoints()
{
    set_.movePoints();
    return true;
}


bool Foam::fv::constantbXiIgnition::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        set_.read(coeffs());
        XiCorrSet_.read(coeffs().subDict("XiCorr"));
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }

    return false;
}


// ************************************************************************* //
