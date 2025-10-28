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
    bXiTimedIgnition(name, modelType, mesh, dict),
    zone_(mesh, dict),
    strength_("strength", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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

    const volScalarField& rhou =
        mesh().lookupObject<volScalarField>(IOobject::groupName("rho", "u"));

    scalarField& Sp = eqn.diag();
    const scalarField& V = mesh().V();

    const labelList& cells = zone_.zone();

    const scalar strength = strength_.value();
    const scalar duration = duration_.value();

    forAll(cells, i)
    {
        const label celli = cells[i];
        const scalar Vc = V[celli];
        Sp[celli] -= Vc*rhou[celli]*strength/(duration*(b[celli] + 0.001));
    }
}


void Foam::fv::constantbXiIgnition::topoChange
(
    const polyTopoChangeMap& map
)
{
    zone_.topoChange(map);
}


void Foam::fv::constantbXiIgnition::mapMesh(const polyMeshMap& map)
{
    zone_.mapMesh(map);
}


void Foam::fv::constantbXiIgnition::distribute
(
    const polyDistributionMap& map
)
{
    zone_.distribute(map);
}


bool Foam::fv::constantbXiIgnition::movePoints()
{
    zone_.movePoints();
    return true;
}


bool Foam::fv::constantbXiIgnition::read(const dictionary& dict)
{
    if (bXiIgnition::read(dict))
    {
        zone_.read(coeffs(dict));
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
