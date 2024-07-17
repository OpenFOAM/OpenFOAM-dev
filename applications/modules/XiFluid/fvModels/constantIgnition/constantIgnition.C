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

#include "constantIgnition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(constantIgnition, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            constantIgnition,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::constantIgnition::readCoeffs()
{
    duration_ = coeffs().lookup<scalar>("duration", mesh().time().userUnits());
    strength_ = coeffs().lookup<scalar>("strength", dimless);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::constantIgnition::constantIgnition
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    XiFluid_(mesh.lookupObject<solvers::XiFluid>(solver::typeName)),
    set_(mesh, coeffs())
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::constantIgnition::addSupFields() const
{
    if (XiFluid_.ignited)
    {
        return wordList({"b"});
    }
    else
    {
        return wordList::null();
    }
}


void Foam::fv::constantIgnition::addSup
(
    const volScalarField& rho,
    const volScalarField& b,
    fvMatrix<scalar>& eqn
) const
{
    const scalar curTime = mesh().time().value();
    const scalar deltaT = mesh().time().deltaTValue();

    const bool igniting =
        XiFluid_.ignited
     && (curTime - deltaT
          < XiFluid_.ignitionStart + max(duration_, deltaT) + small);

    if (!igniting) return;

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


void Foam::fv::constantIgnition::topoChange
(
    const polyTopoChangeMap& map
)
{
    set_.topoChange(map);
}


void Foam::fv::constantIgnition::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::constantIgnition::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
}


bool Foam::fv::constantIgnition::movePoints()
{
    set_.movePoints();
    return true;
}


bool Foam::fv::constantIgnition::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        set_.read(coeffs());
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
