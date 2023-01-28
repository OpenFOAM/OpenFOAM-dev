/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "actuationDiskSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        actuationDiskSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::actuationDiskSource::readCoeffs()
{
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

    UName_ =
        coeffs().lookupOrDefault<word>
        (
            "U",
            IOobject::groupName("U", phaseName_)
        );

    diskDir_ = coeffs().lookup<vector>("diskDir");
    if (mag(diskDir_) < vSmall)
    {
        FatalIOErrorInFunction(coeffs())
           << "disk direction vector is approximately zero"
           << exit(FatalIOError);
    }

    Cp_ = coeffs().lookup<scalar>("Cp");
    Ct_ = coeffs().lookup<scalar>("Ct");
    if (Cp_ <= vSmall || Ct_ <= vSmall)
    {
        FatalIOErrorInFunction(coeffs())
           << "Cp and Ct must be greater than zero"
           << exit(FatalIOError);
    }

    diskArea_ = coeffs().lookup<scalar>("diskArea");
    if (magSqr(diskArea_) <= vSmall)
    {
        FatalIOErrorInFunction(coeffs())
           << "diskArea is approximately zero"
           << exit(FatalIOError);
    }

    upstreamPoint_ = coeffs().lookup<point>("upstreamPoint");
    upstreamCellId_ = mesh().findCell(upstreamPoint_);
    if (returnReduce(upstreamCellId_, maxOp<label>()) == -1)
    {
        FatalIOErrorInFunction(coeffs())
           << "upstream location " << upstreamPoint_  << " not found in mesh"
           << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuationDiskSource::actuationDiskSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    set_(mesh, coeffs()),
    phaseName_(word::null),
    UName_(word::null),
    diskDir_(vector::uniform(NaN)),
    Cp_(NaN),
    Ct_(NaN),
    diskArea_(NaN),
    upstreamPoint_(vector::uniform(NaN)),
    upstreamCellId_(-1)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::actuationDiskSource::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::actuationDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    const scalarField& cellsV = mesh().V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    addActuationDiskAxialInertialResistance
    (
        Usource,
        set_.cells(),
        cellsV,
        geometricOneField(),
        geometricOneField(),
        U
    );
}


void Foam::fv::actuationDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    const scalarField& cellsV = mesh().V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    addActuationDiskAxialInertialResistance
    (
        Usource,
        set_.cells(),
        cellsV,
        geometricOneField(),
        rho,
        U
    );
}


void Foam::fv::actuationDiskSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    const scalarField& cellsV = mesh().V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    addActuationDiskAxialInertialResistance
    (
        Usource,
        set_.cells(),
        cellsV,
        alpha,
        rho,
        U
    );
}


bool Foam::fv::actuationDiskSource::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::actuationDiskSource::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::actuationDiskSource::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::actuationDiskSource::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
}


bool Foam::fv::actuationDiskSource::read(const dictionary& dict)
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
}


// ************************************************************************* //
