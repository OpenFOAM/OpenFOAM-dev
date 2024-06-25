/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "actuationDisk.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuationDisk, 0);
    addToRunTimeSelectionTable(fvModel, actuationDisk, dictionary);
    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvModel,
        actuationDisk,
        dictionary,
        actuationDiskSource,
        "actuationDiskSource"
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::actuationDisk::readCoeffs()
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


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::actuationDisk::addActuationDiskAxialInertialResistance
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const scalar a = 1 - Cp_/Ct_;
    const vector dHat(diskDir_/mag(diskDir_));

    scalar dHatUo(vGreat);
    if (upstreamCellId_ != -1)
    {
        dHatUo = dHat & U[upstreamCellId_];
    }
    reduce(dHatUo, minOp<scalar>());

    const vector T = 2*diskArea_*sqr(dHatUo)*a*(1 - a)*dHat;

    forAll(cells, i)
    {
        Usource[cells[i]] +=
            (alpha[cells[i]]*rho[cells[i]]*(Vcells[cells[i]]/set_.V()))*T;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuationDisk::actuationDisk
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

Foam::wordList Foam::fv::actuationDisk::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::actuationDisk::addSup
(
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    addActuationDiskAxialInertialResistance
    (
        eqn.source(),
        set_.cells(),
        mesh().V(),
        geometricOneField(),
        geometricOneField(),
        U
    );
}


void Foam::fv::actuationDisk::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    addActuationDiskAxialInertialResistance
    (
        eqn.source(),
        set_.cells(),
        mesh().V(),
        geometricOneField(),
        rho,
        U
    );
}


void Foam::fv::actuationDisk::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    addActuationDiskAxialInertialResistance
    (
        eqn.source(),
        set_.cells(),
        mesh().V(),
        alpha,
        rho,
        U
    );
}


bool Foam::fv::actuationDisk::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::actuationDisk::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::actuationDisk::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::actuationDisk::distribute(const polyDistributionMap& map)
{
    set_.distribute(map);
}


bool Foam::fv::actuationDisk::read(const dictionary& dict)
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
