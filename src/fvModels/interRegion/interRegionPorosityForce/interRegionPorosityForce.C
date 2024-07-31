/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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

#include "interRegionPorosityForce.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "porosityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(interRegionPorosityForce, 0);
    addToRunTimeSelectionTable(fvModel, interRegionPorosityForce, dictionary);
    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvModel,
        interRegionPorosityForce,
        dictionary,
        interRegionExplicitPorositySource,
        "interRegionExplicitPorositySource"
    );
}
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::fv::interRegionPorosityForce::readCoeffs(const dictionary& dict)
{
    UName_ = dict.lookupOrDefault<word>("U", "U");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::interRegionPorosityForce::interRegionPorosityForce
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    interRegionModel(name, modelType, mesh, dict),
    UName_(word::null),
    filter_
    (
        volScalarField::Internal::New
        (
            "filter",
            mesh,
            dimensionedScalar(dimless, 0)
        )
    ),
    porosityPtr_(nullptr)
{
    readCoeffs(coeffs(dict));

    const fvMesh& nbrMesh = mesh.time().lookupObject<fvMesh>(nbrRegionName());

    interpolate(scalarField(nbrMesh.nCells(), 1), filter_);

    const word zoneName(name + ":porous");

    if (!mesh.cellZones().found(zoneName))
    {
        // Scan the porous region filter for all cells containing porosity
        labelList porousCells(mesh.nCells());

        label i = 0;
        forAll(filter_, celli)
        {
            if (filter_[celli] > small)
            {
                porousCells[i++] = celli;
            }
        }
        porousCells.setSize(i);

        mesh.cellZones().append
        (
            new cellZone
            (
                zoneName,
                porousCells,
                mesh.cellZones()
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unable to create porous cellZone " << zoneName
            << ": zone already exists"
            << abort(FatalError);
    }

    porosityPtr_ = porosityModel::New
    (
        name,
        mesh,
        coeffs(dict),
        zoneName
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::interRegionPorosityForce::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::interRegionPorosityForce::addSup
(
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());
    porosityPtr_->addResistance(porosityEqn);
    eqn -= filter_*porosityEqn;
}


void Foam::fv::interRegionPorosityForce::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());
    porosityPtr_->addResistance(porosityEqn);
    eqn -= filter_*porosityEqn;
}


void Foam::fv::interRegionPorosityForce::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());
    porosityPtr_->addResistance(porosityEqn);
    eqn -= alpha*filter_*porosityEqn;
}


bool Foam::fv::interRegionPorosityForce::movePoints()
{
    NotImplemented;
    return true;
}


void Foam::fv::interRegionPorosityForce::topoChange(const polyTopoChangeMap&)
{
    NotImplemented;
}


void Foam::fv::interRegionPorosityForce::mapMesh(const polyMeshMap&)
{
    NotImplemented;
}


void Foam::fv::interRegionPorosityForce::distribute(const polyDistributionMap&)
{
    NotImplemented;
}


bool Foam::fv::interRegionPorosityForce::read(const dictionary& dict)
{
    if (interRegionModel::read(dict))
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
