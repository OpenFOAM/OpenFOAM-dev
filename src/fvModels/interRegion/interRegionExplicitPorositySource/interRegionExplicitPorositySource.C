/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "interRegionExplicitPorositySource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "porosityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(interRegionExplicitPorositySource, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        interRegionExplicitPorositySource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::fv::interRegionExplicitPorositySource::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("U", "U");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::interRegionExplicitPorositySource::interRegionExplicitPorositySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionModel(name, modelType, dict, mesh),
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
    readCoeffs();

    const fvMesh& nbrMesh = mesh.time().lookupObject<fvMesh>(nbrRegionName());

    meshInterp().mapTgtToSrc
    (
        scalarField(nbrMesh.nCells(), 1),
        plusEqOp<scalar>(),
        filter_
    );

    const word zoneName(name + ":porous");

    const meshCellZones& cellZones = mesh.cellZones();
    label zoneID = cellZones.findZoneID(zoneName);

    if (zoneID == -1)
    {
        meshCellZones& cz = const_cast<meshCellZones&>(cellZones);

        zoneID = cz.size();

        cz.setSize(zoneID + 1);

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

        cz.set
        (
            zoneID,
            new cellZone
            (
                zoneName,
                porousCells,
                zoneID,
                cellZones
            )
        );

        cz.clearAddressing();
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
        coeffs(),
        zoneName
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList
Foam::fv::interRegionExplicitPorositySource::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::interRegionExplicitPorositySource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());
    porosityPtr_->addResistance(porosityEqn);
    eqn -= filter_*porosityEqn;
}


void Foam::fv::interRegionExplicitPorositySource::updateMesh(const mapPolyMesh&)
{}


void Foam::fv::interRegionExplicitPorositySource::distribute
(
    const mapDistributePolyMesh&
)
{}


bool Foam::fv::interRegionExplicitPorositySource::movePoints()
{
    return true;
}


bool Foam::fv::interRegionExplicitPorositySource::read(const dictionary& dict)
{
    if (interRegionModel::read(dict))
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
