/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
        option,
        interRegionExplicitPorositySource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::interRegionExplicitPorositySource::initialise() const
{
    if (!firstIter_)
    {
        return;
    }

    const word zoneName(name_ + ":porous");

    const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName_);
    const cellZoneMesh& cellZones = nbrMesh.cellZones();
    label zoneID = cellZones.findZoneID(zoneName);

    if (zoneID == -1)
    {
        cellZoneMesh& cz = const_cast<cellZoneMesh&>(cellZones);

        zoneID = cz.size();

        cz.setSize(zoneID + 1);

        cz.set
        (
            zoneID,
            new cellZone
            (
                zoneName,
                nbrMesh.faceNeighbour(), // Neighbour internal cells
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

    porosityPtr_.reset
    (
        porosityModel::New
        (
            name_,
            nbrMesh,
            coeffs_,
            zoneName
        ).ptr()
    ),

    firstIter_ = false;
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
    interRegionOption(name, modelType, dict, mesh),
    porosityPtr_(nullptr),
    firstIter_(true),
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    muName_(coeffs_.lookupOrDefault<word>("mu", "thermo:mu"))
{
    if (active_)
    {
        fieldNames_.setSize(1, UName_);
        applied_.setSize(1, false);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::interRegionExplicitPorositySource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
) const
{
    initialise();

    const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName_);

    const volVectorField& U = eqn.psi();

    volVectorField UNbr
    (
        IOobject
        (
            name_ + ":UNbr",
            nbrMesh.time().timeName(),
            nbrMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nbrMesh,
        dimensionedVector(U.dimensions(), Zero)
    );

    // Map local velocity onto neighbour region
    meshInterp().mapSrcToTgt
    (
        U.primitiveField(),
        plusEqOp<vector>(),
        UNbr.primitiveFieldRef()
    );

    fvMatrix<vector> nbrEqn(UNbr, eqn.dimensions());

    porosityPtr_->addResistance(nbrEqn);

    // Convert source from neighbour to local region
    fvMatrix<vector> porosityEqn(U, eqn.dimensions());
    scalarField& Udiag = porosityEqn.diag();
    vectorField& Usource = porosityEqn.source();

    Udiag.setSize(eqn.diag().size(), 0.0);
    Usource.setSize(eqn.source().size(), Zero);

    meshInterp().mapTgtToSrc(nbrEqn.diag(), plusEqOp<scalar>(), Udiag);
    meshInterp().mapTgtToSrc(nbrEqn.source(), plusEqOp<vector>(), Usource);

    eqn -= porosityEqn;
}


void Foam::fv::interRegionExplicitPorositySource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
) const
{
    initialise();

    const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName_);

    const volVectorField& U = eqn.psi();

    volVectorField UNbr
    (
        IOobject
        (
            name_ + ":UNbr",
            nbrMesh.time().timeName(),
            nbrMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nbrMesh,
        dimensionedVector(U.dimensions(), Zero)
    );

    // Map local velocity onto neighbour region
    meshInterp().mapSrcToTgt
    (
        U.primitiveField(),
        plusEqOp<vector>(),
        UNbr.primitiveFieldRef()
    );

    fvMatrix<vector> nbrEqn(UNbr, eqn.dimensions());

    volScalarField rhoNbr
    (
        IOobject
        (
            "rho:UNbr",
            nbrMesh.time().timeName(),
            nbrMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nbrMesh,
        dimensionedScalar(dimDensity, 0)
    );

    volScalarField muNbr
    (
        IOobject
        (
            "mu:UNbr",
            nbrMesh.time().timeName(),
            nbrMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nbrMesh,
        dimensionedScalar(dimViscosity, 0)
    );

    const volScalarField& mu =
        mesh_.lookupObject<volScalarField>(muName_);

    // Map local rho onto neighbour region
    meshInterp().mapSrcToTgt
    (
        rho.primitiveField(),
        plusEqOp<scalar>(),
        rhoNbr.primitiveFieldRef()
    );

    // Map local mu onto neighbour region
    meshInterp().mapSrcToTgt
    (
        mu.primitiveField(),
        plusEqOp<scalar>(),
        muNbr.primitiveFieldRef()
    );

    porosityPtr_->addResistance(nbrEqn, rhoNbr, muNbr);

    // Convert source from neighbour to local region
    fvMatrix<vector> porosityEqn(U, eqn.dimensions());
    scalarField& Udiag = porosityEqn.diag();
    vectorField& Usource = porosityEqn.source();

    Udiag.setSize(eqn.diag().size(), 0.0);
    Usource.setSize(eqn.source().size(), Zero);

    meshInterp().mapTgtToSrc(nbrEqn.diag(), plusEqOp<scalar>(), Udiag);
    meshInterp().mapTgtToSrc(nbrEqn.source(), plusEqOp<vector>(), Usource);

    eqn -= porosityEqn;
}


bool Foam::fv::interRegionExplicitPorositySource::read(const dictionary& dict)
{
    if (interRegionOption::read(dict))
    {
        coeffs_.readIfPresent("U", UName_);
        coeffs_.readIfPresent("mu", muName_);

        // Reset the porosity model?

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
