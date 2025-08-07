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

#include "propellerDisk.H"
#include "propellerDiskAdjustment.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(propellerDisk, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        propellerDisk,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::propellerDisk::readCoeffs(const dictionary& dict)
{
    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    UName_ = dict.lookupOrDefault<word>
    (
        "U",
        IOobject::groupName("U", phaseName_)
    );

    if (dict.found("centre"))
    {
        centre_ = dict.lookup<vector>("centre");
    }
    else
    {
        const Field<vector> zoneCellCentres(mesh().cellCentres(), zone_.zone());
        const Field<scalar> zoneCellVolumes(mesh().cellVolumes(), zone_.zone());
        centre_ = gSum(zoneCellVolumes*zoneCellCentres)/zone_.V();
    }

    normal_ = normalised(dict.lookup<vector>("normal"));

    n_ = dict.lookup<scalar>("n");
    rotationDir_ = sign(n_);
    n_ = mag(n_);

    dProp_ = dict.lookup<scalar>("dPropeller");
    dHub_ = dict.lookup<scalar>("dHub");

    propellerFunction_ =
        Function1<vector2D>::New("propellerCurve", dimless, dimless, dict);

    log_ = dict.lookupOrDefault<Switch>("log", false);

    if (log_ && Pstream::master())
    {
        logFile_ = new functionObjects::logFile(mesh(), name(), typeName);
        logFile_->writeTimeColumnHeaders
        (
            {
                "n",
                "J",
                "Jcorr",
                "Udisk",
                "Ucorr",
                "Kt",
                "Kq",
                "T/rho",
                "Q/rho",
                "force",
                "moment"
            }
        );
    }

    bool adjustment(dict.lookupOrDefault<Switch>("adjustment", false));

    if (adjustment)
    {
        if (adjustment_.valid())
        {
            adjustment_->readCoeffs(dict);
        }
        else
        {
            adjustment_ = new propellerDiskAdjustment(*this, dict);
        }
    }
    else
    {
        adjustment_.clear();
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::fv::propellerDisk::diskThickness(const vector& centre) const
{
    const Field<vector> zoneCellCentres(mesh().cellCentres(), zone_.zone());
    const scalar r2 = gMax(magSqr(zoneCellCentres - centre));

    return zone_.V()/(constant::mathematical::pi*r2);
}


Foam::scalar Foam::fv::propellerDisk::J
(
    const vectorField& U,
    const vector& nHat
) const
{
    const labelList& cells = zone_.zone();
    const scalarField& V = mesh().V();

    scalar VUn = 0;
    forAll(cells, i)
    {
        VUn += V[cells[i]]*(nHat & U[cells[i]]);
    }
    reduce(VUn, sumOp<scalar>());

    const scalar Uref = VUn/zone_.V();

    return mag(Uref/(n_*dProp_));
}


Foam::scalar Foam::fv::propellerDisk::n() const
{
    if (adjustment_.valid())
    {
        return adjustment_->n();
    }
    else
    {
        return n_;
    }
}


void Foam::fv::propellerDisk::correctn(const scalar T) const
{
    if (adjustment_.valid())
    {
        adjustment_->correctn(T);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::propellerDisk::propellerDisk
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    zone_(mesh, dict),
    phaseName_(word::null),
    UName_(word::null),
    force_(Zero),
    moment_(Zero),
    log_(false)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::propellerDisk::~propellerDisk()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::propellerDisk::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        zone_.read(coeffs(dict));
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


Foam::wordList Foam::fv::propellerDisk::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::propellerDisk::addSup
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
        alpha,
        rho,
        U
    );
}


void Foam::fv::propellerDisk::addSup
(
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    addActuationDiskAxialInertialResistance
    (
        eqn.source(),
        geometricOneField(),
        geometricOneField(),
        U
    );
}


void Foam::fv::propellerDisk::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    addActuationDiskAxialInertialResistance
    (
        eqn.source(),
        geometricOneField(),
        rho,
        U
    );
}


bool Foam::fv::propellerDisk::movePoints()
{
    zone_.movePoints();
    return true;
}


void Foam::fv::propellerDisk::topoChange(const polyTopoChangeMap& map)
{
    zone_.topoChange(map);
}


void Foam::fv::propellerDisk::mapMesh(const polyMeshMap& map)
{
    zone_.mapMesh(map);
}


void Foam::fv::propellerDisk::distribute
(
    const polyDistributionMap& map
)
{
    zone_.distribute(map);
}


// ************************************************************************* //
