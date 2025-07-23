/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "propellerDiskAdjustment.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::propellerDiskAdjustment::readCoeffs(const dictionary& dict)
{
    forces_ = new functionObjects::forces
    (
        "forces",
        propellerDisk_.mesh().time(),
        dict.subDict("forces")
    );

    deltaTStar_ = dict.lookup<scalar>("deltaTStar");
    Tmin_ = dict.lookup<scalar>("Tmin");
    sfc_ = dict.lookup<scalar>("sfc");
    startTime_ = dict.lookup<scalar>("startTime");
    nFraction_ = dict.lookup<scalar>("nFraction");
    resistanceFraction_ =
        dict.lookupOrDefault<scalar>("resistanceFraction", 1.0);
    resistanceDirection_ = dict.lookupOrDefault<vector>
    (
        "resistanceDirection",
        propellerDisk_.normal_
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::fv::propellerDiskAdjustment::resistance() const
{
    forces_->calcForcesMoments();
    return -resistanceFraction_*(forces_->forceEff() & resistanceDirection_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::propellerDiskAdjustment::propellerDiskAdjustment
(
    const propellerDisk& pd,
    const dictionary& dict
)
:
    propellerDisk_(pd),
    n_
    (
        IOobject
        (
            propellerDisk_.typedName("n"),
            propellerDisk_.mesh().time().name(),
            propellerDisk_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar("n", dimless/dimTime, mag(dict.lookup<scalar>("n")))
    )
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::propellerDiskAdjustment::correctn(const scalar T) const
{
    const fvMesh& mesh = propellerDisk_.mesh();

    if (mesh.time().value() > startTime_)
    {
        const scalar& n = n_.value();

        const fvCellZone& zone = propellerDisk_.zone_;
        const labelList& zoneCells = zone.zone();

        // Average propeller disk density
        const scalarField zoneCellVolumes(mesh.cellVolumes(), zoneCells);
        const scalarField zoneRho(forces_->rho(), zoneCells);
        const scalar rho = gSum(zoneCellVolumes*zoneRho)/zone.V();

        // Get ship resistance
        const scalar res = resistance();

        const scalar beta = max(Tmin_, T*rho)/(sqr(n));
        const scalar diffForce = (res - sfc_ - T*rho);

        scalar deltaN =
            mesh.time().deltaT().value()*(diffForce/(2*beta*n*deltaTStar_));
        const scalar magDeltaN = min(mag(deltaN), nFraction_*n);
        deltaN = sign(deltaN)*magDeltaN;

        n_ = n_.oldTime() + dimensionedScalar(dimRate, deltaN);
    }
}


// ************************************************************************* //
