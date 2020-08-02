/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "radialActuationDiskSource.H"
#include "volFields.H"
#include "fvMatrix.H"
#include "fvm.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::radialActuationDiskSource::
addRadialActuationDiskAxialInertialResistance
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    scalar a = 1.0 - Cp_/Ct_;
    scalarField Tr(cells.size());
    const vector uniDiskDir = diskDir_/mag(diskDir_);

    tensor E(Zero);
    E.xx() = uniDiskDir.x();
    E.yy() = uniDiskDir.y();
    E.zz() = uniDiskDir.z();

    const Field<vector> zoneCellCentres(mesh().cellCentres(), cells);
    const Field<scalar> zoneCellVolumes(mesh().cellVolumes(), cells);

    const vector avgCentre = gSum(zoneCellVolumes*zoneCellCentres)/V();
    const scalar maxR = gMax(mag(zoneCellCentres - avgCentre));

    scalar intCoeffs =
        radialCoeffs_[0]
      + radialCoeffs_[1]*sqr(maxR)/2.0
      + radialCoeffs_[2]*pow4(maxR)/3.0;

    vector upU = vector(vGreat, vGreat, vGreat);
    scalar upRho = vGreat;
    if (upstreamCellId_ != -1)
    {
        upU =  U[upstreamCellId_];
        upRho = rho[upstreamCellId_];
    }
    reduce(upU, minOp<vector>());
    reduce(upRho, minOp<scalar>());

    scalar T = 2.0*upRho*diskArea_*mag(upU)*a*(1.0 - a);
    forAll(cells, i)
    {
        scalar r2 = magSqr(mesh().cellCentres()[cells[i]] - avgCentre);

        Tr[i] =
            T
           *(radialCoeffs_[0] + radialCoeffs_[1]*r2 + radialCoeffs_[2]*sqr(r2))
           /intCoeffs;

        Usource[cells[i]] += ((Vcells[cells[i]]/V())*Tr[i]*E) & upU;
    }

    if (debug)
    {
        Info<< "Source name: " << name() << nl
            << "Average centre: " << avgCentre << nl
            << "Maximum radius: " << maxR << endl;
    }
}


// ************************************************************************* //
