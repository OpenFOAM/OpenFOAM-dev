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

#include "cellPointLagrangianAddressor.H"
#include "remote.H"
#include "globalIndex.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellPointLagrangianAddressor, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellPointLagrangianAddressor::cellPointLagrangianAddressor
(
    const polyMesh& mesh
)
:
    DemandDrivenMeshObject
    <
        polyMesh,
        MoveableMeshObject,
        cellPointLagrangianAddressor
    >(mesh),
    pointCoupledPoint_(mesh.nPoints(), -1),
    coupledPointPoint_(),
    coupledPointCoupledPoint_()
{
    // Enumerate all the coupled points, and create forward and reverse maps
    // from the mesh point indices
    {
        DynamicList<label> coupledPointPoint;
        forAll(mesh.boundaryMesh(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];

            if (pp.coupled())
            {
                forAll(pp.meshPoints(), patchPointi)
                {
                    const label pointi = pp.meshPoints()[patchPointi];

                    if (pointCoupledPoint_[pointi] == -1)
                    {
                        pointCoupledPoint_[pointi] = coupledPointPoint.size();
                        coupledPointPoint.append(pointi);
                    }
                }
            }
        }
        coupledPointPoint_.transfer(coupledPointPoint);
    }

    // For every coupled point, create a list of the points it is connected to
    {
        List<DynamicList<remote>> coupledPointCoupledPoint
        (
            coupledPointPoint_.size(),
            DynamicList<remote>(1)
        );

        forAll(coupledPointPoint_, coupledPointi)
        {
            coupledPointCoupledPoint[coupledPointi].append
            (
                remote(coupledPointi)
            );
        }

        syncTools::syncPointList
        (
            mesh,
            coupledPointPoint_,
            coupledPointCoupledPoint,
            [](DynamicList<remote>& x, const DynamicList<remote>& y)
            {
                x.append(y);
            },
            DynamicList<remote>()
        );

        forAll(coupledPointCoupledPoint, coupledPointi)
        {
            label iNew = 0;
            forAll(coupledPointCoupledPoint[coupledPointi], iOld)
            {
                if
                (
                    coupledPointCoupledPoint[coupledPointi][iOld]
                 != remote(coupledPointi)
                )
                {
                    coupledPointCoupledPoint[coupledPointi][iNew] =
                        coupledPointCoupledPoint[coupledPointi][iOld];
                    iNew ++;
                }
            }
            coupledPointCoupledPoint[coupledPointi].resize(iNew);
        }

        CompactListList<remote> cpcp(coupledPointCoupledPoint);
        coupledPointCoupledPoint_.transfer(cpcp);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellPointLagrangianAddressor::~cellPointLagrangianAddressor()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::cellPointLagrangianAddressor::movePoints()
{
    return true;
}


void Foam::cellPointLagrangianAddressor::sync
(
    labelList& pointIndex,
    DynamicList<label>& indexPoint
) const
{
    List<DynamicList<label>> remoteCoupledPoints(Pstream::nProcs());

    forAll(indexPoint, indexi)
    {
        const label pointi = indexPoint[indexi];
        const label coupledPointi = pointCoupledPoint_[pointi];

        if (coupledPointi == -1) continue;

        const UList<remote> coupledPoints =
            coupledPointCoupledPoint_[coupledPointi];

        forAll(coupledPoints, i)
        {
            remoteCoupledPoints[coupledPoints[i].proci].append
            (
                coupledPoints[i].elementi
            );
        }
    }

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        forAll(remoteCoupledPoints, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                UOPstream(proci, pBufs)() << remoteCoupledPoints[proci];
            }
        }

        pBufs.finishedSends();

        forAll(remoteCoupledPoints, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                UIPstream(proci, pBufs)() >> remoteCoupledPoints[proci];
            }
        }
    }

    forAll(remoteCoupledPoints, proci)
    {
        forAll(remoteCoupledPoints[proci], i)
        {
            const label coupledPointi = remoteCoupledPoints[proci][i];
            const label pointi = coupledPointPoint_[coupledPointi];

            if (pointIndex[pointi] == -1)
            {
                pointIndex[pointi] = indexPoint.size();
                indexPoint.append(pointi);
            }
        }
    }
}


// ************************************************************************* //
