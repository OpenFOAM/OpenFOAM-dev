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

#include "volumeInjection.H"
#include "CompactListList.H"
#include "LagrangianFields.H"
#include "tetIndices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(volumeInjection, 0);
    addToRunTimeSelectionTable(LagrangianModel, volumeInjection, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Lagrangian::volumeInjection::readCoeffs(const dictionary& modelDict)
{
    set_.read(modelDict);

    haveNumber_ = modelDict.found("number");
    const bool haveNumberDensity = modelDict.found("numberDensity");

    if (haveNumber_ == haveNumberDensity)
    {
        FatalIOErrorInFunction(modelDict)
            << "keywords number and numberDensity are both "
            << (haveNumber_ ? "" : "un") << "defined in "
            << "dictionary " << modelDict.name()
            << exit(FatalIOError);
    }

    numberOrNumberDensity_ =
        haveNumber_
      ? modelDict.lookup<scalar>("number", dimless)
      : modelDict.lookup<scalar>("numberDensity", dimless/dimVolume);

    time_ =
        modelDict.lookupOrDefault<scalar>
        (
            "time",
            mesh().time().userUnits(),
            mesh().time().beginTime().value()
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::volumeInjection::volumeInjection
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianInjection(name, mesh),
    set_(mesh.mesh()),
    haveNumber_(false),
    numberOrNumberDensity_(NaN),
    time_(NaN),
    globalRndGen_("globalRndGen", stateDict, name, true),
    localRndGen_("localRndGen", stateDict, name, false),
    timeIndex_(-1)
{
    readCoeffs(modelDict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::Lagrangian::volumeInjection::correct()
{
    // In theory the set can change as a result of motion. We can't handle this
    // entirely correctly as we can't construct the set at an intermediate
    // time. We assume that the most likely scenario is injection at the
    // beginning of the time-step. So, when we are inside the injection
    // time-step we delay update of the set until after injection.

    const scalar t1 = mesh().time().value();
    const scalar t0 = t1 - mesh().time().deltaT().value();

    if (t0 <= time_ && time_ < t1) return;

    if (mesh().mesh().moving())
    {
        set_.movePoints();
    }
}


Foam::LagrangianSubMesh Foam::Lagrangian::volumeInjection::modify
(
    LagrangianMesh& mesh,
    const LagrangianSubMesh&
) const
{
    const scalar t1 = mesh.time().value();
    const scalar t0 = t1 - mesh.time().deltaT().value();

    if (!(t0 <= time_ && time_ < t1)) return mesh.subNone();

    const scalar fraction = (time_ - t0)/(t1 - t0);

    // Restart the generators if necessary and set the time index up to date
    localRndGen_.start(timeIndex_ == db().time().timeIndex());
    globalRndGen_.start(timeIndex_ == db().time().timeIndex());
    timeIndex_ = db().time().timeIndex();

    // Reference the mesh faces and the set cells
    const faceList& faces = mesh.mesh().faces();
    const labelUList setCellCells = set_.cells();
    const UIndirectList<cell> setCells(mesh.mesh().cells(), setCellCells);

    // Create point and cell-centre fields at the current fraction
    const pointField points
    (
        (1 - fraction)*mesh.mesh().oldPoints()
      + fraction*mesh.mesh().points()
    );
    const pointField cellCentres
    (
        (1 - fraction)*mesh.mesh().oldCellCentres()
      + fraction*mesh.mesh().cellCentres()
    );

    // Count the number of tetrahedra in each set cell
    labelList setCellNTets(setCells.size(), label(0));
    forAll(setCells, setCelli)
    {
        forAll(setCells[setCelli], cellFacei)
        {
            setCellNTets[setCelli] +=
                faces[setCells[setCelli][cellFacei]].nTriangles();
        }
    }

    // Construct indexing for the cell-tetrahedra and cumulative sums of the
    // volumes in the cells and cell-tetrahedra
    CompactListList<labelPair> setCellTetFaceAndFaceTri
    (
        setCellNTets,
        labelPair(-1, -1)
    );
    scalarList setCellSumVolume(setCells.size(), scalar(0));
    CompactListList<scalar> setCellTetSumVolume(setCellNTets, scalar(0));
    forAll(setCells, setCelli)
    {
        const cell& c = setCells[setCelli];

        label cellTeti = 0;

        forAll(c, cFacei)
        {
            const face& f = faces[c[cFacei]];

            for (label fTrii = 1; fTrii <= f.nTriangles(); ++ fTrii)
            {
                setCellTetFaceAndFaceTri[setCelli][cellTeti] =
                    labelPair(c[cFacei], fTrii);

                setCellSumVolume[setCelli] +=
                    tetIndices(setCellCells[setCelli], c[cFacei], fTrii)
                   .tet(mesh.mesh(), points, cellCentres)
                   .mag();

                setCellTetSumVolume[setCelli][cellTeti] =
                    setCellSumVolume[setCelli];

                cellTeti ++;
            }
        }
    }
    for (label setCelli = 1; setCelli < setCells.size(); ++ setCelli)
    {
        setCellSumVolume[setCelli] += setCellSumVolume[setCelli - 1];
    }

    // Construct a cumulative sum of the volumes across the processes
    scalarList procSumVolume(Pstream::nProcs(), -vGreat);
    procSumVolume[Pstream::myProcNo()] = setCellSumVolume.last();
    Pstream::listCombineGather(procSumVolume, maxEqOp<scalar>());
    Pstream::listCombineScatter(procSumVolume);
    for (label proci = 1; proci < Pstream::nProcs(); proci ++)
    {
        procSumVolume[proci] += procSumVolume[proci - 1];
    }

    // Determine the number of particles to inject
    const label numberInt =
        haveNumber_
      ? round(numberOrNumberDensity_)
      : round(numberOrNumberDensity_*procSumVolume.last());

    // Binary search for a given volume within a supplied cumulative volume sum.
    // Once the interval has been found, subtract the lower value from that
    // volume so that a subsequent binary search can be done on a sub-set of
    // volumes. Note that the supplied cumulative volume sum does not include
    // the leading zero.
    auto findVolume = []
    (
        const scalarUList& volumes,
        scalar& volume
    )
    {
        label i0 = -1, i1 = volumes.size() - 1;

        while (i0 + 1 != i1)
        {
            const label i = (i0 + i1)/2;

            (volume < volumes[i] ? i1 : i0) = i;
        }

        if (i0 != -1)
        {
            volume -= volumes[i0];
        }

        return i1;
    };

    // Initialise storage for the injection geometry and topology. This is
    // dynamic as we don't know how much will end up on each processor yet.
    DynamicList<barycentric> injectCoordinates(numberInt/Pstream::nProcs());
    DynamicList<label> injectCells(numberInt/Pstream::nProcs());
    DynamicList<label> injectFaces(numberInt/Pstream::nProcs());
    DynamicList<label> injectFaceTris(numberInt/Pstream::nProcs());

    // Create a (global) list of volumes at which to inject. Each volume is
    // searched for in the cumulative lists to identify a tetrahedron into
    // which to inject.
    scalarField volume(globalRndGen_.scalar01(numberInt)*procSumVolume.last());
    forAll(volume, volumei)
    {
        const label proci = findVolume(procSumVolume, volume[volumei]);

        if (Pstream::myProcNo() == proci)
        {
            const label setCelli =
                findVolume(setCellSumVolume, volume[volumei]);
            const label celli = setCellCells[setCelli];
            const label cellTeti =
                findVolume(setCellTetSumVolume[setCelli], volume[volumei]);

            injectCoordinates.append(barycentric01(localRndGen_));
            injectCells.append(celli);
            injectFaces.append
            (
                setCellTetFaceAndFaceTri[setCelli][cellTeti].first()
            );
            injectFaceTris.append
            (
                setCellTetFaceAndFaceTri[setCelli][cellTeti].second()
            );
        }
    }

    // Create the particles at the identified locations
    LagrangianSubMesh injectionMesh =
        mesh.inject
        (
            *this,
            barycentricField(injectCoordinates),
            labelField(injectCells),
            labelField(injectFaces),
            labelField(injectFaceTris),
            LagrangianMesh::fractionName,
            scalarField(injectCoordinates.size(), fraction)
        );

    // See the note in volumeInjection::correct
    if (mesh.mesh().moving())
    {
        set_.movePoints();
    }

    return injectionMesh;
}


void Foam::Lagrangian::volumeInjection::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::Lagrangian::volumeInjection::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::Lagrangian::volumeInjection::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
}


bool Foam::Lagrangian::volumeInjection::read(const dictionary& modelDict)
{
    if (LagrangianInjection::read(modelDict))
    {
        readCoeffs(modelDict);
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::Lagrangian::volumeInjection::writeState(Ostream& os) const
{
    LagrangianInjection::writeState(os);

    writeEntry(os, "globalRndGen", globalRndGen_);
}


void Foam::Lagrangian::volumeInjection::writeProcessorState(Ostream& os) const
{
    LagrangianInjection::writeProcessorState(os);

    writeEntry(os, "localRndGen", localRndGen_);
}


// ************************************************************************* //
