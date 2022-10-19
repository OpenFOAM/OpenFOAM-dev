/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "thermalBaffle.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "absorptionEmissionModel.H"
#include "zeroGradientFvPatchFields.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "FaceCellWave.H"
#include "LayerInfoData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalBaffle, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<volScalarField::Internal> thermalBaffle::calcThickness() const
{
    if (intCoupledPatchIDs_.size() != 2)
    {
        FatalErrorInFunction
            << "Mesh \"" << regionMesh().name()
            << "\" does not have exactly two coupled patches"
            << exit(FatalError);
    }

    // Initial faces from which to wave
    DynamicList<label> initialFaces;

    // Initialise faces on the first coupled patch with their centres as data
    initialFaces.clear();
    DynamicList<LayerInfoData<point>> initialFaceInfoPoints;
    {
        const polyPatch& pp =
            regionMesh().boundaryMesh()[intCoupledPatchIDs_[0]];

        initialFaces.setCapacity(initialFaces.size() + pp.size());
        initialFaceInfoPoints.setCapacity(initialFaces.size() + pp.size());

        forAll(pp, ppFacei)
        {
            const point& c = pp.faceCentres()[ppFacei];

            initialFaces.append(pp.start() + ppFacei);
            initialFaceInfoPoints.append(LayerInfoData<point>(0, -1, c));
        }
    }

    // Wave across the mesh layers
    List<LayerInfoData<point>> faceInfoPoints(regionMesh().nFaces());
    List<LayerInfoData<point>> cellInfoPoints(regionMesh().nCells());
    FaceCellWave<LayerInfoData<point>>
    (
        regionMesh(),
        initialFaces,
        initialFaceInfoPoints,
        faceInfoPoints,
        cellInfoPoints,
        regionMesh().globalData().nTotalCells() + 1
    );

    // Calculate the distances between the opposite patch and load into data to
    // wave back in the opposite direction
    initialFaces.clear();
    DynamicList<LayerInfoData<scalar>> initialFaceInfoDeltas;
    {
        const polyPatch& pp =
            regionMesh().boundaryMesh()[intCoupledPatchIDs_[1]];

        forAll(pp, ppFacei)
        {
            const point& c = pp.faceCentres()[ppFacei];

            static nil td;

            if (faceInfoPoints[pp.start() + ppFacei].valid(td))
            {
                const scalar d =
                    mag(c - faceInfoPoints[pp.start() + ppFacei].data());

                initialFaces.append(pp.start() + ppFacei);
                initialFaceInfoDeltas.append(LayerInfoData<scalar>(0, -1, d));
            }
        }
    }

    // Wave back across the layers
    List<LayerInfoData<scalar>> faceInfoDeltas(regionMesh().nFaces());
    List<LayerInfoData<scalar>> cellInfoDeltas(regionMesh().nCells());
    FaceCellWave<LayerInfoData<scalar>>
    (
        regionMesh(),
        initialFaces,
        initialFaceInfoDeltas,
        faceInfoDeltas,
        cellInfoDeltas,
        regionMesh().globalData().nTotalCells() + 1
    );

    // Unpack distances into a dimensioned field and return
    tmp<volScalarField::Internal> tThickness
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "thickness",
                regionMesh().pointsInstance(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            regionMesh(),
            dimLength
        )
    );
    volScalarField::Internal& thickness = tThickness.ref();

    forAll(cellInfoDeltas, celli)
    {
        static nil td;

        if (!cellInfoDeltas[celli].valid(td))
        {
            FatalErrorInFunction
                << "Mesh \"" << regionMesh().name()
                << "\" is not layered between its coupled patches"
                << exit(FatalError);
        }

        thickness[celli] = cellInfoDeltas[celli].data();
    }

    return tThickness;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermalBaffle::read()
{
    solution().lookup("nNonOrthCorr") >> nNonOrthCorr_;
    return regionModel::read();
}


bool thermalBaffle::read(const dictionary& dict)
{
    solution().lookup("nNonOrthCorr") >> nNonOrthCorr_;
    return regionModel::read(dict);
}


void thermalBaffle::solveEnergy()
{
    DebugInFunction << endl;

    volScalarField& e = thermo_->he();
    const volScalarField& rho = thermo_->rho();

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, e)
      + thermophysicalTransport_->divq(e)
     ==
        Q_ + Qs_/thickness_
    );

    hEqn.relax();
    hEqn.solve();

    thermo_->correct();

    Info<< "T min/max   = " << min(thermo_->T()) << ", "
        << max(thermo_->T()) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalBaffle::thermalBaffle
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    regionModel(mesh, "thermalBaffle", modelType, dict, true),
    thickness_(calcThickness()),
    nNonOrthCorr_(solution().lookup<label>("nNonOrthCorr")),
    thermo_(solidThermo::New(regionMesh())),
    thermophysicalTransport_(solidThermophysicalTransportModel::New(thermo())),
    Qs_
    (
        IOobject
        (
            "Qs",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimArea/dimTime, Zero)
    ),
    Q_
    (
        IOobject
        (
            "Q",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    ),
    radiation_
    (
        radiationModel::New
        (
            dict.subDict("radiation"),
            thermo_->T()
        )
    )
{
    thermo_->correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermalBaffle::~thermalBaffle()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermalBaffle::preEvolveRegion()
{}


void thermalBaffle::evolveRegion()
{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveEnergy();
    }
}


const solidThermo& thermalBaffle::thermo() const
{
    return thermo_;
}


void thermalBaffle::info()
{
    const labelList& coupledPatches = intCoupledPatchIDs();

    forAll(coupledPatches, i)
    {
        const label patchi = coupledPatches[i];
        const fvPatchScalarField& phe = thermo_->he().boundaryField()[patchi];
        const word patchName = regionMesh().boundary()[patchi].name();

        Info<< indent << "Q : " << patchName << indent
            <<
            gSum
            (
                mag(regionMesh().Sf().boundaryField()[patchi])
               *phe.snGrad()
               *thermo_->alphahe(patchi)
            )
            << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace regionModels
} // end namespace Foam

// ************************************************************************* //
