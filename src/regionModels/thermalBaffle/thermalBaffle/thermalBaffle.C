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
#include "DeltaInfoData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalBaffle, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<volScalarField::Internal> thermalBaffle::calcDelta() const
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
    DynamicList<DeltaInfoData<point>> initialFaceInfoPoints;
    {
        const polyPatch& pp =
            regionMesh().boundaryMesh()[intCoupledPatchIDs_[0]];

        initialFaces.setCapacity(initialFaces.size() + pp.size());
        initialFaceInfoPoints.setCapacity(initialFaces.size() + pp.size());

        forAll(pp, ppFacei)
        {
            const point& c = pp.faceCentres()[ppFacei];

            initialFaces.append(pp.start() + ppFacei);
            initialFaceInfoPoints.append(DeltaInfoData<point>(-1, c));
        }
    }

    // Wave across the mesh layers
    List<DeltaInfoData<point>> faceInfoPoints(regionMesh().nFaces());
    List<DeltaInfoData<point>> cellInfoPoints(regionMesh().nCells());
    FaceCellWave<DeltaInfoData<point>>
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
    DynamicList<DeltaInfoData<scalar>> initialFaceInfoDeltas;
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
                initialFaceInfoDeltas.append(DeltaInfoData<scalar>(-1, d));
            }
        }
    }

    // Wave back across the layers
    List<DeltaInfoData<scalar>> faceInfoDeltas(regionMesh().nFaces());
    List<DeltaInfoData<scalar>> cellInfoDeltas(regionMesh().nCells());
    FaceCellWave<DeltaInfoData<scalar>>
    (
        regionMesh(),
        initialFaces,
        initialFaceInfoDeltas,
        faceInfoDeltas,
        cellInfoDeltas,
        regionMesh().globalData().nTotalCells() + 1
    );

    // Unpack distances into a dimensioned field and return
    tmp<volScalarField::Internal> tDelta
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
    volScalarField::Internal& delta = tDelta.ref();

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

        delta[celli] = cellInfoDeltas[celli].data();
    }

    return tDelta;
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

    // Modify thermo density and diffusivity to take into account the thickness
    volScalarField dByT
    (
        volScalarField::New
        (
            "dByT",
            regionMesh(),
            dimless,
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        )
    );
    dByT.ref() = delta_/thickness_;
    dByT.correctBoundaryConditions();
    const volScalarField rho("rho", thermo_->rho()*dByT);
    const volScalarField alphahe("alphahe", thermo_->alphahe()*dByT);

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, he_)
      - fvm::laplacian(alphahe, he_)
     ==
        Q_ + Qs_/thickness_
    );

    if (regionMesh().moving())
    {
        surfaceScalarField phiMesh
        (
            fvc::interpolate(rho*he_)*regionMesh().phi()
        );

        hEqn -= fvc::div(phiMesh);
    }

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
    delta_(calcDelta()),
    thickness_
    (
        IOobject
        (
            "thickness",
            regionMesh().pointsInstance(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        delta_
    ),
    nNonOrthCorr_(solution().lookup<label>("nNonOrthCorr")),
    thermo_(solidThermo::New(regionMesh())),
    he_(thermo_->he()),
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


const tmp<volScalarField> thermalBaffle::Cp() const
{
    return thermo_->Cp();
}


const volScalarField& thermalBaffle::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


const volScalarField& thermalBaffle::rho() const
{
    return thermo_->rho();
}


const volScalarField& thermalBaffle::kappa() const
{
    return thermo_->kappa();
}


const volScalarField& thermalBaffle::T() const
{
    return thermo_->T();
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
        const fvPatchScalarField& phe = he_.boundaryField()[patchi];
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
