/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvMatrices.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace thermalBaffleModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalBaffle, 0);

addToRunTimeSelectionTable(thermalBaffleModel, thermalBaffle, mesh);
addToRunTimeSelectionTable(thermalBaffleModel, thermalBaffle, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermalBaffle::read()
{
    this->solution().lookup("nNonOrthCorr") >> nNonOrthCorr_;
    return regionModel1D::read();
}


bool thermalBaffle::read(const dictionary& dict)
{
    this->solution().lookup("nNonOrthCorr") >> nNonOrthCorr_;
    return regionModel1D::read(dict);
}


void thermalBaffle::solveEnergy()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();

    tmp<volScalarField> tQ
    (
        new volScalarField
        (
            IOobject
            (
                "tQ",
                regionMesh().time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            regionMesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    volScalarField& Q = tQ.ref();

    volScalarField rho("rho", thermo_->rho());
    volScalarField alpha("alpha", thermo_->alpha());


    // If region is one-dimension variable thickness
    if (oneD_ && !constantThickness_)
    {
        // Scale K and rhoCp and fill Q in the internal baffle region.
        const label patchi = intCoupledPatchIDs_[0];
        const polyPatch& ppCoupled = rbm[patchi];

        forAll(ppCoupled, localFacei)
        {
            const labelList& cells = boundaryFaceCells_[localFacei];
            forAll(cells, i)
            {
                const label cellId = cells[i];

                Q[cellId] =
                    Qs_.boundaryField()[patchi][localFacei]
                   /thickness_[localFacei];

                rho[cellId] *= delta_.value()/thickness_[localFacei];

                alpha[cellId] *= delta_.value()/thickness_[localFacei];
            }
        }
    }
    else
    {
        Q = Q_;
    }

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h_)
      - fvm::laplacian(alpha, h_)
     ==
        Q
    );

    if (moveMesh_)
    {
        surfaceScalarField phiMesh
        (
            fvc::interpolate(rho*h_)*regionMesh().phi()
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
    thermalBaffleModel(modelType, mesh, dict),
    nNonOrthCorr_(readLabel(solution().lookup("nNonOrthCorr"))),
    thermo_(solidThermo::New(regionMesh(), dict)),
    h_(thermo_->he()),
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
        dimensionedScalar
        (
            "zero",
            dimEnergy/dimArea/dimTime,
            Zero
        )
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
        dimensionedScalar
        (
            "zero",
            dimEnergy/dimVolume/dimTime,
            Zero
        )
    ),
    radiation_
    (
        radiation::radiationModel::New
        (
            dict.subDict("radiation"),
            thermo_->T()
        )
    )
{
    init();
    thermo_->correct();
}


thermalBaffle::thermalBaffle
(
    const word& modelType,
    const fvMesh& mesh
)
:
    thermalBaffleModel(modelType, mesh),
    nNonOrthCorr_(readLabel(solution().lookup("nNonOrthCorr"))),
    thermo_(solidThermo::New(regionMesh())),
    h_(thermo_->he()),
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
        dimensionedScalar
        (
            "zero",
            dimEnergy/dimArea/dimTime,
            Zero
        )
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
        dimensionedScalar
        (
            "zero",
            dimEnergy/dimVolume/dimTime,
            Zero
        )
    ),
    radiation_
    (
        radiation::radiationModel::New
        (
            thermo_->T()
        )
    )
{
    init();
    thermo_->correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermalBaffle::~thermalBaffle()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermalBaffle::init()
{
    if (oneD_ && !constantThickness_)
    {
        label patchi = intCoupledPatchIDs_[0];
        const label Qsb = Qs_.boundaryField()[patchi].size();

        if (Qsb!= thickness_.size())
        {
            FatalErrorInFunction
                << "the boundary field of Qs is "
                << Qsb << " and " << nl
                << "the field 'thickness' is " << thickness_.size() << nl
                << exit(FatalError);
        }
    }
}


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
        const fvPatchScalarField& ph = h_.boundaryField()[patchi];
        const word patchName = regionMesh().boundary()[patchi].name();
        Info<< indent << "Q : " << patchName << indent <<
            gSum
            (
                mag(regionMesh().Sf().boundaryField()[patchi])
              * ph.snGrad()
              * thermo_->alpha().boundaryField()[patchi]
            ) << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace thermalBaffleModels
} // end namespace regionModels
} // end namespace Foam

// ************************************************************************* //
