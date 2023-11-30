/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "anisotropic.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SolidThermophysicalTransportModel>
void Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::
setZonesPatchFaces() const
{
    if (!zoneCoordinateSystems_.size()) return;

    // Find all the patch faces adjacent to zones

    const fvMesh& mesh = this->thermo().mesh();
    const fvBoundaryMesh& patches = mesh.boundary();
    const labelList& own = mesh.faceOwner();

    zonesPatchFaces_.setSize(zoneCoordinateSystems_.size());
    label zonei = 0;

    forAllConstIter
    (
        PtrDictionary<coordinateSystem>,
        zoneCoordinateSystems_,
        iter
    )
    {
        zonesPatchFaces_[zonei].setSize(patches.size());

        const labelList& zoneCells = mesh.cellZones()[iter().name()];

        // Cell in zone
        boolList cellInZone(mesh.nCells(), false);

        forAll(zoneCells, i)
        {
            cellInZone[zoneCells[i]] = true;
        }

        forAll(patches, patchi)
        {
            const fvPatch& pp = patches[patchi];

            zonesPatchFaces_[zonei][patchi].setSize(pp.size());

            label nZonePatchFaces = 0;

            forAll(pp, patchFacei)
            {
                const label facei = pp.start() + patchFacei;

                if (cellInZone[own[facei]])
                {
                    zonesPatchFaces_[zonei][patchi][nZonePatchFaces++] =
                        patchFacei;
                }
            }

            zonesPatchFaces_[zonei][patchi].setSize(nZonePatchFaces);
        }

        zonei++;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SolidThermophysicalTransportModel>
Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::anisotropic
(
    const alphaField& alpha,
    const solidThermo& thermo
)
:
    SolidThermophysicalTransportModel(typeName, alpha, thermo),
    TopoChangeableMeshObject(*this, thermo.mesh()),
    coordinateSystem_(coordinateSystem::New(thermo.mesh(), this->coeffDict())),
    boundaryAligned_
    (
        this->coeffDict().template lookupOrDefault<Switch>
        (
            "boundaryAligned",
            false
        )
    ),
    aligned_(thermo.mesh().boundary().size(), true)
{
    if (this->coeffDict().found("zones"))
    {
        const dictionary& zonesDict(this->coeffDict().subDict("zones"));

        Info<< "    Reading coordinate system for zones:" << endl;

        forAllConstIter(dictionary, zonesDict, iter)
        {
            if (iter().isDict())
            {
                const word& zoneName = iter().keyword();
                const dictionary& dict = iter().dict();

                Info<< "        " << zoneName << endl;

                zoneCoordinateSystems_.insert
                (
                    zoneName,
                    coordinateSystem::New(zoneName, dict).ptr()
                );
            }
        }

        // Find all the patch faces adjacent to zones
        setZonesPatchFaces();
    }

    const fvMesh& mesh = thermo.mesh();
    const fvBoundaryMesh& bMesh = mesh.boundary();

    bool aligned = true;

    forAll(bMesh, patchi)
    {
        if
        (
            !bMesh[patchi].coupled()
          && returnReduce(bMesh[patchi].size(), sumOp<label>())
        )
        {
            const vectorField n(bMesh[patchi].nf());
            const vectorField nKappa(n & Kappa(patchi));

            // Calculate a patch alignment factor for Kappa
            const scalar alignmentFactor =
                gSum
                (
                    mesh.magSf().boundaryField()[patchi]
                   *mag(nKappa - n*(nKappa & n))
                )/gSum(mesh.magSf().boundaryField()[patchi]*mag(nKappa));

            if (alignmentFactor > 1e-3)
            {
                aligned_[patchi] = false;
                aligned = false;

                Info<< "    Kappa is not aligned with patch "
                    << bMesh[patchi].name()
                    << " by an alignment factor of " << alignmentFactor
                    << " (0=aligned, 1=unaligned)" << nl
                    << "        heat-flux correction will be applied." << endl;
            }
        }
    }

    if (!aligned && boundaryAligned_)
    {
        aligned_ = true;
        aligned = true;

        Info<<
            "    boundaryAligned is set true, "
            "boundary alignment of kappa will be enforced." << endl;
    }

    // If Kappa is not aligned with any patch enable grad(T) caching
    // because the patch grad(T) will be required for the heat-flux correction
    if (!aligned)
    {
        mesh.solution().enableCache("grad(T)");
    }

    Info << endl;

    mesh.schemes().setFluxRequired(thermo.T().name());

    this->printCoeffs(typeName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SolidThermophysicalTransportModel>
bool Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::read()
{

    return true;
}


template<class SolidThermophysicalTransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::Kappa() const
{
    const solidThermo& thermo = this->thermo();
    const fvMesh& mesh = thermo.mesh();

    const volVectorField& materialKappa = thermo.Kappa();

    tmp<volSymmTensorField> tKappa
    (
        volSymmTensorField::New
        (
            "Kappa",
            mesh,
            dimensionedSymmTensor(materialKappa.dimensions(), Zero)
        )
    );
    volSymmTensorField& Kappa = tKappa.ref();

    Kappa.primitiveFieldRef() =
        coordinateSystem_.R(mesh.C()).transformDiagTensor(materialKappa);

    forAll(Kappa.boundaryField(), patchi)
    {
        Kappa.boundaryFieldRef()[patchi] =
            coordinateSystem_.R(mesh.boundary()[patchi].Cf())
           .transformDiagTensor(materialKappa.boundaryField()[patchi]);
    }

    label zonei = 0;

    forAllConstIter
    (
        PtrDictionary<coordinateSystem>,
        zoneCoordinateSystems_,
        iter
    )
    {
        const labelList& zoneCells = mesh.cellZones()[iter().name()];
        const coordinateSystem& cs = iter();

        symmTensorField& KappaIf = Kappa;

        forAll(zoneCells, i)
        {
            const label celli = zoneCells[i];

            KappaIf[celli] =
                cs.R().transformDiagTensor
                (
                    mesh.C()[celli],
                    materialKappa[celli]
                );
        }

        forAll(zonesPatchFaces_[zonei], patchi)
        {
            symmTensorField& KappaPf = Kappa.boundaryFieldRef()[patchi];

            const vectorField& materialKappaPf =
                materialKappa.boundaryField()[patchi];

            const vectorField& CPf = mesh.boundary()[patchi].Cf();

            forAll(zonesPatchFaces_[zonei][patchi], zonePatchFacei)
            {
                const label patchFacei =
                    zonesPatchFaces_[zonei][patchi][zonePatchFacei];

                KappaPf[patchFacei] =
                    cs.R().transformDiagTensor
                    (
                        CPf[patchFacei],
                        materialKappaPf[patchFacei]
                    );
            }
        }

        zonei++;
    }

    return tKappa;
}


template<class SolidThermophysicalTransportModel>
Foam::tmp<Foam::symmTensorField>
Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::Kappa
(
    const label patchi
) const
{
    const solidThermo& thermo = this->thermo();
    const vectorField& CPf = thermo.mesh().boundary()[patchi].Cf();

    const vectorField& materialKappaPf(thermo.Kappa().boundaryField()[patchi]);

    tmp<symmTensorField> tKappa
    (
        coordinateSystem_.R(CPf).transformDiagTensor(materialKappaPf)
    );
    symmTensorField& KappaPf = tKappa.ref();

    label zonei = 0;

    forAllConstIter
    (
        PtrDictionary<coordinateSystem>,
        zoneCoordinateSystems_,
        iter
    )
    {
        const coordinateSystem& cs = iter();

        forAll(zonesPatchFaces_[zonei][patchi], zonePatchFacei)
        {
            const label patchFacei =
                zonesPatchFaces_[zonei][patchi][zonePatchFacei];

            KappaPf[patchFacei] =
                cs.R().transformDiagTensor
                (
                    CPf[patchFacei],
                    materialKappaPf[patchFacei]
                );
        }

        zonei++;
    }

    return tKappa;
}


template<class SolidThermophysicalTransportModel>
Foam::tmp<Foam::volScalarField>
Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::kappa() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


template<class SolidThermophysicalTransportModel>
Foam::tmp<Foam::scalarField>
Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::kappa
(
    const label patchi
) const
{
    const vectorField n(this->thermo().mesh().boundary()[patchi].nf());
    return n & Kappa(patchi) & n;
}


template<class SolidThermophysicalTransportModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::q() const
{
    const solidThermo& thermo = this->thermo();
    const fvMesh& mesh = thermo.mesh();

    return surfaceScalarField::New
    (
        "q",
       -fvm::laplacian(this->alpha()*Kappa(), thermo.T())().flux()/mesh.magSf()
    );
}


template<class SolidThermophysicalTransportModel>
Foam::tmp<Foam::scalarField>
Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::qCorr
(
    const label patchi
) const
{
    if (!aligned_[patchi])
    {
        tmp<volVectorField> gradT(fvc::grad(this->thermo().T()));

        const vectorField n(this->thermo().mesh().boundary()[patchi].nf());
        const vectorField nKappa(n & Kappa(patchi));

        return
           -(
                this->alpha().boundaryField()[patchi]
               *((nKappa - n*(nKappa & n)) & gradT().boundaryField()[patchi])
            );
    }
    else
    {
        return tmp<scalarField>(nullptr);
    }
}


template<class SolidThermophysicalTransportModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::divq
(
    volScalarField& e
) const
{
    const solidThermo& thermo = this->thermo();
    const volSymmTensorField Kappa(this->Kappa());
    const surfaceVectorField& Sf = thermo.mesh().Sf();
    const surfaceScalarField& magSf = thermo.mesh().magSf();

    // Return heat flux source as an implicit energy correction
    // to the temperature gradient flux
    return
       -fvc::laplacian(this->alpha()*Kappa, thermo.T())
       -fvm::laplacianCorrection
        (
            (Sf & fvc::interpolate(this->alpha()*Kappa/thermo.Cv()) & Sf)
           /sqr(magSf),
            e
        );
}


template<class SolidThermophysicalTransportModel>
void Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::predict()
{
    SolidThermophysicalTransportModel::predict();

    // Recalculate zonesPatchFaces if they have been deleted
    // following mesh changes
    if (!zonesPatchFaces_.size())
    {
        setZonesPatchFaces();
    }
}


template<class SolidThermophysicalTransportModel>
bool Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::movePoints()
{
    return true;
}


template<class SolidThermophysicalTransportModel>
void Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::topoChange
(
    const polyTopoChangeMap& map
)
{
    // Delete the cached zonesPatchFaces, will be re-created in predict
    zonesPatchFaces_.clear();
}


template<class SolidThermophysicalTransportModel>
void Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::mapMesh
(
    const polyMeshMap& map
)
{
    // Delete the cached zonesPatchFaces, will be re-created in predict
    zonesPatchFaces_.clear();
}


template<class SolidThermophysicalTransportModel>
void Foam::solidThermophysicalTransportModels::
anisotropic<SolidThermophysicalTransportModel>::distribute
(
    const polyDistributionMap& map
)
{
    // Delete the cached zonesPatchFaces, will be re-created in predict
    zonesPatchFaces_.clear();
}


// ************************************************************************* //
