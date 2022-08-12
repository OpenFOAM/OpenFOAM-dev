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

#include "MRFZone.H"
#include "MRFPatchField.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "faceSet.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MRFZone, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::MRFZone::setMRFFaces()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Is face in MRF zone:
    //     false: not in MRF zone
    //     true: in MRF zone
    boolList faceInMRF(mesh_.nFaces(), false);

    // Determine faces in cell zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (without constructing cells)

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Cells in zone
    boolList zoneCell(mesh_.nCells(), false);

    const labelList& cellLabels = cellSet_.cells();
    forAll(cellLabels, i)
    {
        zoneCell[cellLabels[i]] = true;
    }

    label nZoneFaces = 0;

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (zoneCell[own[facei]] || zoneCell[nei[facei]])
        {
            faceInMRF[facei] = true;
            nZoneFaces++;
        }
    }

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                const label facei = pp.start() + i;

                if (zoneCell[own[facei]])
                {
                    faceInMRF[facei] = true;
                    nZoneFaces++;
                }
            }
        }
    }

    // Synchronise the faceInMRF across processor patches
    syncTools::syncFaceList(mesh_, faceInMRF, orEqOp<bool>());

    // Sort into lists per patch.

    internalFaces_.setSize(mesh_.nFaces());
    label nInternal = 0;

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceInMRF[facei])
        {
            internalFaces_[nInternal++] = facei;
        }
    }
    internalFaces_.setSize(nInternal);

    labelList nPatchFaces(patches.size(), 0);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, patchFacei)
        {
            const label facei = pp.start() + patchFacei;

            if (faceInMRF[facei])
            {
                nPatchFaces[patchi]++;
            }
        }
    }

    patchFaces_.setSize(patches.size());
    forAll(patchFaces_, patchi)
    {
        patchFaces_[patchi].setSize(nPatchFaces[patchi]);
    }
    nPatchFaces = 0;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, patchFacei)
        {
            const label facei = pp.start() + patchFacei;

            if (faceInMRF[facei])
            {
                patchFaces_[patchi][nPatchFaces[patchi]++] = patchFacei;
            }
        }
    }


    if (debug)
    {
        faceSet internalFaces(mesh_, "internalFaces", internalFaces_);

        Pout<< "Writing " << internalFaces.size()
            << " internal faces in MRF zone to faceSet "
            << internalFaces.name() << endl;

        internalFaces.write();

        faceSet patchFaces(mesh_, "patchFaces", 100);
        forAll(patchFaces_, patchi)
        {
            forAll(patchFaces_[patchi], i)
            {
                const label patchFacei = patchFaces_[patchi][i];
                patchFaces.insert(patches[patchi].start()+patchFacei);
            }
        }

        Pout<< "Writing " << patchFaces.size()
            << " faces in MRF zone with special handling to faceSet "
            << patchFaces.name() << endl;

        patchFaces.write();
    }
}


void Foam::MRFZone::checkMRFBCs(const volVectorField& U) const
{
    static bool checked = false;

    if (!checked && U.nOldTimes())
    {
        const volVectorField::Boundary& Ubf = U.boundaryField();

        forAll(Ubf, patchi)
        {
            if (isA<MRFPatchField>(Ubf[patchi]))
            {
                return;
            }
        }

        FatalErrorInFunction
            << "Field " << U.name()
            << " does not provide any MRF specific boundary conditions "
               "for MRF region " << name() << nl
            << "    Walls rotating in the MRF region should have the "
               "MRFnoSlip boundary condition."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZone::MRFZone
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    name_(name),
    coeffs_(dict),
    cellSet_(mesh, coeffs_),
    origin_(coeffs_.lookup("origin")),
    axis_(coeffs_.lookup("axis")),
    omega_(coeffs_)
{
    axis_ = axis_/mag(axis_);
    setMRFFaces();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::MRFZone::Omega() const
{
    return omega_.value(mesh_.time().userTimeValue())*axis_;
}


void Foam::MRFZone::addCoriolis
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    checkMRFBCs(U);

    const labelList& cells = cellSet_.cells();
    vectorField& ddtUc = ddtU.primitiveFieldRef();
    const vectorField& Uc = U;

    const vector Omega = this->Omega();

    forAll(cells, i)
    {
        const label celli = cells[i];
        ddtUc[celli] += (Omega ^ Uc[celli]);
    }
}


void Foam::MRFZone::addCentrifugalAcceleration
(
    volVectorField& centrifugalAcceleration
) const
{
    const labelList& cells = cellSet_.cells();
    const volVectorField& C = mesh_.C();
    vectorField& cac = centrifugalAcceleration.primitiveFieldRef();

    const vector Omega = this->Omega();

    forAll(cells, i)
    {
        const label celli = cells[i];
        cac[celli] -= Omega ^ (Omega ^ (C[celli] - origin_));
    }

    volVectorField::Boundary& caf = centrifugalAcceleration.boundaryFieldRef();

    forAll(patchFaces_, patchi)
    {
        forAll(patchFaces_[patchi], i)
        {
            const label patchFacei = patchFaces_[patchi][i];
            caf[patchi][patchFacei] -=
                Omega
              ^ (Omega ^ (C.boundaryField()[patchi][patchFacei] - origin_));
        }
    }
}


void Foam::MRFZone::makeRelative(volVectorField& U) const
{
    const volVectorField& C = mesh_.C();
    const labelList& cells = cellSet_.cells();

    const vector Omega = this->Omega();

    forAll(cells, i)
    {
        const label celli = cells[i];
        U[celli] -= (Omega ^ (C[celli] - origin_));
    }

    volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(patchFaces_, patchi)
    {
        forAll(patchFaces_[patchi], i)
        {
            const label patchFacei = patchFaces_[patchi][i];
            Ubf[patchi][patchFacei] -=
                (Omega
              ^ (C.boundaryField()[patchi][patchFacei] - origin_));
        }
    }
}


void Foam::MRFZone::makeRelative(surfaceScalarField& phi) const
{
    makeRelativeRhoFlux(geometricOneField(), phi);
}


void Foam::MRFZone::makeRelative(FieldField<fvsPatchField, scalar>& phi) const
{
    makeRelativeRhoFlux(oneFieldField(), phi);
}


void Foam::MRFZone::makeRelative(Field<scalar>& phi, const label patchi) const
{
    makeRelativeRhoFlux(oneField(), phi, patchi);
}


void Foam::MRFZone::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeRelativeRhoFlux(rho, phi);
}


void Foam::MRFZone::makeAbsolute(volVectorField& U) const
{
    const volVectorField& C = mesh_.C();
    const labelList& cells = cellSet_.cells();

    const vector Omega = this->Omega();

    forAll(cells, i)
    {
        const label celli = cells[i];
        U[celli] += (Omega ^ (C[celli] - origin_));
    }

    volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(patchFaces_, patchi)
    {
        forAll(patchFaces_[patchi], i)
        {
            const label patchFacei = patchFaces_[patchi][i];
            Ubf[patchi][patchFacei] +=
                (Omega ^ (C.boundaryField()[patchi][patchFacei] - origin_));
        }
    }
}


void Foam::MRFZone::makeAbsolute(surfaceScalarField& phi) const
{
    makeAbsoluteRhoFlux(geometricOneField(), phi);
}


void Foam::MRFZone::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeAbsoluteRhoFlux(rho, phi);
}


void Foam::MRFZone::makeAbsolute(Field<vector>& Up, const label patchi) const
{
    const vector Omega = this->Omega();

    Up += (Omega ^ (mesh_.Cf().boundaryField()[patchi] - origin_));
}


bool Foam::MRFZone::read(const dictionary& dict)
{
    coeffs_ = dict;
    cellSet_.read(coeffs_);
    setMRFFaces();

    return true;
}


void Foam::MRFZone::update()
{
    if (mesh_.topoChanged())
    {
        setMRFFaces();
    }
}


// ************************************************************************* //
