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

#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "SubField.H"
#include "cyclicFvPatchFields.H"
#include "cyclicAMIFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMesh::makeSf() const
{
    if (debug)
    {
        InfoInFunction << "Assembling face areas" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (SfPtr_)
    {
        FatalErrorInFunction
            << "face areas already exist"
            << abort(FatalError);
    }

    SfPtr_ = new slicedSurfaceVectorField
    (
        IOobject
        (
            "Sf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimArea,
        faceAreas()
    );
}


void Foam::fvMesh::makeMagSf() const
{
    if (debug)
    {
        InfoInFunction << "Assembling mag face areas" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (magSfPtr_)
    {
        FatalErrorInFunction
            << "mag face areas already exist"
            << abort(FatalError);
    }

    magSfPtr_ = new slicedSurfaceScalarField
    (
        IOobject
        (
            "magSf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimArea,
        magFaceAreas()
    );
}


void Foam::fvMesh::makeC() const
{
    if (debug)
    {
        InfoInFunction << "Assembling cell centres" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (CPtr_)
    {
        FatalErrorInFunction
            << "cell centres already exist"
            << abort(FatalError);
    }

    // Construct as slices. Only preserve processor (not e.g. cyclic)

    CPtr_ = new slicedVolVectorField
    (
        IOobject
        (
            "C",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimLength,
        cellCentres(),
        faceCentres(),
        true,               // preserveCouples
        true                // preserveProcOnly
    );
}


void Foam::fvMesh::makeCf() const
{
    if (debug)
    {
        InfoInFunction << "Assembling face centres" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (CfPtr_)
    {
        FatalErrorInFunction
            << "face centres already exist"
            << abort(FatalError);
    }

    CfPtr_ = new slicedSurfaceVectorField
    (
        IOobject
        (
            "Cf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimLength,
        faceCentres()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField::Internal& Foam::fvMesh::V() const
{
    if (!VPtr_)
    {
        if (debug)
        {
            InfoInFunction
                << "Constructing from primitiveMesh::cellVolumes()" << endl;
        }

        VPtr_ = new slicedVolScalarField::Internal
        (
            IOobject
            (
                "V",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimVolume,
            cellVolumes()
        );
    }

    return *static_cast<slicedVolScalarField::Internal*>(VPtr_);
}


const Foam::volScalarField::Internal& Foam::fvMesh::V0() const
{
    if (!V0Ptr_)
    {
        FatalErrorInFunction
            << "V0 is not available"
            << abort(FatalError);
    }

    return *V0Ptr_;
}


Foam::volScalarField::Internal& Foam::fvMesh::setV0()
{
    if (!V0Ptr_)
    {
        FatalErrorInFunction
            << "V0 is not available"
            << abort(FatalError);
    }

    return *V0Ptr_;
}


const Foam::volScalarField::Internal& Foam::fvMesh::V00() const
{
    if (!V00Ptr_)
    {
        if (debug)
        {
            InfoInFunction << "Constructing from V0" << endl;
        }

        V00Ptr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "V00",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            V0()
        );
    }

    return *V00Ptr_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fvMesh::Vsc() const
{
    if (moving() && time().subCycling())
    {
        const TimeState& ts = time();
        const TimeState& ts0 = time().prevTimeState();

        scalar tFrac =
        (
            ts.value() - (ts0.value() - ts0.deltaTValue())
        )/ts0.deltaTValue();

        if (tFrac < (1 - small))
        {
            return V0() + tFrac*(V() - V0());
        }
        else
        {
            return V();
        }
    }
    else
    {
        return V();
    }
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fvMesh::Vsc0() const
{
    if (moving() && time().subCycling())
    {
        const TimeState& ts = time();
        const TimeState& ts0 = time().prevTimeState();

        scalar t0Frac =
        (
            (ts.value() - ts.deltaTValue())
          - (ts0.value() - ts0.deltaTValue())
        )/ts0.deltaTValue();

        if (t0Frac > small)
        {
            return V0() + t0Frac*(V() - V0());
        }
        else
        {
            return V0();
        }
    }
    else
    {
        return V0();
    }
}


const Foam::surfaceVectorField& Foam::fvMesh::Sf() const
{
    if (!SfPtr_)
    {
        makeSf();
    }

    return *SfPtr_;
}


const Foam::surfaceScalarField& Foam::fvMesh::magSf() const
{
    if (!magSfPtr_)
    {
        makeMagSf();
    }

    return *magSfPtr_;
}


const Foam::volVectorField& Foam::fvMesh::C() const
{
    if (!CPtr_)
    {
        makeC();
    }

    return *CPtr_;
}


const Foam::surfaceVectorField& Foam::fvMesh::Cf() const
{
    if (!CfPtr_)
    {
        makeCf();
    }

    return *CfPtr_;
}


Foam::tmp<Foam::surfaceVectorField> Foam::fvMesh::delta() const
{
    if (debug)
    {
        InfoInFunction << "Calculating face deltas" << endl;
    }

    tmp<surfaceVectorField> tdelta
    (
        surfaceVectorField::New
        (
            "delta",
            *this,
            dimLength
        )
    );
    surfaceVectorField& delta = tdelta.ref();

    const volVectorField& C = this->C();
    const labelUList& owner = this->owner();
    const labelUList& neighbour = this->neighbour();

    forAll(owner, facei)
    {
        delta[facei] = C[neighbour[facei]] - C[owner[facei]];
    }

    surfaceVectorField::Boundary& deltabf =
        delta.boundaryFieldRef();

    forAll(deltabf, patchi)
    {
        deltabf[patchi] = boundary()[patchi].delta();
    }

    return tdelta;
}


const Foam::surfaceScalarField& Foam::fvMesh::phi() const
{
    if (!phiPtr_)
    {
        FatalErrorInFunction
            << "mesh flux field does not exist, is the mesh actually moving?"
            << abort(FatalError);
    }

    // Set zero current time
    // mesh motion fluxes if the time has been incremented
    if (!time().subCycling() && phiPtr_->timeIndex() != time().timeIndex())
    {
        (*phiPtr_) = dimensionedScalar(dimVolume/dimTime, 0);
    }

    return *phiPtr_;
}


Foam::surfaceScalarField& Foam::fvMesh::setPhi()
{
    if (!phiPtr_)
    {
        FatalErrorInFunction
            << "mesh flux field does not exist, is the mesh actually moving?"
            << abort(FatalError);
    }

    return *phiPtr_;
}


// ************************************************************************* //
