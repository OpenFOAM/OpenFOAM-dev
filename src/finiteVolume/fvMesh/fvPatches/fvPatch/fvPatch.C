/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "fvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvBoundaryMesh.H"
#include "fvMesh.H"
#include "primitiveMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvPatch, 0);
    defineRunTimeSelectionTable(fvPatch, polyPatch);
    addToRunTimeSelectionTable(fvPatch, fvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvPatch::fvPatch(const polyPatch& p, const fvBoundaryMesh& bm)
:
    poly_(p),
    boundaryMesh_(bm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvPatch::~fvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::objectRegistry& Foam::fvPatch::db() const
{
    return boundaryMesh_.mesh();
}


const Foam::fvMesh& Foam::fvPatch::mesh() const
{
    return boundaryMesh_.mesh();
}


const Foam::Time& Foam::fvPatch::time() const
{
    return boundaryMesh_.mesh().time();
}


const Foam::labelUList& Foam::fvPatch::faceCells() const
{
    return poly_.faceCells();
}


const Foam::vectorField& Foam::fvPatch::Cf() const
{
    return mesh().Cf().boundaryField()[index()];
}


const Foam::DimensionedField<Foam::vector, Foam::fvPatch>&
Foam::fvPatch::C() const
{
    if (!CPtr_.valid())
    {
        CPtr_ = new SlicedDimensionedField<vector, fvPatch>
        (
            IOobject
            (
                "C",
                mesh().time().name(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimensions::length,
            Cf()
        );
    }
    else
    {
        CPtr_->reset(Cf());
    }

    return *CPtr_;
}



Foam::tmp<Foam::vectorField> Foam::fvPatch::Cn() const
{
    tmp<vectorField> tcc(new vectorField(size()));
    vectorField& cc = tcc.ref();

    const labelUList& faceCells = this->faceCells();

    // get reference to global cell centres
    const vectorField& gcc = mesh().cellCentres();

    forAll(faceCells, facei)
    {
        cc[facei] = gcc[faceCells[facei]];
    }

    return tcc;
}


Foam::tmp<Foam::vectorField> Foam::fvPatch::nf() const
{
    return Sf()/magSf();
}


const Foam::vectorField& Foam::fvPatch::Sf() const
{
    return mesh().Sf().boundaryField()[index()];
}


const Foam::scalarField& Foam::fvPatch::magSf() const
{
    return mesh().magSf().boundaryField()[index()];
}


Foam::tmp<Foam::vectorField> Foam::fvPatch::delta() const
{
    // Use patch-normal delta for all non-coupled BCs
    const vectorField nHat(nf());
    return nHat*(nHat & (Cf() - Cn()));
}


Foam::tmp<Foam::scalarField> Foam::fvPatch::polyFaceFraction() const
{
    return
        mesh().conformal()
      ? tmp<scalarField>(new scalarField(size(), scalar(1)))
      : magSf()
       /scalarField
        (
            mesh().magFaceAreas(),
            mesh().polyFacesBf()[index()]
        );
}


void Foam::fvPatch::makeWeights(scalarField& w) const
{
    w = 1.0;
}


const Foam::scalarField& Foam::fvPatch::deltaCoeffs() const
{
    return mesh().deltaCoeffs().boundaryField()[index()];
}


const Foam::scalarField& Foam::fvPatch::weights() const
{
    return mesh().weights().boundaryField()[index()];
}


const Foam::fvMesh& Foam::fvPatch::operator()() const
{
    return mesh();
}


// ************************************************************************* //
