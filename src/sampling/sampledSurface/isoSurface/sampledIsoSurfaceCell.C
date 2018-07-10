/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "sampledIsoSurfaceCell.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "isoSurfaceCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledIsoSurfaceCell, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledIsoSurfaceCell,
        word,
        isoSurfaceCell
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::sampledIsoSurfaceCell::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // no update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    prevTimeIndex_ = fvm.time().timeIndex();

    // Clear any stored topo
    facesPtr_.clear();

    // Clear derived data
    sampledSurface::clearGeom();

    // Optionally read volScalarField
    autoPtr<volScalarField> readFieldPtr_;

    // 1. see if field in database
    // 2. see if field can be read
    const volScalarField* cellFldPtr = nullptr;
    if (fvm.foundObject<volScalarField>(isoField_))
    {
        if (debug)
        {
            InfoInFunction << "Lookup " << isoField_ << endl;
        }

        cellFldPtr = &fvm.lookupObject<volScalarField>(isoField_);
    }
    else
    {
        // Bit of a hack. Read field and store.

        if (debug)
        {
            InfoInFunction
                << "Reading " << isoField_
                << " from time " <<fvm.time().timeName()
                << endl;
        }

        readFieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    isoField_,
                    fvm.time().timeName(),
                    fvm,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fvm
            )
        );

        cellFldPtr = readFieldPtr_.operator->();
    }
    const volScalarField& cellFld = *cellFldPtr;

    tmp<pointScalarField> pointFld
    (
        volPointInterpolation::New(fvm).interpolate(cellFld)
    );

    if (average_)
    {
        //- From point field and interpolated cell.
        scalarField cellAvg(fvm.nCells(), scalar(0));
        labelField nPointCells(fvm.nCells(), 0);
        {
            for (label pointi = 0; pointi < fvm.nPoints(); pointi++)
            {
                const labelList& pCells = fvm.pointCells(pointi);

                forAll(pCells, i)
                {
                    label celli = pCells[i];

                    cellAvg[celli] += pointFld().primitiveField()[pointi];
                    nPointCells[celli]++;
                }
            }
        }
        forAll(cellAvg, celli)
        {
            cellAvg[celli] /= nPointCells[celli];
        }

        const isoSurfaceCell iso
        (
            fvm,
            cellAvg,
            pointFld().primitiveField(),
            isoVal_,
            regularise_
        );

        const_cast<sampledIsoSurfaceCell&>
        (
            *this
        ).triSurface::operator=(iso);
        meshCells_ = iso.meshCells();
    }
    else
    {
        //- Direct from cell field and point field. Gives bad continuity.
        const isoSurfaceCell iso
        (
            fvm,
            cellFld.primitiveField(),
            pointFld().primitiveField(),
            isoVal_,
            regularise_
        );

        const_cast<sampledIsoSurfaceCell&>
        (
            *this
        ).triSurface::operator=(iso);
        meshCells_ = iso.meshCells();
    }


    if (debug)
    {
        Pout<< "sampledIsoSurfaceCell::updateGeometry() : constructed iso:"
            << nl
            << "    regularise     : " << regularise_ << nl
            << "    average        : " << average_ << nl
            << "    isoField       : " << isoField_ << nl
            << "    isoValue       : " << isoVal_ << nl
            << "    points         : " << points().size() << nl
            << "    tris           : " << triSurface::size() << nl
            << "    cut cells      : " << meshCells_.size() << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledIsoSurfaceCell::sampledIsoSurfaceCell
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    isoField_(dict.lookup("isoField")),
    isoVal_(readScalar(dict.lookup("isoValue"))),
    regularise_(dict.lookupOrDefault("regularise", true)),
    average_(dict.lookupOrDefault("average", true)),
    zoneKey_(keyType::null),
    facesPtr_(nullptr),
    prevTimeIndex_(-1),
    meshCells_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledIsoSurfaceCell::~sampledIsoSurfaceCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledIsoSurfaceCell::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledIsoSurfaceCell::expire()
{
    facesPtr_.clear();

    // Clear derived data
    sampledSurface::clearGeom();

    // already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledIsoSurfaceCell::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurfaceCell::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurfaceCell::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurfaceCell::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurfaceCell::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurfaceCell::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurfaceCell::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurfaceCell::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurfaceCell::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurfaceCell::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurfaceCell::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledIsoSurfaceCell::print(Ostream& os) const
{
    os  << "sampledIsoSurfaceCell: " << name() << " :"
        << "  field:" << isoField_
        << "  value:" << isoVal_;
        //<< "  faces:" << faces().size()   // possibly no geom yet
        //<< "  points:" << points().size();
}


// ************************************************************************* //
