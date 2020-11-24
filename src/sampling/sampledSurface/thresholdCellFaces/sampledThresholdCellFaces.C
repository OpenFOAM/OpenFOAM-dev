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

#include "sampledThresholdCellFaces.H"
#include "thresholdCellFaces.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{
    defineTypeNameAndDebug(thresholdCellFaces, 0);
    addToRunTimeSelectionTable(sampledSurface, thresholdCellFaces, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::sampledSurfaces::thresholdCellFaces::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // no update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    prevTimeIndex_ = fvm.time().timeIndex();

    // Optionally read volScalarField
    autoPtr<volScalarField> readFieldPtr_;

    // 1. see if field in database
    // 2. see if field can be read
    const volScalarField* cellFldPtr = nullptr;
    if (fvm.foundObject<volScalarField>(fieldName_))
    {
        if (debug)
        {
            InfoInFunction<< "Lookup " << fieldName_ << endl;
        }

        cellFldPtr = &fvm.lookupObject<volScalarField>(fieldName_);
    }
    else
    {
        // Bit of a hack. Read field and store.

        if (debug)
        {
            InfoInFunction
                << "Reading " << fieldName_
                << " from time " << fvm.time().timeName()
                << endl;
        }

        readFieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName_,
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


    Foam::thresholdCellFaces surf
    (
        fvm,
        cellFld.primitiveField(),
        lowerThreshold_,
        upperThreshold_,
        triangulate_
    );

    const_cast<thresholdCellFaces&>
    (
        *this
    ).MeshedSurface<face>::transfer(surf);
    meshCells_.transfer(surf.meshCells());

    // clear derived data
    sampledSurface::clearGeom();

    if (debug)
    {
        Pout<< "thresholdCellFaces::updateGeometry() : constructed"
            << nl
            << "    field         : " << fieldName_ << nl
            << "    lowerLimit    : " << lowerThreshold_ << nl
            << "    upperLimit    : " << upperThreshold_ << nl
            << "    point         : " << points().size() << nl
            << "    faces         : " << faces().size() << nl
            << "    cut cells     : " << meshCells_.size() << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::thresholdCellFaces::thresholdCellFaces
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    fieldName_(dict.lookup("field")),
    lowerThreshold_(dict.lookupOrDefault<scalar>("lowerLimit", -vGreat)),
    upperThreshold_(dict.lookupOrDefault<scalar>("upperLimit", vGreat)),
    triangulate_(dict.lookupOrDefault("triangulate", false)),
    prevTimeIndex_(-1),
    meshCells_(0)
{
    if (!dict.found("lowerLimit") && !dict.found("upperLimit"))
    {
        FatalErrorInFunction
            << "require at least one of 'lowerLimit' or 'upperLimit'" << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::thresholdCellFaces::~thresholdCellFaces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledSurfaces::thresholdCellFaces::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledSurfaces::thresholdCellFaces::expire()
{
    // already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledSurfaces::thresholdCellFaces::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField>
Foam::sampledSurfaces::thresholdCellFaces::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledSurfaces::thresholdCellFaces::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledSurfaces::thresholdCellFaces::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledSurfaces::thresholdCellFaces::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledSurfaces::thresholdCellFaces::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledSurfaces::thresholdCellFaces::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledSurfaces::thresholdCellFaces::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledSurfaces::thresholdCellFaces::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledSurfaces::thresholdCellFaces::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledSurfaces::thresholdCellFaces::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledSurfaces::thresholdCellFaces::print(Ostream& os) const
{
    os  << "thresholdCellFaces: " << name() << " :"
        << "  field:" << fieldName_
        << "  lowerLimit:" << lowerThreshold_
        << "  upperLimit:" << upperThreshold_;
        //<< "  faces:" << faces().size()   // possibly no geom yet
        //<< "  points:" << points().size();
}


// ************************************************************************* //
