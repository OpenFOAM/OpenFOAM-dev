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

#include "sampledPatchInternalField.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "volFields.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledPatchInternalField, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledPatchInternalField,
        word,
        patchInternalField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledPatchInternalField::sampledPatchInternalField
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledPatch(name, mesh, dict),
    mappers_(patchIDs().size())
{
    mappedPatchBase::offsetMode mode = mappedPatchBase::NORMAL;
    if (dict.found("offsetMode"))
    {
        mode = mappedPatchBase::offsetModeNames_.read
        (
            dict.lookup("offsetMode")
        );
    }

    switch (mode)
    {
        case mappedPatchBase::NORMAL:
        {
            const scalar distance = readScalar(dict.lookup("distance"));
            forAll(patchIDs(), i)
            {
                mappers_.set
                (
                    i,
                    new mappedPatchBase
                    (
                        mesh.boundaryMesh()[patchIDs()[i]],
                        mesh.name(),                        // sampleRegion
                        mappedPatchBase::NEARESTCELL,       // sampleMode
                        word::null,                         // samplePatch
                        -distance                  // sample inside my domain
                    )
                );
            }
        }
        break;

        case mappedPatchBase::UNIFORM:
        {
            const point offset(dict.lookup("offset"));
            forAll(patchIDs(), i)
            {
                mappers_.set
                (
                    i,
                    new mappedPatchBase
                    (
                        mesh.boundaryMesh()[patchIDs()[i]],
                        mesh.name(),                        // sampleRegion
                        mappedPatchBase::NEARESTCELL,       // sampleMode
                        word::null,                         // samplePatch
                        offset                  // sample inside my domain
                    )
                );
            }
        }
        break;

        case mappedPatchBase::NONUNIFORM:
        {
            const pointField offsets(dict.lookup("offsets"));
            forAll(patchIDs(), i)
            {
                mappers_.set
                (
                    i,
                    new mappedPatchBase
                    (
                        mesh.boundaryMesh()[patchIDs()[i]],
                        mesh.name(),                        // sampleRegion
                        mappedPatchBase::NEARESTCELL,       // sampleMode
                        word::null,                         // samplePatch
                        offsets                  // sample inside my domain
                    )
                );
            }
        }
        break;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledPatchInternalField::~sampledPatchInternalField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::sampledPatchInternalField::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::sampledPatchInternalField::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledPatchInternalField::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPatchInternalField::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::sampledPatchInternalField::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::sampledPatchInternalField::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledPatchInternalField::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledPatchInternalField::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPatchInternalField::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledPatchInternalField::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledPatchInternalField::print(Ostream& os) const
{
    os  << "sampledPatchInternalField: " << name() << " :"
        << "  patches:" << patchNames()
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
