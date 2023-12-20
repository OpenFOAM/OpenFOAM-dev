/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{
    defineTypeNameAndDebug(patchInternalField, 0);
    addToRunTimeSelectionTable(sampledSurface, patchInternalField, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::patchInternalField::patchInternalField
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    patch(name, mesh, dict),
    mappers_(patchIndices().size())
{
    // Negate the distance so that we sample cells inside the patch
    dictionary mappersDict(dict);
    if (dict.found("distance"))
    {
        mappersDict.set("distance", -mappersDict.lookup<scalar>("distance"));
    }

    forAll(patchIndices(), i)
    {
        mappers_.set
        (
            i,
            new mappedInternalPatchBase
            (
                mesh.boundaryMesh()[patchIndices()[i]],
                mappersDict
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::patchInternalField::~patchInternalField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::sampledSurfaces::patchInternalField::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledSurfaces::patchInternalField::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledSurfaces::patchInternalField::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledSurfaces::patchInternalField::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledSurfaces::patchInternalField::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledSurfaces::patchInternalField::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledSurfaces::patchInternalField::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledSurfaces::patchInternalField::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledSurfaces::patchInternalField::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledSurfaces::patchInternalField::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledSurfaces::patchInternalField::print(Ostream& os) const
{
    os  << "patchInternalField: " << name() << " :"
        << "  patches:" << patchNames()
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
