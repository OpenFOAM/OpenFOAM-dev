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

#include "sampledIsoSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{
    defineTypeNameAndDebug(isoSurface, 0);
    addToRunTimeSelectionTable(sampledSurface, isoSurface, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::cutPolyIsoSurface>
Foam::sampledSurfaces::isoSurface::calcIsoSurf() const
{
    // Lookup the field
    const volScalarField& vField =
        mesh().lookupObject<volScalarField>(isoFieldName_);

    // Interpolate the field to the points
    tmp<pointScalarField> pField
    (
        volPointInterpolation::New(vField.mesh()).interpolate(vField)
    );

    // Construct iso-surfaces for the given iso-values
    PtrList<cutPolyIsoSurface> isoSurfs(isoValues_.size());
    forAll(isoValues_, i)
    {
        isoSurfs.set
        (
            i,
            new cutPolyIsoSurface(mesh(), pField, isoValues_[i], zoneIDs())
        );
    }

    // Return or combine into a single surface
    if (isoValues_.size() == 1)
    {
        return autoPtr<cutPolyIsoSurface>(isoSurfs.set(0, nullptr));
    }
    else
    {
        return autoPtr<cutPolyIsoSurface>(new cutPolyIsoSurface(isoSurfs));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::isoSurface::isoSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledIsoSurfaceSurface(name, mesh, dict),
    isoFieldName_(dict.lookup("isoField")),
    isoValues_
    (
        dict.found("isoValues")
      ? scalarField(dict.lookup("isoValues"))
      : scalarField(1, dict.lookup<scalar>("isoValue"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::isoSurface::~isoSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::sampledSurfaces::isoSurface::fields() const
{
    return wordList(isoFieldName_);
}


bool Foam::sampledSurfaces::isoSurface::needsUpdate() const
{
    return timeIndex() != mesh().time().timeIndex();
}


void Foam::sampledSurfaces::isoSurface::print(Ostream& os) const
{
    os  << "isoSurface: " << name() << " :"
        << "  field:" << isoFieldName_;
    if (isoValues_.size() == 1)
    {
        os << "  value:" << isoValues_[0];
    }
    else
    {
        os << "  values:" << isoValues_;
    }
    os  << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
