/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "manualInjection.H"
#include "LagrangianFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(manualInjection, 0);
    addToRunTimeSelectionTable(LagrangianModel, manualInjection, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Lagrangian::manualInjection::readCoeffs(const dictionary& modelDict)
{
    fileName_ = modelDict.lookup<fileName>("file");

    units_.readIfPresent("units", modelDict);

    time_ =
        modelDict.lookupOrDefault<scalar>
        (
            "time",
            mesh().time().userUnits(),
            mesh().time().beginTime().value()
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::manualInjection::manualInjection
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianInjection(name, mesh),
    fileName_(fileName::null),
    units_(dimLength),
    time_(NaN)
{
    readCoeffs(modelDict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::LagrangianSubMesh Foam::Lagrangian::manualInjection::modify
(
    LagrangianMesh& mesh,
    const LagrangianSubMesh&
) const
{
    const scalar t1 = mesh.time().value();
    const scalar t0 = t1 - mesh.time().deltaT().value();

    if (!(t0 <= time_ && time_ < t1)) return mesh.subNone();

    // Read the positions file
    GlobalIOField<point> positions
    (
        IOobject
        (
            fileName_,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    units_.makeStandard(positions);

    const label number = positions.size();

    scalarField fraction(number, (time_ - t0)/(t1 - t0));

    // Locate within the mesh
    barycentricField coordinates(number);
    labelField celli(number, -1), facei(number), faceTrii(number);
    const List<LagrangianMesh::location> locations =
        mesh.locate
        (
            positions,
            coordinates,
            celli,
            facei,
            faceTrii,
            fraction
        );

    // Check for any failures
    checkLocation(locations, positions);

    // Remove particles not on this process
    filter(coordinates, celli, facei, faceTrii, fraction);

    // Inject particles
    return
        mesh.inject
        (
            *this,
            coordinates,
            celli,
            facei,
            faceTrii,
            LagrangianMesh::fractionName,
            fraction
        );
}


bool Foam::Lagrangian::manualInjection::read(const dictionary& modelDict)
{
    if (LagrangianInjection::read(modelDict))
    {
        readCoeffs(modelDict);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
