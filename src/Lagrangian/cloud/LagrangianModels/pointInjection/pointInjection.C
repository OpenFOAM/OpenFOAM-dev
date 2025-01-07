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

#include "pointInjection.H"
#include "addToRunTimeSelectionTable.H"
#include "LagrangianFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(pointInjection, 0);
    addToRunTimeSelectionTable(LagrangianModel, pointInjection, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Lagrangian::pointInjection::readCoeffs(const dictionary& modelDict)
{
    point_.reset
    (
        Function1<vector>::New
        (
            "point",
            mesh().time().userUnits(),
            dimLength,
            modelDict
        ).ptr()
    );

    numberRate_.reset
    (
        Function1<scalar>::New
        (
            "numberRate",
            mesh().time().userUnits(),
            dimRate,
            modelDict
        ).ptr()
    );

    numberDeferred_ = 0;

    injectionLocation_ = injectionLocation::unset;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::pointInjection::pointInjection
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianInjection(name, mesh),
    point_(nullptr),
    numberRate_(nullptr),
    numberDeferred_(stateDict.lookupOrDefault<scalar>("numberDeferred", 0)),
    rndGen_("rndGen", stateDict, name, true),
    timeIndex_(-1),
    injectionLocation_(injectionLocation::unset),
    coordinates_(barycentric::uniform(NaN)),
    celli_(-1),
    facei_(-1),
    faceTrii_(-1)
{
    readCoeffs(modelDict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::Lagrangian::pointInjection::correct()
{
    // If anything is moving then all coordinates have to be constructed on the
    // fly, so don't store anything
    if (mesh().mesh().moving() || !point_->constant())
    {
        injectionLocation_ =
            injectionLocation::multiplePoints;
    }

    // If nothing is moving, and coordinates are not already stored, then
    // construct and store the coordinates for repeated usage
    if
    (
        injectionLocation_ == injectionLocation::unset
     || injectionLocation_ == injectionLocation::multiplePoints
    )
    {
        const scalar t1 = mesh().time().value();
        const scalar t0 = t1 - mesh().time().deltaT().value();

        // Evaluate the position
        const point p = point_->value(t0);

        // Locate within the mesh
        const LagrangianMesh::location l =
            mesh().locate(p, coordinates_, celli_, facei_, faceTrii_, 0);

        // Check for failure
        checkLocation(l, p);

        // Set the location flag
        injectionLocation_ =
            celli_ >= 0
          ? injectionLocation::fixedPointOnThisProcessor
          : injectionLocation::fixedPointOnAnotherProcessor;
    }
}


Foam::LagrangianSubMesh Foam::Lagrangian::pointInjection::modify
(
    LagrangianMesh& mesh,
    const LagrangianSubMesh&
) const
{
    const scalar t1 = mesh.time().value();
    const scalar t0 = t1 - mesh.time().deltaT().value();

    // Restart the generator if necessary and set the time index up to date
    rndGen_.start(timeIndex_ == db().time().timeIndex());
    timeIndex_ = db().time().timeIndex();

    // Calculate the number of particles to inject. Round down to get an
    // integer number. Store the excess to apply at a later time.
    const scalar number = numberRate_->integral(t0, t1) + numberDeferred_;
    const label numberInt = floor(number);
    numberDeferred_ = number - numberInt;

    // Inject at random times throughout the time-step
    scalarField fraction(rndGen_.scalar01(numberInt));

    // Create particles in the Lagrangian mesh
    switch (injectionLocation_)
    {
        case injectionLocation::unset:
            FatalErrorInFunction
                << "Injection location has not been set"
                << exit(FatalError);
            return LagrangianSubMesh(mesh, LagrangianGroup::none, -1, -1);

        case injectionLocation::fixedPointOnThisProcessor:
        {
            // Use stored coordinates to inject the particles
            return
                mesh.inject
                (
                    *this,
                    barycentricField(numberInt, coordinates_),
                    labelField(numberInt, celli_),
                    labelField(numberInt, facei_),
                    labelField(numberInt, faceTrii_),
                    LagrangianMesh::fractionName,
                    fraction
                );
        }

        case injectionLocation::fixedPointOnAnotherProcessor:
        {
            // Inject nothing. We need this because source conditions evaluated
            // on the process containing the point may communicate; e.g., to
            // sum the flow rate across the processes.
            return
                mesh.inject
                (
                    *this,
                    barycentricField(),
                    labelField(),
                    labelField(),
                    labelField(),
                    LagrangianMesh::fractionName,
                    scalarField()
                );
        }

        case injectionLocation::multiplePoints:
        {
            // Evaluate the variable positions
            const pointField positions(point_->value(t0 + fraction*(t1 - t0)));

            // Locate within the moving mesh
            barycentricField coordinates(numberInt);
            labelField celli(numberInt), facei(numberInt), faceTrii(numberInt);
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
    }

    return LagrangianSubMesh(mesh, LagrangianGroup::none, -1, -1);
}


void Foam::Lagrangian::pointInjection::topoChange(const polyTopoChangeMap&)
{
    injectionLocation_ = injectionLocation::unset;
}


void Foam::Lagrangian::pointInjection::mapMesh(const polyMeshMap&)
{
    injectionLocation_ = injectionLocation::unset;
}


void Foam::Lagrangian::pointInjection::distribute(const polyDistributionMap&)
{
    injectionLocation_ = injectionLocation::unset;
}


bool Foam::Lagrangian::pointInjection::read(const dictionary& modelDict)
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


void Foam::Lagrangian::pointInjection::writeState(Ostream& os) const
{
    LagrangianInjection::writeState(os);

    writeEntry(os, "numberDeferred", numberDeferred_);
    writeEntry(os, "rndGen", rndGen_);
}


// ************************************************************************* //
