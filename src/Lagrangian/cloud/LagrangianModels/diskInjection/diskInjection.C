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

#include "LagrangianMesh.H"
#include "diskInjection.H"
#include "addToRunTimeSelectionTable.H"
#include "LagrangianFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(diskInjection, 0);
    addToRunTimeSelectionTable(LagrangianModel, diskInjection, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Lagrangian::diskInjection::readCoeffs(const dictionary& modelDict)
{
    centre_.reset
    (
        Function1<vector>::New
        (
            "centre",
            mesh().time().userUnits(),
            dimLength,
            modelDict
        ).ptr()
    );

    axis_.reset
    (
        Function1<vector>::New
        (
            "axis",
            mesh().time().userUnits(),
            dimless,
            modelDict
        ).ptr()
    );

    const bool haveDiemeter = modelDict.found("diameter");
    const bool haveInnerDiemeter = modelDict.found("innerDiameter");
    const bool haveOuterDiemeter = modelDict.found("outerDiameter");

    if (haveDiemeter == (haveInnerDiemeter || haveOuterDiemeter))
    {
        FatalIOErrorInFunction(modelDict)
            << "keywords diameter and innerDiameter/outerDiameter are both "
            << (haveDiemeter ? "" : "un") << "defined in "
            << "dictionary " << modelDict.name()
            << exit(FatalIOError);
    }

    if (haveInnerDiemeter != haveOuterDiemeter)
    {
        FatalIOErrorInFunction(modelDict)
            << "keywords innerDiameter and outerDiameter are not both defined "
            << "in dictionary " << modelDict.name()
            << exit(FatalIOError);
    }

    if (haveDiemeter)
    {
        innerDiameter_ = 0;
        outerDiameter_ = modelDict.lookup<scalar>("diameter", dimLength);
    }
    else
    {
        innerDiameter_ = modelDict.lookup<scalar>("innerDiameter", dimLength);
        outerDiameter_ = modelDict.lookup<scalar>("outerDiameter", dimLength);
    }

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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::diskInjection::diskInjection
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianInjection(name, mesh),
    centre_(nullptr),
    axis_(nullptr),
    innerDiameter_(NaN),
    outerDiameter_(NaN),
    numberRate_(nullptr),
    numberDeferred_(stateDict.lookupOrDefault<scalar>("numberDeferred", 0)),
    rndGen_("rndGen", stateDict, name, true),
    timeIndex_(-1)
{
    readCoeffs(modelDict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::scalarField& Foam::Lagrangian::diskInjection::rFrac() const
{
    if (!rFracPtr_.valid())
    {
        FatalErrorInFunction
            << "Radius fraction requested outside of the injection"
            << exit(FatalError);
    }

    return rFracPtr_();
}


const Foam::vectorField& Foam::Lagrangian::diskInjection::axis() const
{
    if (!axisPtr_.valid())
    {
        FatalErrorInFunction
            << "Axis requested outside of the injection"
            << exit(FatalError);
    }

    return axisPtr_();
}


const Foam::vectorField& Foam::Lagrangian::diskInjection::radial() const
{
    if (!radialPtr_.valid())
    {
        FatalErrorInFunction
            << "Axis requested outside of the injection"
            << exit(FatalError);
    }

    return radialPtr_();
}


Foam::LagrangianSubMesh Foam::Lagrangian::diskInjection::modify
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

    // Evaluate the variable centre and axis, and create radial vectors to
    // complete the local coordinate system
    const pointField centre(centre_->value(t0 + fraction*(t1 - t0)));
    axisPtr_.set(normalised(axis_->value(t0 + fraction*(t1 - t0))).ptr());
    const vectorField radial1(normalised(perpendicular(axis())));
    const vectorField radial2(axis() ^ radial1);

    // Create random radii within and angles around the disk. Store temporarily
    // so that source conditions can use it.
    rFracPtr_.set(rndGen_.scalar01(numberInt).ptr());
    const scalarField r
    (
        sqrt
        (
            (1 - rFrac())*sqr(innerDiameter_/2)
          + rFrac()*sqr(outerDiameter_/2)
        )
    );
    const scalarField phi
    (
        constant::mathematical::twoPi*rndGen_.scalar01(numberInt)
    );

    // Evaluate the radial vector
    radialPtr_.set((cos(phi)*radial1 + sin(phi)*radial2).ptr());

    // Evaluate the positions
    const pointField positions(centre + r*radial());

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
    LagrangianSubMesh injectionMesh =
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

    // Clean up
    rFracPtr_.clear();
    axisPtr_.clear();
    radialPtr_.clear();

    return injectionMesh;
}


bool Foam::Lagrangian::diskInjection::read(const dictionary& modelDict)
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


void Foam::Lagrangian::diskInjection::writeState(Ostream& os) const
{
    LagrangianInjection::writeState(os);

    writeEntry(os, "numberDeferred", numberDeferred_);
    writeEntry(os, "rndGen", rndGen_);
}


// ************************************************************************* //
