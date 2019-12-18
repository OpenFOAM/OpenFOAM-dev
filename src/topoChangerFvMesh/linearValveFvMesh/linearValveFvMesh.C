/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "linearValveFvMesh.H"
#include "Time.H"
#include "slidingInterface.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearValveFvMesh, 0);

    addToRunTimeSelectionTable(topoChangerFvMesh, linearValveFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::linearValveFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if
    (
        pointZones().size()
     || faceZones().size()
     || cellZones().size()
     || topoChanger_.size()
    )
    {
        InfoInFunction
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    // Add zones
    List<pointZone*> pz(1);

    // Add an empty zone for cut points

    pz[0] = new pointZone
    (
        "cutPointZone",
        labelList(0),
        0,
        pointZones()
    );


    // Do face zones for slider

    List<faceZone*> fz(3);

    // Inner slider
    const word innerSliderName(motionDict_.subDict("slider").lookup("inside"));
    const polyPatch& innerSlider = boundaryMesh()[innerSliderName];

    labelList isf(innerSlider.size());

    forAll(isf, i)
    {
        isf[i] = innerSlider.start() + i;
    }

    fz[0] = new faceZone
    (
        "insideSliderZone",
        isf,
        boolList(innerSlider.size(), false),
        0,
        faceZones()
    );

    // Outer slider
    const word outerSliderName(motionDict_.subDict("slider").lookup("outside"));
    const polyPatch& outerSlider = boundaryMesh()[outerSliderName];

    labelList osf(outerSlider.size());

    forAll(osf, i)
    {
        osf[i] = outerSlider.start() + i;
    }

    fz[1] = new faceZone
    (
        "outsideSliderZone",
        osf,
        boolList(outerSlider.size(), false),
        1,
        faceZones()
    );

    // Add empty zone for cut faces
    fz[2] = new faceZone
    (
        "cutFaceZone",
        labelList(0),
        boolList(0, false),
        2,
        faceZones()
    );

    List<cellZone*> cz(0);

    Info<< "Adding point, face and cell zones" << endl;
    addZones(pz, fz, cz);

    // Add a topology modifier
    Info<< "Adding topology modifiers" << endl;
    topoChanger_.setSize(1);
    topoChanger_.set
    (
        0,
        new slidingInterface
        (
            "mixerSlider",
            0,
            topoChanger_,
            outerSliderName + "Zone",
            innerSliderName + "Zone",
            "cutPointZone",
            "cutFaceZone",
            outerSliderName,
            innerSliderName,
            slidingInterface::INTEGRAL,
            true                          // Attach-detach action
        )
    );
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;

    // Write mesh
    write();
}


void Foam::linearValveFvMesh::makeSlidersDead()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable layering
    forAll(topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].disable();
        }
        else
        {
            FatalErrorInFunction
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::linearValveFvMesh::makeSlidersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable sliding interface
    forAll(topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else
        {
            FatalErrorInFunction
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


bool Foam::linearValveFvMesh::attached() const
{
    const polyTopoChanger& topoChanges = topoChanger_;

    bool result = false;

    forAll(topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            result =
                result
             || refCast<const slidingInterface>(topoChanges[modI]).attached();
        }
    }

    // Check thal all sliders are in sync (debug only)
    forAll(topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            if
            (
                result
             != refCast<const slidingInterface>(topoChanges[modI]).attached()
            )
            {
                FatalErrorInFunction
                    << "Slider " << modI
                    << " named " << topoChanges[modI].name()
                    << " out of sync: Should be" << result
                    << abort(FatalError);
            }
        }
    }

    if (result)
    {
        Info<< "linearValveFvMesh: attached!" << endl;
    }
    else
    {
        Info<< "linearValveFvMesh: detached!" << endl;
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearValveFvMesh::linearValveFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io),
    motionDict_(dict().optionalSubDict(typeName + "Coeffs")),
    msPtr_(motionSolver::New(*this), dict())
{
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearValveFvMesh::~linearValveFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::linearValveFvMesh::update()
{
    // Detaching the interface
    if (attached())
    {
        Info<< "Decoupling sliding interfaces" << endl;
        makeSlidersLive();

        // Changing topology by hand
        resetMorph();
        setMorphTimeIndex(3*time().timeIndex());
        updateMesh();

        msPtr_->updateMesh();
    }
    else
    {
        Info<< "Sliding interfaces decoupled" << endl;
    }

    // Perform mesh motion
    makeSlidersDead();

    // Changing topology by hand
    resetMorph();
    setMorphTimeIndex(3*time().timeIndex() + 1);
    updateMesh();

    msPtr_->updateMesh();

    if (topoChangeMap.valid())
    {
        if (topoChangeMap().hasMotionPoints())
        {
            Info<< "Topology change; executing pre-motion" << endl;
            movePoints(topoChangeMap().preMotionPoints());
        }
    }

    // Solve for motion
    msPtr_->solve();

    movePoints(msPtr_->curPoints());

    // Attach the interface
    Info<< "Coupling sliding interfaces" << endl;
    makeSlidersLive();
    resetMorph();
    setMorphTimeIndex(3*time().timeIndex() + 2);
    updateMesh();

    Info<< "Moving points post slider attach" << endl;

    msPtr_->updateMesh();

    Info<< "Sliding interfaces coupled: " << attached() << endl;
}


// ************************************************************************* //
