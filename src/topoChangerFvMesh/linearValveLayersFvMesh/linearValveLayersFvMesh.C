/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "linearValveLayersFvMesh.H"
#include "Time.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "pointField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearValveLayersFvMesh, 0);
    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        linearValveLayersFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::linearValveLayersFvMesh::addZonesAndModifiers()
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
        Info<< "void linearValveLayersFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    // Add zones
    List<pointZone*> pz(1);
    List<faceZone*> fz(4);
    List<cellZone*> cz(0);


    // Add an empty zone for cut points

    pz[0] = new pointZone
    (
        "cutPointZone",
        labelList(0),
        0,
        pointZones()
    );


    // Do face zones for slider

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

    // Add face zone for layer addition
    const word layerPatchName
    (
        motionDict_.subDict("layer").lookup("patch")
    );

    const polyPatch& layerPatch = boundaryMesh()[layerPatchName];

    labelList lpf(layerPatch.size());

    forAll(lpf, i)
    {
        lpf[i] = layerPatch.start() + i;
    }

    fz[3] = new faceZone
    (
        "valveLayerZone",
        lpf,
        boolList(layerPatch.size(), true),
        0,
        faceZones()
    );


    Info<< "Adding point and face zones" << endl;
    addZones(pz, fz, cz);

    // Add a topology modifier

    List<polyMeshModifier*> tm(2);

    tm[0] = new slidingInterface
    (
        "valveSlider",
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
    );

    tm[1] =
        new layerAdditionRemoval
        (
            "valveLayer",
            1,
            topoChanger_,
            "valveLayerZone",
            readScalar
            (
                motionDict_.subDict("layer").lookup("minThickness")
            ),
            readScalar
            (
                motionDict_.subDict("layer").lookup("maxThickness")
            )
        );


    Info<< "Adding topology modifiers" << endl;
    addTopologyModifiers(tm);

    // Write mesh
    write();
}


void Foam::linearValveLayersFvMesh::makeLayersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable layering
    forAll(topoChanges, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].disable();
        }
        else
        {
            FatalErrorIn("void linearValveLayersFvMesh::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::linearValveLayersFvMesh::makeSlidersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable sliding interface
    forAll(topoChanges, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanges[modI]))
        {
            topoChanges[modI].disable();
        }
        else if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else
        {
            FatalErrorIn("void linearValveLayersFvMesh::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


bool Foam::linearValveLayersFvMesh::attached() const
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
                FatalErrorIn("bool linearValveLayersFvMesh::attached() const")
                    << "Slider " << modI << " named "
                    << topoChanges[modI].name()
                    << " out of sync: Should be" << result
                    << abort(FatalError);
            }
        }
    }

    return result;
}


Foam::tmp<Foam::pointField> Foam::linearValveLayersFvMesh::newPoints() const
{
    tmp<pointField> tnewPoints
    (
        new pointField(points())
    );

    pointField& np = tnewPoints();

    const word layerPatchName
    (
        motionDict_.subDict("layer").lookup("patch")
    );

    const polyPatch& layerPatch = boundaryMesh()[layerPatchName];

    const labelList& patchPoints = layerPatch.meshPoints();

    const vector vel
    (
        motionDict_.lookup("pistonVelocity")
    );

    forAll(patchPoints, ppI)
    {
        np[patchPoints[ppI]] += vel*time().deltaTValue();
    }

    return tnewPoints;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::linearValveLayersFvMesh::linearValveLayersFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io),
    motionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    )
{
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearValveLayersFvMesh::~linearValveLayersFvMesh()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::linearValveLayersFvMesh::update()
{
    // Detaching the interface
    if (attached())
    {
        Info<< "Decoupling sliding interfaces" << endl;
        makeSlidersLive();

        // Changing topology
        resetMorph();
        setMorphTimeIndex(3*time().timeIndex());
        updateMesh();
    }
    else
    {
        Info<< "Sliding interfaces decoupled" << endl;
    }

    // Perform layer action and mesh motion
    makeLayersLive();

    // Changing topology
    resetMorph();
    setMorphTimeIndex(3*time().timeIndex() + 1);
    updateMesh();

    if (topoChangeMap.valid())
    {
        if (topoChangeMap().hasMotionPoints())
        {
            Info<< "Topology change; executing pre-motion" << endl;
            movePoints(topoChangeMap().preMotionPoints());
        }
    }

    // Move points
    movePoints(newPoints());

    // Attach the interface
    Info<< "Coupling sliding interfaces" << endl;
    makeSlidersLive();

    // Changing topology
    resetMorph();
    setMorphTimeIndex(3*time().timeIndex() + 2);
    updateMesh();

    Info<< "Moving points post slider attach" << endl;
//     const pointField p = allPoints();
//     movePoints(p);

    Info<< "Sliding interfaces coupled: " << attached() << endl;
}


// ************************************************************************* //

