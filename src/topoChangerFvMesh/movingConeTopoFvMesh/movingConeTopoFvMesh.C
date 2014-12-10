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

#include "movingConeTopoFvMesh.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "layerAdditionRemoval.H"
#include "addToRunTimeSelectionTable.H"
#include "meshTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingConeTopoFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        movingConeTopoFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::movingConeTopoFvMesh::vertexMarkup
(
    const pointField& p,
    const scalar curLeft,
    const scalar curRight
) const
{
    Info<< "Updating vertex markup.  curLeft: "
        << curLeft << " curRight: " << curRight << endl;

    tmp<scalarField> tvertexMarkup(new scalarField(p.size()));
    scalarField& vertexMarkup = tvertexMarkup();

    forAll(p, pI)
    {
        if (p[pI].x() < curLeft - SMALL)
        {
            vertexMarkup[pI] = -1;
        }
        else if (p[pI].x() > curRight + SMALL)
        {
            vertexMarkup[pI] = 1;
        }
        else
        {
            vertexMarkup[pI] = 0;
        }
    }

    return tvertexMarkup;
}


void Foam::movingConeTopoFvMesh::addZonesAndModifiers()
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
        Info<< "void movingConeTopoFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    const vectorField& fc = faceCentres();
    const vectorField& fa = faceAreas();

    labelList zone1(fc.size());
    boolList flipZone1(fc.size(), false);
    label nZoneFaces1 = 0;

    labelList zone2(fc.size());
    boolList flipZone2(fc.size(), false);
    label nZoneFaces2 = 0;

    forAll(fc, faceI)
    {
        if
        (
            fc[faceI].x() > -0.003501
         && fc[faceI].x() < -0.003499
        )
        {
            if ((fa[faceI] & vector(1, 0, 0)) < 0)
            {
                flipZone1[nZoneFaces1] = true;
            }

            zone1[nZoneFaces1] = faceI;
            Info<< "face " << faceI << " for zone 1.  Flip: "
                << flipZone1[nZoneFaces1] << endl;
            nZoneFaces1++;
        }
        else if
        (
            fc[faceI].x() > -0.00701
         && fc[faceI].x() < -0.00699
        )
        {
            zone2[nZoneFaces2] = faceI;

            if ((fa[faceI] & vector(1, 0, 0)) > 0)
            {
                flipZone2[nZoneFaces2] = true;
            }

            Info<< "face " << faceI << " for zone 2.  Flip: "
                << flipZone2[nZoneFaces2] << endl;
            nZoneFaces2++;
        }
    }

    zone1.setSize(nZoneFaces1);
    flipZone1.setSize(nZoneFaces1);

    zone2.setSize(nZoneFaces2);
    flipZone2.setSize(nZoneFaces2);

    Info<< "zone: " << zone1 << endl;
    Info<< "zone: " << zone2 << endl;

    List<pointZone*> pz(0);
    List<faceZone*> fz(2);
    List<cellZone*> cz(0);

    label nFz = 0;

    fz[nFz] =
        new faceZone
        (
            "rightExtrusionFaces",
            zone1,
            flipZone1,
            nFz,
            faceZones()
        );
    nFz++;

    fz[nFz] =
        new faceZone
        (
            "leftExtrusionFaces",
            zone2,
            flipZone2,
            nFz,
            faceZones()
        );
    nFz++;

    fz.setSize(nFz);

    Info<< "Adding mesh zones." << endl;
    addZones(pz, fz, cz);


    // Add layer addition/removal interfaces

    List<polyMeshModifier*> tm(2);
    label nMods = 0;

    tm[nMods] =
        new layerAdditionRemoval
        (
            "right",
            nMods,
            topoChanger_,
            "rightExtrusionFaces",
            readScalar
            (
                motionDict_.subDict("right").lookup("minThickness")
            ),
            readScalar
            (
                motionDict_.subDict("right").lookup("maxThickness")
            )
        );
    nMods++;

    tm[nMods] = new layerAdditionRemoval
    (
        "left",
        nMods,
        topoChanger_,
        "leftExtrusionFaces",
        readScalar
        (
            motionDict_.subDict("left").lookup("minThickness")
        ),
        readScalar
        (
            motionDict_.subDict("left").lookup("maxThickness")
        )
    );
    nMods++;
    tm.setSize(nMods);

    Info<< "Adding " << nMods << " mesh modifiers" << endl;
    topoChanger_.addTopologyModifiers(tm);

    write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::movingConeTopoFvMesh::movingConeTopoFvMesh(const IOobject& io)
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
    ),
    motionVelAmplitude_(motionDict_.lookup("motionVelAmplitude")),
    motionVelPeriod_(readScalar(motionDict_.lookup("motionVelPeriod"))),
    curMotionVel_
    (
        motionVelAmplitude_*
        Foam::sin(time().value()*M_PI/motionVelPeriod_)
    ),
    leftEdge_(readScalar(motionDict_.lookup("leftEdge"))),
    curLeft_(readScalar(motionDict_.lookup("leftObstacleEdge"))),
    curRight_(readScalar(motionDict_.lookup("rightObstacleEdge")))
{
    Pout<< "Initial time:" << time().value()
        << " Initial curMotionVel_:" << curMotionVel_
        << endl;

    addZonesAndModifiers();

    curLeft_ = average
    (
        faceZones()
        [
            faceZones().findZoneID("leftExtrusionFaces")
        ]().localPoints()
    ).x() - SMALL;

    curRight_ = average
    (
        faceZones()
        [
            faceZones().findZoneID("rightExtrusionFaces")
        ]().localPoints()
    ).x() + SMALL;

    motionMask_ = vertexMarkup
    (
        points(),
        curLeft_,
        curRight_
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingConeTopoFvMesh::~movingConeTopoFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::movingConeTopoFvMesh::update()
{
    // Do mesh changes (use inflation - put new points in topoChangeMap)
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh(true);

    // Calculate the new point positions depending on whether the
    // topological change has happened or not
    pointField newPoints;

    vector curMotionVel_ =
        motionVelAmplitude_*
        Foam::sin(time().value()*M_PI/motionVelPeriod_);

    Pout<< "time:" << time().value() << " curMotionVel_:" << curMotionVel_
        << " curLeft:" << curLeft_ << " curRight:" << curRight_
        << endl;

    if (topoChangeMap.valid())
    {
        Info<< "Topology change. Calculating motion points" << endl;

        if (topoChangeMap().hasMotionPoints())
        {
            Info<< "Topology change. Has premotion points" << endl;
            //Info<< "preMotionPoints:" << topoChangeMap().preMotionPoints()
            //    << endl;

            //mkDir(time().timePath());
            //{
            //    OFstream str(time().timePath()/"meshPoints.obj");
            //    Pout<< "Writing mesh with meshPoints to " << str.name()
            //        << endl;
            //
            //    const pointField& currentPoints = points();
            //    label vertI = 0;
            //    forAll(currentPoints, pointI)
            //    {
            //        meshTools::writeOBJ(str, currentPoints[pointI]);
            //        vertI++;
            //    }
            //    forAll(edges(), edgeI)
            //    {
            //        const edge& e = edges()[edgeI];
            //        str << "l " << e[0]+1 << ' ' << e[1]+1 << nl;
            //    }
            //}
            //{
            //    OFstream str(time().timePath()/"preMotionPoints.obj");
            //    Pout<< "Writing mesh with preMotionPoints to " << str.name()
            //        << endl;
            //
            //    const pointField& newPoints =
            //        topoChangeMap().preMotionPoints();
            //    label vertI = 0;
            //    forAll(newPoints, pointI)
            //    {
            //        meshTools::writeOBJ(str, newPoints[pointI]);
            //        vertI++;
            //    }
            //    forAll(edges(), edgeI)
            //    {
            //        const edge& e = edges()[edgeI];
            //        str << "l " << e[0]+1 << ' ' << e[1]+1 << nl;
            //    }
            //}


            motionMask_ =
                vertexMarkup
                (
                    topoChangeMap().preMotionPoints(),
                    curLeft_,
                    curRight_
                );

            // Move points inside the motionMask
            newPoints =
                topoChangeMap().preMotionPoints()
              + (
                    pos(0.5 - mag(motionMask_)) // cells above the body
                )*curMotionVel_*time().deltaT().value();
        }
        else
        {
            Info<< "Topology change. Already set mesh points" << endl;

            motionMask_ =
                vertexMarkup
                (
                    points(),
                    curLeft_,
                    curRight_
                );

            // Move points inside the motionMask
            newPoints =
                points()
              + (
                    pos(0.5 - mag(motionMask_)) // cells above the body
                )*curMotionVel_*time().deltaT().value();
        }
    }
    else
    {
        Info<< "No topology change" << endl;
        // Set the mesh motion
        newPoints =
            points()
          + (
                pos(0.5 - mag(motionMask_)) // cells above the body
           )*curMotionVel_*time().deltaT().value();
    }

    // The mesh now contains the cells with zero volume
    Info << "Executing mesh motion" << endl;
    movePoints(newPoints);
    //  The mesh now has got non-zero volume cells

    curLeft_ = average
    (
        faceZones()
        [
            faceZones().findZoneID("leftExtrusionFaces")
        ]().localPoints()
    ).x() - SMALL;

    curRight_ = average
    (
        faceZones()
        [
            faceZones().findZoneID("rightExtrusionFaces")
        ]().localPoints()
    ).x() + SMALL;


    return true;
}


// ************************************************************************* //
