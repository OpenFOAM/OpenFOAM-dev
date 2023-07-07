/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "snappyHexMeshConfiguration.H"
#include "Tuple3.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::snappyHexMeshConfiguration::writeSnappySwitches()
{
    dictionary dict("switches");

    dict.add("castellatedMesh", "on", true);
    dict.add("snap", "on", true);
    dict.add
    (
        "addLayers",
        layers_ == 0 ? "off" : "on",
        true
    );

    dict.write(os_, false);
    os_ << endl;
}


void Foam::snappyHexMeshConfiguration::writeGeometrySurface(const label surfID)
{
    beginDict(os_, surfaces_[surfID].name());

    os_ << indent << "type triSurfaceMesh;" << nl
        << indent << "file " << surfaces_[surfID].file() << ";" << endl;

    const wordList& inletRegions = surfaces_[surfID].inletRegions();
    const wordList& outletRegions = surfaces_[surfID].outletRegions();

    if (!inletRegions.empty() || !outletRegions.empty())
    {
        beginDict(os_, "regions");

        forAll(inletRegions, i)
        {
            const word& region = inletRegions[i];
            const word patch =
            (
                region == "<inletRegion>" || inletRegions.size() == 1
              ? "inlet"
              : region
            );

            os_ << indent << region << " { name " << patch << "; }" << endl;
        }

        forAll(outletRegions, i)
        {
            const word& region = outletRegions[i];
            const word patch =
            (
                region == "<outletRegion>" || outletRegions.size() == 1
              ? "outlet"
              : region
            );

            os_ << indent << region << " { name " << patch << "; }" << endl;
        }

        endDict(os_, false);
    }

    endDict(os_, false);
}


void Foam::snappyHexMeshConfiguration::writeSearchableBox(const label i)
{
    beginDict(os_, "box" + std::to_string(i));

    os_ << indent << "type searchableBox;" << nl
        << indent << "min " << refinementBoxes_[i].first() << ";" << nl
        << indent << "max " << refinementBoxes_[i].second() << ";" << endl;

    endDict(os_);
}


void Foam::snappyHexMeshConfiguration::writeSnappyGeometry()
{
    beginDict(os_, "geometry");

    forAll(surfaces_, i)
    {
        if (i != 0)
        {
            os_ << endl;
        }

        writeGeometrySurface(i);
    }

    forAll(refinementBoxes_, i)
    {
        os_ << endl;
        writeSearchableBox(i);
    }

    endDict(os_);
}


void Foam::snappyHexMeshConfiguration::writeFeatures()
{
    beginList(os_, "features");

    if (!implicitFeatures_)
    {
        forAll(surfaces_, i)
        {
            fileName eMeshFile(surfaces_[i].name() + ".eMesh");
            os_ << indent << "{ file " << eMeshFile
                << "; level 1; }" << endl;
        }
    }

    endList(os_);
}


void Foam::snappyHexMeshConfiguration::writeRefinementSurfacesLevel
(
    const label rl
)
{
    os_ << indent << "level (" << rl << " " << rl << ");" << endl;
}


void Foam::snappyHexMeshConfiguration::writeRefinementSurfacesLevel()
{
    writeRefinementSurfacesLevel(refinementLevel_);
}


void Foam::snappyHexMeshConfiguration::writeRefinementSurfacesLevel
(
    const word& name
)
{
    label rl(refinementLevel_);

    forAll(surfaceLevels_, i)
    {
        if (surfaceLevels_[i].first() == name)
        {
            rl = surfaceLevels_[i].second();
            break;
        }
    }

    writeRefinementSurfacesLevel(rl);
}


void Foam::snappyHexMeshConfiguration::writePatchInfo
(
    const word& type,
    const word& group
)
{
    if (group.empty())
    {
        os_ << indent << "patchInfo { type " << type << "; }" << endl;
    }
    else
    {
        beginDict(os_, "patchInfo");
        os_ << indent << "type "<< type << ";" << endl;
        os_ << indent << "inGroups (" << group << ");" << endl;
        endDict(os_, false);
    }
}


void Foam::snappyHexMeshConfiguration::writeRefinementSurfacesRegion
(
    const word regionName,
    const List<word>& regions
)
{
    switch (regions.size())
    {
        case 0:
        {
            return;
        }
        case 1:
        {
            os_ << indent << regions[0] << endl;
            break;
        }
        default:
        {
            os_ << indent << "\"" << regionName << ".*\"" << endl;
        }
    }

    beginDict(os_);
    writeRefinementSurfacesLevel();

    word group(regionName);
    if (regions.size() == 1 && regions[0] == regionName)
    {
        group = "";
    }

    writePatchInfo("patch", group);

    endDict(os_, false);
}


void Foam::snappyHexMeshConfiguration::writeRefinementSurfacesRegions
(
    const wordList& inletRegions,
    const wordList& outletRegions
)
{
    if (inletRegions.empty() && outletRegions.empty())
    {
        return;
    }

    os_ << endl;
    beginDict(os_, "regions");
    writeRefinementSurfacesRegion("inlet", inletRegions);
    writeRefinementSurfacesRegion("outlet", outletRegions);
    endDict(os_, false);
}


void Foam::snappyHexMeshConfiguration::writeRefinementSurfaces()
{
    beginDict(os_, "refinementSurfaces");

    forAll(surfaces_, i)
    {
        if (i != 0)
        {
            os_ << endl;
        }

        const word& name = surfaces_[i].name();

        beginDict(os_, name);

        writeRefinementSurfacesLevel(name);

        const wordList& inletRegions = surfaces_[i].inletRegions();
        const wordList& outletRegions = surfaces_[i].outletRegions();

        switch (surfaces_[i].type())
        {
            case surfaceType::wall:
            {
                writePatchInfo("wall");
                writeRefinementSurfacesRegions(inletRegions, outletRegions);
                break;
            }

            case surfaceType::external:
            {
                writePatchInfo("wall", "externalWall");
                writeRefinementSurfacesRegions(inletRegions, outletRegions);
                break;
            }

            case surfaceType::cellZone:
            case surfaceType::rotatingZone:
            {
                os_ << indent << "faceZone "
                    << surfaces_[i].name() << ";" << nl
                    << indent << "cellZone "
                    << surfaces_[i].name() << ";" << nl
                    << indent << "mode inside;" << endl;

                break;
            }

            case surfaceType::baffle:
            {
                writePatchInfo("wall");

                os_ << indent << "faceZone "
                    << surfaces_[i].name() << ";" << nl
                    << indent << "faceType boundary;" << endl;

                break;
            }
        }

        endDict(os_, false);
    }

    endDict(os_);
}


void Foam::snappyHexMeshConfiguration::writeRefinementRegion
(
    const word& name,
    const label level
)
{
    beginDict(os_, name);

    os_ << indent << "mode    inside;" << nl
        << indent << "level   " << level << ";" << endl;

    endDict(os_, false);
}


void Foam::snappyHexMeshConfiguration::writeRefinementRegions()
{
    if
    (
        refinementRegions_.empty()
     && refinementBoxes_.empty()
     && refinementDists_.empty()
    )
    {
        os_ << indent
            << "// delete \"-disabled\" below to enable refinementRegions"
            << endl;

        beginDict(os_, "refinementRegions-disabled");
        writeRefinementRegion("<surface>", refinementLevel_);
        endDict(os_);
    }
    else
    {
        beginDict(os_, "refinementRegions");

        forAll(refinementRegions_, i)
        {
            writeRefinementRegion
            (
                refinementRegions_[i].first(),
                refinementRegions_[i].second()
            );
        }

        forAll(refinementBoxes_, i)
        {
            writeRefinementRegion
            (
                "box" + std::to_string(i),
                refinementBoxes_[i].third()
            );
        }

        forAll(refinementDists_, i)
        {
            beginDict(os_, refinementDists_[i].first());

            os_ << indent << "mode    distance;" << nl
                << indent << "levels   (("
                << refinementDists_[i].second() << " "
                << refinementDists_[i].third() << "));" << endl;

            endDict(os_, false);
        }

        endDict(os_);
    }
}


void Foam::snappyHexMeshConfiguration::writeCastellatedMeshControls()
{
    beginDict(os_, "castellatedMeshControls");

    writeFeatures();
    writeRefinementSurfaces();
    writeRefinementRegions();

    // Needs customising
    os_ << indent << "insidePoint "
        << insidePoint_ << ";" << endl;

    os_ << indent << "nCellsBetweenLevels "
        << nCellsBetweenLevels_
        << ";" << endl;

    endDict(os_);
}


void Foam::snappyHexMeshConfiguration::writeSnapControls()
{
    beginDict(os_, "snapControls");

    os_ << indent << "explicitFeatureSnap    "
        << (implicitFeatures_ ? "off" : "on") << ";" << endl;
    os_ << indent << "implicitFeatureSnap    "
        << (implicitFeatures_ ? "on" : "off") << ";" << endl;

    endDict(os_);
}


void Foam::snappyHexMeshConfiguration::writeAddLayersControls()
{
    if (layers_ == 0)
    {
        return;
    }

    beginDict(os_, "addLayersControls");

    beginDict(os_, "layers");

    forAll(surfaces_, i)
    {
        switch (surfaces_[i].type())
        {
            case surfaceType::wall:
            case surfaceType::external:
            case surfaceType::baffle:
            {
                os_ << indent << "\"" << surfaces_[i].name()
                    << ".*\" { nSurfaceLayers "
                    << layers_ << "; }" << endl;
                break;
            }

            default: break;
        }
    }

    endDict(os_);

    os_ << indent << "relativeSizes       on; "
        << "// off, usually with firstLayerThickness" << nl
        << indent << "expansionRatio      1.2;" << nl
        << indent << "finalLayerThickness 0.5;" << nl
        << indent << "minThickness        1e-3;" << nl
        << indent << "firstLayerThickness-disabled 0.01;" << nl << nl
        << indent << "maxThicknessToMedialRatio-disabled 0.3;" << endl;

    endDict(os_);
}


void Foam::snappyHexMeshConfiguration::writeWriteFlags()
{
    os_ << "// delete \"-disabled\" to output mesh data, e.g. for layers"
        << endl;

    beginList(os_, "writeFlags-disabled");
    os_ << indent << "scalarLevels" << endl;
    os_ << indent << "layerSets" << endl;
    os_ << indent << "layerFields" << endl;
    endList(os_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::snappyHexMeshConfiguration::snappyHexMeshConfiguration
(
    const fileName& name,
    const fileName& dir,
    const Time& time,
    const meshingSurfaceList& surfaces,
    const label refinementLevel,
    const List<Tuple2<word, label>>& surfaceLevels,
    const List<Tuple2<word, label>>& refinementRegions,
    const List<Tuple3<vector, vector, label>>& refinementBoxes,
    const List<Tuple3<word, scalar, label>>& refinementDists,
    const bool implicitFeatures,
    const label layers,
    const point& insidePoint,
    const label nCellsBetweenLevels
)
:
    caseFileConfiguration(name, dir, time),
    surfaces_(surfaces),
    refinementLevel_(refinementLevel),
    surfaceLevels_(surfaceLevels),
    refinementRegions_(refinementRegions),
    refinementBoxes_(refinementBoxes),
    refinementDists_(refinementDists),
    implicitFeatures_(implicitFeatures),
    layers_(layers),
    insidePoint_(insidePoint),
    nCellsBetweenLevels_(nCellsBetweenLevels)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::snappyHexMeshConfiguration::~snappyHexMeshConfiguration()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::snappyHexMeshConfiguration::write()
{
    dict_.writeHeader(os_, word("dictionary"));
    os_ << "#includeEtc \"caseDicts/mesh/generation/snappyHexMeshDict.cfg\""
        << nl << endl;
    writeSnappySwitches();
    writeSnappyGeometry();
    writeCastellatedMeshControls();
    writeSnapControls();
    writeAddLayersControls();
    writeWriteFlags();

    os_ << "mergeTolerance 1e-6;";

    dict_.writeEndDivider(os_);
}


// ************************************************************************* //
