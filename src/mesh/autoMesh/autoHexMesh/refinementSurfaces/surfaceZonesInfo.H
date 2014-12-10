/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Class
    Foam::surfaceZonesInfo

Description

SourceFiles
    surfaceZonesInfo.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceZonesInfo_H
#define surfaceZonesInfo_H

#include "NamedEnum.H"
#include "point.H"
#include "word.H"
#include "PtrList.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurface;
class searchableSurfaces;
class polyMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                      Class surfaceZonesInfo Declaration
\*---------------------------------------------------------------------------*/

class surfaceZonesInfo
{
public:

    //- Types of selection of area
    enum areaSelectionAlgo
    {
        INSIDE,
        OUTSIDE,
        INSIDEPOINT,
        NONE
    };

    static const NamedEnum<areaSelectionAlgo, 4> areaSelectionAlgoNames;

    //- What to do with faceZone faces
    enum faceZoneType
    {
        INTERNAL,
        BAFFLE,
        BOUNDARY
    };

    static const NamedEnum<faceZoneType, 3> faceZoneTypeNames;


private:

    // Private data

        //- Per 'interface' surface : name of faceZone to put faces into
        word faceZoneName_;

        //- Per 'interface' surface : name of cellZone to put cells into
        word cellZoneName_;

        //- Per 'interface' surface : (only used if surface is closed)
        //  How to select zone cells : surface inside or outside or given
        //  inside location.
        areaSelectionAlgo zoneInside_;

        //- If zoneInside=location gives the corresponding inside point
        point zoneInsidePoint_;

        //- Per 'interface' surface :
        //  What to do with outside
        faceZoneType faceType_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const surfaceZonesInfo&);


public:

    // Constructors

        //- Construct from surfaces and dictionary
        surfaceZonesInfo
        (
            const searchableSurface& surface,
            const dictionary& surfacesDict
        );

        //- Construct from components
        surfaceZonesInfo
        (
            const word& faceZoneNames,
            const word& cellZoneNames,
            const areaSelectionAlgo& zoneInside,
            const point& zoneInsidePoints,
            const faceZoneType& faceType
        );

        //- Copy constructor
        surfaceZonesInfo(const surfaceZonesInfo&);

        //- Return clone
        autoPtr<surfaceZonesInfo> clone() const
        {
            return autoPtr<surfaceZonesInfo>(new surfaceZonesInfo(*this));
        }


    // Member Functions

        // Access

            //- Per 'interface' surface : empty or name of faceZone to put
            //  faces into
            const word& faceZoneName() const
            {
                return faceZoneName_;
            }

            //- Per 'interface' surface : empty or name of cellZone to put
            //  cells into
            const word& cellZoneName() const
            {
                return cellZoneName_;
            }

            const areaSelectionAlgo& zoneInside() const
            {
                return zoneInside_;
            }

            //- Get specified inside locations for surfaces with a cellZone
            const point& zoneInsidePoint() const
            {
                return zoneInsidePoint_;
            }

            //- How to handle face of surfaces with a faceZone
            const faceZoneType& faceType() const
            {
                return faceType_;
            }


        // Query

            //- Get indices of unnamed surfaces (surfaces without faceZoneName)
            static labelList getUnnamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of named surfaces (surfaces with faceZoneName)
            static labelList getNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of surfaces with a cellZone that are closed and
            //  have 'inside' or 'outside' selection.
            static labelList getClosedNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const searchableSurfaces& allGeometry,
                const labelList& surfaces
            );

            //- Get indices of surfaces with a cellZone that are unclosed
            static labelList getUnclosedNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const searchableSurfaces& allGeometry,
                const labelList& surfaces
            );

            //- Get indices of surfaces with a cellZone that are closed.
            static labelList getAllClosedNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const searchableSurfaces& allGeometry,
                const labelList& surfaces
            );

            //- Get indices of surfaces with a cellZone that have 'insidePoint'
            //  section.
            static labelList getInsidePointNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            static labelList addCellZonesToMesh
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const labelList& namedSurfaces,
                polyMesh& mesh
            );

            static labelList addFaceZonesToMesh
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const labelList& namedSurfaces,
                polyMesh& mesh
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
