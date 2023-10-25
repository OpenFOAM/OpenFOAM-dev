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

Description
    Functions for checking mesh topology and geometry

SourceFiles
    checkTopology.C
    checkGeometry.C
    checkMesh.C

\*---------------------------------------------------------------------------*/

#ifndef meshCheck_H
#define meshCheck_H

#include "primitiveMeshCheck.H"
#include "polyMeshCheck.H"
#include "wedgePolyPatch.H"
#include "mergeAndWrite.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace meshCheck
{
    //- Find wedge with opposite orientation. Note: does not actually check
    //  that it is opposite, only that it has opposite normal and same axis.
    label findOppositeWedge(const polyMesh&, const wedgePolyPatch&);

    //- Check wedge orientation
    bool checkWedges
    (
        const polyMesh&,
        const bool report,
        const Vector<label>&,
        labelHashSet*
    );

    //- Check 0th vertex on coupled faces
    bool checkCoupledPoints(const polyMesh&, const bool report, labelHashSet*);

    //- Check the topology
    label checkTopology
    (
        const polyMesh& mesh,
        const bool allTopology,
        const autoPtr<surfaceWriter>& surfWriter,
        const autoPtr<setWriter>& setWriter
    );

    //- Check the geometry
    label checkGeometry
    (
        const polyMesh& mesh,
        const bool allGeometry,
        const scalar nonOrthThreshold,
        const scalar skewThreshold,
        const autoPtr<surfaceWriter>&,
        const autoPtr<setWriter>&
    );

    //- Check (subset of mesh including baffles) with mesh settings in dict.
    //  Collects incorrect faces in set.
    //  Returns true if one or more faces in error.
    bool checkMesh
    (
        const bool report,
        const polyMesh& mesh,
        const dictionary& dict,
        const labelList& checkFaces,
        const List<labelPair>& baffles,
        labelHashSet& wrongFaces
    );

    //- Check (subset of mesh) with mesh settings in dict.
    //  Collects incorrect faces in set.
    //  Returns true if one or more faces in error.
    bool checkMesh
    (
        const bool report,
        const polyMesh& mesh,
        const dictionary& dict,
        const labelList& checkFaces,
        labelHashSet& wrongFaces
    );

    //- Check mesh with mesh settings in dict.
    //  Collects incorrect faces in set.
    //  Returns true if one or more faces in error.
    bool checkMesh
    (
        const bool report,
        const polyMesh& mesh,
        const dictionary& dict,
        labelHashSet& wrongFaces
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
