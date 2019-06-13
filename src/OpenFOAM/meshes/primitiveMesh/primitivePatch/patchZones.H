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

Class
    Foam::patchZones

Description
    Calculates zone number for every face of patch.

    Gets constructed from patch and is a labeList with zone number for
    every patch face.

SourceFiles
    patchZones.C

\*---------------------------------------------------------------------------*/

#ifndef patchZones_H
#define patchZones_H

#include "labelList.H"
#include "pointField.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class patchZones Declaration
\*---------------------------------------------------------------------------*/

class patchZones
:
    public labelList
{
    // Private Data

        //- Reference to patch
        const primitivePatch& pp_;

        //- Reference to protected edges
        const boolList& borderEdge_;

        //- Max number of zones
        label nZones_;


    // Private Member Functions

        //- Gets labels of changed faces and propagates them to the edges.
        //  Returns labels of edges changed. Does not cross edges marked in
        //  regionEdge
        labelList faceToEdge
        (
            const labelList& changedFaces,
            labelList& edgeRegion
        );

        //- Reverse of faceToEdge: gets edges and returns faces
        labelList edgeToFace(const labelList& changedEdges);

        //- Fill *this with current zone for every face reachable
        //  from facei without crossing edge marked in borderEdge.
        void markZone(label facei);

public:

    ClassName("patchZones");


    // Constructors

        //- Fills *this with zone of face. Zone is area
        //  reachable by edge crossing without crossing borderEdge
        //  (bool for every edge in patch).
        patchZones(const primitivePatch& pp, const boolList& borderEdge);


    // Member Functions

        //- Number of zones
        label nZones() const
        {
            return nZones_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
