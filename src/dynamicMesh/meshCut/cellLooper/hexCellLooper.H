/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::hexCellLooper

Description
    Implementation of cellLooper.

    This one walks hexes in a topological way:
      - cross edge to other face
      - cross face by walking edge-point-edge across to reach the other side.
    (edges are always cut through the middle)

    For anything else (tet, prism, .. poly) it will use geomCellLooper
    (which does a purely geometric cut using a plane through cell centre)

SourceFiles
    hexCellLooper.C

\*---------------------------------------------------------------------------*/

#ifndef hexCellLooper_H
#define hexCellLooper_H

#include "geomCellLooper.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class cellModel;

/*---------------------------------------------------------------------------*\
                        Class hexCellLooper Declaration
\*---------------------------------------------------------------------------*/

class hexCellLooper
:
    public geomCellLooper
{

protected:

    // Protected data

        //- Reference to hex cell shape.
        const cellModel& hex_;

private:

    // Private Member Functions

        //- Walk across faces of hex. Update loop/loopWeights with edges cut.
        //  Return true if successful walk. (should be always!)
        bool walkHex
        (
            const label celli,
            const label startFacei,
            const label startEdgeI,

            labelList& loop,
            scalarField& loopWeights
        ) const;

        //- Convert loop into face and points
        void makeFace
        (
            const labelList& loop,
            const scalarField& loopWeights,

            labelList& faceVerts,
            pointField& facePoints
        ) const;


public:

    //- Runtime type information
    TypeName("hexCellLooper");


    // Constructors

        //- Construct from components
        hexCellLooper(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        hexCellLooper(const hexCellLooper&) = delete;


    //- Destructor
    virtual ~hexCellLooper();


    // Member Functions

        //- Create cut along circumference of celli. Gets current mesh cuts.
        //  Cut along circumference is expressed as loop of cuts plus weights
        //  for cuts along edges (only valid for edge cuts).
        //  Return true if successful cut.
        virtual bool cut
        (
            const vector& refDir,
            const label celli,
            const boolList& vertIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,

            labelList& loop,
            scalarField& loopWeights
        ) const;


        //- Same but now also base point of cut provided (instead of always
        //  cell centre)
        virtual bool cut
        (
            const plane& cutPlane,
            const label celli,
            const boolList& vertIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,

            labelList& loop,
            scalarField& loopWeights
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const hexCellLooper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
