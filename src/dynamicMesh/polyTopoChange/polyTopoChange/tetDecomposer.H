/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::tetDecomposer

Description
    Decomposes polyMesh into tets.

SourceFiles
    tetDecomposer.C

\*---------------------------------------------------------------------------*/

#ifndef tetDecomposer_H
#define tetDecomposer_H

#include "DynamicList.H"
#include "PackedBoolList.H"
#include "boolList.H"
#include "typeInfo.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class polyTopoChange;
class face;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class tetDecomposer Declaration
\*---------------------------------------------------------------------------*/

class tetDecomposer
{
public:

    // Public data types

        enum decompositionType
        {
            FACE_CENTRE_TRIS, //- Faces decomposed into triangles
                              //  using face-centre

            FACE_DIAG_TRIS    //- Faces decomposed into triangles diagonally
        };
        static const NamedEnum<decompositionType, 2> decompositionTypeNames;


private:

    // Private data

        const polyMesh& mesh_;

        //- From cell to tet point
        labelList cellToPoint_;

        //- From face to tet point
        labelList faceToPoint_;


        // Per face, per point (faceCentre) or triangle (faceDiag)
        // the added tet on the owner side
        labelListList faceOwnerCells_;

        // Per face, per point (faceCentre) or triangle (faceDiag)
        // the added tet on the neighbour side
        labelListList faceNeighbourCells_;


    // Private Member Functions

        //- Modify a face
        void modifyFace
        (
            polyTopoChange& meshMod,
            const face& f,
            const label facei,
            const label own,
            const label nei,
            const label patchi,
            const label zoneI,
            const bool zoneFlip
        ) const;

        //- Add a face
        void addFace
        (
            polyTopoChange& meshMod,
            const face& f,
            const label own,
            const label nei,
            const label masterPointID,
            const label masterEdgeID,
            const label masterFaceID,
            const label patchi,
            const label zoneI,
            const bool zoneFlip
        ) const;

        //- Work out triangle index given the starting vertex in the face
        label triIndex(const label facei, const label fp) const;

        //- Disallow default bitwise copy construct
        tetDecomposer(const tetDecomposer&);

        //- Disallow default bitwise assignment
        void operator=(const tetDecomposer&);


public:

    //- Runtime type information
    ClassName("tetDecomposer");


    // Constructors

        //- Construct from mesh
        tetDecomposer(const polyMesh&);


    // Member Functions

        // Access

            //- From cell to tet point
            const labelList& cellToPoint() const
            {
                return cellToPoint_;
            }

            //- From face to tet point
            const labelList& faceToPoint() const
            {
                return faceToPoint_;
            }


            //- Per face, per point (faceCentre) or triangle (faceDiag)
            //  the added tet on the owner side
            const labelListList& faceOwnerCells() const
            {
                return faceOwnerCells_;
            }

            //- Per face, per point (faceCentre) or triangle (faceDiag)
            //  the added tet on the neighbour side
            const labelListList& faceNeighbourCells() const
            {
                return faceNeighbourCells_;
            }


        // Edit

            //- Insert all changes into meshMod to convert the polyMesh into
            //  tets.
            void setRefinement
            (
                const decompositionType decomposeType,
                polyTopoChange& meshMod
            );

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
