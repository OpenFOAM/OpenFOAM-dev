/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::MPLICcellStorage

Description
    Provides local cell addressing for geometry and data for MPLIC class.

SourceFiles
    MPLICcellStorage.C

\*---------------------------------------------------------------------------*/

#ifndef MPLICcellStorage_H
#define MPLICcellStorage_H

#include "primitiveMesh.H"
#include "UIndirectList.H"
#include "uindirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MPLICcell Declaration
\*---------------------------------------------------------------------------*/

class MPLICcellStorage
{
    // Private data

        //- Reference to the mesh points
        const pointField& points_;

        //- Reference to the mesh faces
        const faceList& faces_;

        //- Reference to the mesh edges
        const edgeList& edges_;

        //- Reference to the mesh face edges
        const labelListList& edgeFaces_;

        //- Reference to the cell points
        const labelList& cPoints_;

        //- Reference to face list
        const labelList& cFaces_;

        //- Reference to cell edges
        const labelList& cEdges_;

        //- Alpha at the cell vertices
        const scalarField& pointsAlpha_;

        //- Velocity at the cell vertices
        const vectorField& pointsU_;

        //- Cell centre value of alpha
        const scalar cellAlpha_;

        //- Cell centre value of velocity
        const vector& celllU_;

        //- Owner? For all faces on the cell
        const boolList owns_;

        //- Cell volume
        const scalar volume_;

        //- Cell centre
        const vector& centre_;

        //- Face area vectors
        const UIndirectList<vector> Sf_;

        //- Face centres
        const UIndirectList<vector> Cf_;

        //- Face areas
        const UIndirectList<scalar> magSf_;

        //- Cell alpha min
        const scalar cellAlphaMin_;

        //- Cell alpha max
        const scalar cellAlphaMax_;

        //- Face alpha mins
        const scalarField facesAlphaMin_;

        //- Face alpha maxs
        const scalarField facesAlphaMax_;

        //- Global face edges
        const labelListList& faceEdges_;


    // Private Member Functions

        //- Return the face owner list for the faces of the cell
        boolList calcIsOwner
        (
            const primitiveMesh& mesh, const label celli
        ) const;

        //- Calculate minimum point alpha value in the cell
        scalar calcAlphaMin() const;

        //- Calculate maximum point alpha value in the cell
        scalar calcAlphaMax() const;

        //- Calculate minimum point alpha value on the cell faces
        scalarField calcFacesAlphaMin() const;

        //- Calculate maximum point alpha value on the cell faces
        scalarField calcFacesAlphaMax() const;


    public:

        // Constructors

            //- Construct from components
            MPLICcellStorage
            (
                const primitiveMesh& mesh,
                const scalarField& pointsAlpha,
                const vectorField& pointsU,
                const scalar cellAlpha,
                const vector& cellU,
                const label celli
            );


        // Public Member Functions

            //- Return reference to mesh points
            inline const pointField& points() const;

            //- Return reference to the mesh faces
            inline const faceList& faces() const;

            //- Return reference to the mesh edges
            inline const edgeList& edges() const;

            //- Return reference to the mesh face edges
            inline const labelListList& faceEdges() const;

            //- Return reference to the cell points
            inline const labelList& cellPoints() const;

            //- Return reference to face list
            inline const labelList& cellFaces() const;

            //- Return reference to cell edges
            inline const labelList& cellEdges() const;

            //- Return isOwners
            inline const boolList& isOwner() const;

            //- Return point alphas
            inline const scalarField& pointsAlpha() const;

            //- Return point velocities
            inline const vectorField& pointsU() const;

            //- Return cell alpha
            inline scalar cellAlpha() const;

            //- Return cell velocity
            inline const vector& cellU() const;

            //- Return cell volume
            inline scalar V() const;

            //- Return cell centre
            inline const vector& C() const;

            //- Return face surface area vectors
            inline const UIndirectList<vector>& Sf() const;

            //- Return face centres
            inline const UIndirectList<vector>& Cf() const;

            //- Return face areas
            inline const UIndirectList<scalar>& magSf() const;

            //- Return maximum point alpha value in the cell
            inline scalar cellAlphaMax() const;

            //- Return minimum point alpha value in the cell
            inline scalar cellAlphaMin() const;

            //- Return minimum point alpha value on the cell faces
            inline const scalarField& facesAlphaMin() const;

            //- Return maximum point alpha value on the cell faces
            inline const scalarField& facesAlphaMax() const;

            //- Return number of faces in the the cell
            inline label size() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MPLICcellStorageI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
