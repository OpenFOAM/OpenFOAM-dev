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

Class
    Foam::cellLooper

Description
    Abstract base class. Concrete implementations know how to cut a cell
    (i.e. determine a loop around the circumference).

    Loop around the cell is given as the vertices to be cut and edges to
    be cut (and a weight between 0 and 1 giving where the cut is to be
    made). Main routine is 'cut' which gets called for every cell and
    gets the current cut situation and expects to return a loop on the
    cell circumference.

    Calling function needs to determine whether cellLooper is compatible with
    existing set of cuts.

    Also contains various utility functions which implementations might want to
    use.

SourceFiles
    cellLooper.C

\*---------------------------------------------------------------------------*/

#ifndef cellLooper_H
#define cellLooper_H

#include "edgeVertex.H"
#include "vector.H"
#include "boolList.H"
#include "scalarField.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class plane;

/*---------------------------------------------------------------------------*\
                           Class cellLooper Declaration
\*---------------------------------------------------------------------------*/

class cellLooper
:
    public edgeVertex
{
    // Private data


protected:

    // Protected Member Functions

        //- Get faces (on cell) connected to vertI which are not using edgeI
        labelList getVertFacesNonEdge
        (
            const label cellI,
            const label edgeI,
            const label vertI
        ) const;

        //- Get first edge connected to vertI and on faceI
        label getFirstVertEdge
        (
            const label faceI,
            const label vertI
        ) const;

        //- Get edges (on cell) connected to vertI which are not on faceI
        labelList getVertEdgesNonFace
        (
            const label cellI,
            const label faceI,
            const label vertI
        ) const;

        //- Return edge from cellEdges that is most perpendicular
        //  to refinement direction.
        label getMisAlignedEdge(const vector& refDir, const label cellI) const;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellLooper(const cellLooper&);

        //- Disallow default bitwise assignment
        void operator=(const cellLooper&);


public:

    //- Runtime type information
    TypeName("cellLooper");


    // Declare run-time constructor selection table

        // For the direct constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            cellLooper,
            word,
            (
                const polyMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        cellLooper(const polyMesh& mesh);

        //- Clone
        autoPtr<cellLooper> clone() const
        {
            notImplemented("autoPtr<tcellLooper> clone() const");
            return autoPtr<cellLooper>(NULL);
        }


    // Selectors

        //- Return a reference to the selected cellLooper
        static autoPtr<cellLooper> New
        (
            const word& type,
            const polyMesh& mesh
        );


    //- Destructor
    virtual ~cellLooper();


    // Member Functions

        //- Create cut along circumference of cellI. Gets current mesh cuts
        //  vertIsCut, edgeIsCut, edgeWeight).
        //  Cut along circumference is expressed as cellVertCut,
        //  cellEdgeToWeight. Returns true if succesfull. Still might not
        //  be compatible with existing cuts but this should be handled by
        //  caller).
        virtual bool cut
        (
            const vector& refDir,
            const label cellI,
            const boolList& vertIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,

            labelList& loop,
            scalarField& loopWeights
        ) const = 0;

        //- Same but now also base point of cut provided (instead of always
        //  cell centre)
        virtual bool cut
        (
            const plane& cutPlane,
            const label cellI,
            const boolList& vertIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,

            labelList& loop,
            scalarField& loopWeights
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
