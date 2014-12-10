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
    Foam::polyAddCell

Description
    Class containing data for cell addition.

\*---------------------------------------------------------------------------*/

#ifndef polyAddCell_H
#define polyAddCell_H

#include "label.H"
#include "topoAction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyAddCell Declaration
\*---------------------------------------------------------------------------*/

class polyAddCell
:
    public topoAction
{
    // Private data

        //- Master point ID for cells blown up from points
        label masterPointID_;

        //- Master edge ID for cells blown up from edges
        label masterEdgeID_;

        //- Master face ID for cells blown up from faces
        label masterFaceID_;

        //- Master cell ID for cells blown up from cells
        label masterCellID_;

        //- Cell zone ID
        label zoneID_;


public:

    // Static data members

        //- Runtime type information
        TypeName("addCell");


    // Constructors

        //- Construct null.  Used for constructing lists
        polyAddCell()
        :
            masterPointID_(-1),
            masterEdgeID_(-1),
            masterFaceID_(-1),
            masterCellID_(-1),
            zoneID_(-1)
        {}

        //- Construct from components
        polyAddCell
        (
            const label masterPointID,
            const label masterEdgeID,
            const label masterFaceID,
            const label masterCellID,
            const label zoneID
        )
        :
            masterPointID_(masterPointID),
            masterEdgeID_(masterEdgeID),
            masterFaceID_(masterFaceID),
            masterCellID_(masterCellID),
            zoneID_(zoneID)
        {}

        //- Construct and return a clone
        virtual autoPtr<topoAction> clone() const
        {
            return autoPtr<topoAction>(new polyAddCell(*this));
        }


    // Default Destructor


    // Member Functions

        //- Is the cell mastered by a point
        bool isPointMaster() const
        {
            return masterPointID_ >= 0;
        }

        //- Is the cell mastered by an edge
        bool isEdgeMaster() const
        {
            return masterEdgeID_ >= 0;
        }

        //- Is the cell mastered by another face
        bool isFaceMaster() const
        {
            return masterFaceID_ >= 0;
        }

        //- Is the cell mastered by another cell
        bool isCellMaster() const
        {
            return masterCellID_ >= 0;
        }

        //- Is the cell appended with no master
        bool appended() const
        {
            return
                !isPointMaster() && !isEdgeMaster()
             && !isFaceMaster() && !isCellMaster();
        }

        //- Return master point ID
        label masterPointID() const
        {
            return masterPointID_;
        }

        //- Return master edge ID
        label masterEdgeID() const
        {
            return masterEdgeID_;
        }

        //- Return master face ID
        label masterFaceID() const
        {
            return masterFaceID_;
        }

        //- Return master cell ID
        label masterCellID() const
        {
            return masterCellID_;
        }

        //- Does the cell belong to a zone?
        bool isInZone() const
        {
            return zoneID_ >= 0;
        }

        //- Cell zone ID
        label zoneID() const
        {
            return zoneID_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
