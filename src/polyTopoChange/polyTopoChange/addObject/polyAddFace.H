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

Class
    Foam::polyAddFace

Description
    A face addition data class. A face can be inflated either from a
    point or from another face and can either be in internal or a
    boundary face.

\*---------------------------------------------------------------------------*/

#ifndef polyAddFace_H
#define polyAddFace_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "label.H"
#include "face.H"
#include "topoAction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class polyAddFace Declaration
\*---------------------------------------------------------------------------*/

class polyAddFace
:
    public topoAction
{
    // Private Data

        //- Face identifier
        face face_;

        //- Face owner
        label owner_;

        //- Face neighbour
        label neighbour_;

        //- Master point ID for faces blown up from points
        label masterPointIndex_;

        //- Master edge ID for faces blown up from edges
        label masterEdgeIndex_;

        //- Master face ID for faces blown up from faces
        label masterFaceIndex_;

        //- Does the face flux need to be flipped
        bool flipFaceFlux_;

        //- Boundary patch ID
        label patchIndex_;

        //- Face zone ID
        label zoneIndex_;

        //- Face zone flip
        bool zoneFlip_;


public:

    // Static Data Members

        //- Runtime type information
        TypeName("addFace");


    // Constructors

        //- Construct null.  Used for constructing lists
        polyAddFace()
        :
            face_(0),
            owner_(-1),
            neighbour_(-1),
            masterPointIndex_(-1),
            masterEdgeIndex_(-1),
            masterFaceIndex_(-1),
            flipFaceFlux_(false),
            patchIndex_(-1),
            zoneIndex_(-1),
            zoneFlip_(false)
        {}


        //- Construct from components
        polyAddFace
        (
            const face& f,
            const label owner,
            const label neighbour,
            const label masterPointID,
            const label masterEdgeID,
            const label masterFaceID,
            const bool flipFaceFlux,
            const label patchID,
            const label zoneID,
            const bool zoneFlip
        )
        :
            face_(f),
            owner_(owner),
            neighbour_(neighbour),
            masterPointIndex_(masterPointID),
            masterEdgeIndex_(masterEdgeID),
            masterFaceIndex_(masterFaceID),
            flipFaceFlux_(flipFaceFlux),
            patchIndex_(patchID),
            zoneIndex_(zoneID),
            zoneFlip_(zoneFlip)
        {
            if (face_.size() < 3)
            {
                FatalErrorInFunction
                    << "This is not allowed.\n"
                    << "Face: " << face_
                    << " masterPointID:" << masterPointIndex_
                    << " masterEdgeID:" << masterEdgeIndex_
                    << " masterFaceID:" << masterFaceIndex_
                    << " patchID:" << patchIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << abort(FatalError);
            }

            if (min(face_) < 0)
            {
                FatalErrorInFunction
                    << "This is not allowed.\n"
                    << "Face: " << face_
                    << " masterPointID:" << masterPointIndex_
                    << " masterEdgeID:" << masterEdgeIndex_
                    << " masterFaceID:" << masterFaceIndex_
                    << " patchID:" << patchIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << abort(FatalError);
            }

            if (min(owner_, neighbour_) >= 0 && owner_ == neighbour_)
            {
                FatalErrorInFunction
                    << "This is not allowed.\n"
                    << "Face: " << face_
                    << " masterPointID:" << masterPointIndex_
                    << " masterEdgeID:" << masterEdgeIndex_
                    << " masterFaceID:" << masterFaceIndex_
                    << " patchID:" << patchIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << abort(FatalError);
            }

            if (neighbour_ >= 0 && patchID >= 0)
            {
                FatalErrorInFunction
                    << ".  This is not allowed.\n"
                    << "Face: " << face_
                    << " masterPointID:" << masterPointIndex_
                    << " masterEdgeID:" << masterEdgeIndex_
                    << " masterFaceID:" << masterFaceIndex_
                    << " patchID:" << patchIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << abort(FatalError);
            }

            if (owner_ < 0 && zoneID < 0)
            {
                FatalErrorInFunction
                    << "This is not allowed.\n"
                    << "Face: " << face_
                    << "Face: " << face_
                    << " masterPointID:" << masterPointIndex_
                    << " masterEdgeID:" << masterEdgeIndex_
                    << " masterFaceID:" << masterFaceIndex_
                    << " patchID:" << patchIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << abort(FatalError);
            }

            if (zoneIndex_ == -1 && zoneFlip)
            {
                FatalErrorInFunction
                    << "belong to zone.  This is not allowed.\n"
                    << "Face: " << face_
                    << " masterPointID:" << masterPointIndex_
                    << " masterEdgeID:" << masterEdgeIndex_
                    << " masterFaceID:" << masterFaceIndex_
                    << " patchID:" << patchIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << abort(FatalError);
            }
        }

        //- Construct and return a clone
        virtual autoPtr<topoAction> clone() const
        {
            return autoPtr<topoAction>(new polyAddFace(*this));
        }


    // Default Destructor

    // Member Functions

        //- Return face
        const face& newFace() const
        {
            return face_;
        }

        //- Return owner cell
        label owner() const
        {
            return owner_;
        }

        //- Return neighbour cell
        label neighbour() const
        {
            return neighbour_;
        }

        //- Is the face mastered by a point
        bool isPointMaster() const
        {
            return masterPointIndex_ >= 0;
        }

        //- Is the face mastered by an edge
        bool isEdgeMaster() const
        {
            return masterEdgeIndex_ >= 0;
        }

        //- Is the face mastered by another face
        bool isFaceMaster() const
        {
            return masterFaceIndex_ >= 0;
        }

        //- Is the face appended with no master
        bool appended() const
        {
            return !isPointMaster() && !isEdgeMaster() && !isFaceMaster();
        }

        //- Return master point ID
        label masterPointIndex() const
        {
            return masterPointIndex_;
        }

        //- Return master edge ID
        label masterEdgeIndex() const
        {
            return masterEdgeIndex_;
        }

        //- Return master face ID
        label masterFaceIndex() const
        {
            return masterFaceIndex_;
        }

        //- Does the face flux need to be flipped
        bool flipFaceFlux() const
        {
            return flipFaceFlux_;
        }

        //- Does the face belong to a boundary patch?
        bool isInPatch() const
        {
            return patchIndex_ >= 0;
        }

        //- Boundary patch ID
        label patchIndex() const
        {
            return patchIndex_;
        }

        //- Does the face belong to a zone?
        bool isInZone() const
        {
            return zoneIndex_ >= 0;
        }

        //- Is the face only a zone face (i.e. not belonging to a cell)
        bool onlyInZone() const
        {
            return zoneIndex_ >= 0 && owner_ < 0 && neighbour_ < 0;
        }

        //- Face zone ID
        label zoneIndex() const
        {
            return zoneIndex_;
        }

        //- Face zone flip
        label zoneFlip() const
        {
            return zoneFlip_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
