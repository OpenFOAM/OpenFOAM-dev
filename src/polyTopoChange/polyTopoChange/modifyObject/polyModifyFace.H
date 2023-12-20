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
    Foam::polyModifyFace

Description
    Class describing modification of a face.


\*---------------------------------------------------------------------------*/

#ifndef polyModifyFace_H
#define polyModifyFace_H

#include "label.H"
#include "face.H"
#include "topoAction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class polyModifyFace Declaration
\*---------------------------------------------------------------------------*/

class polyModifyFace
:
    public topoAction
{
    // Private Data

        //- Face
        face face_;

        //- Master face ID
        label faceIndex_;

        //- Face owner
        label owner_;

        //- Face neighbour
        label neighbour_;

        //- Does the face flux need to be flipped
        bool flipFaceFlux_;

        //- Boundary patch ID
        label patchIndex_;

        //- Remove from current zone
        bool removeFromZone_;

        //- Face zone ID
        label zoneIndex_;

        //- Face zone flip
        bool zoneFlip_;


public:

    // Static Data Members

        //- Runtime type information
        TypeName("modifyFace");


    // Constructors

        //- Construct null.  Used in constructing lists
        polyModifyFace()
        :
            face_(0),
            faceIndex_(-1),
            owner_(-1),
            neighbour_(-1),
            flipFaceFlux_(false),
            patchIndex_(-1),
            removeFromZone_(false),
            zoneIndex_(-1),
            zoneFlip_(false)
        {}

        //- Construct from components
        polyModifyFace
        (
            const face& f,
            const label faceID,
            const label owner,
            const label neighbour,
            const bool flipFaceFlux,
            const label patchID,
            const bool removeFromZone,
            const label zoneID,
            const bool zoneFlip
        )
        :
            face_(f),
            faceIndex_(faceID),
            owner_(owner),
            neighbour_(neighbour),
            flipFaceFlux_(flipFaceFlux),
            patchIndex_(patchID),
            removeFromZone_(removeFromZone),
            zoneIndex_(zoneID),
            zoneFlip_(zoneFlip)
        {
            if (face_.size() < 3)
            {
                FatalErrorInFunction
                    << "Invalid face: less than 3 points. This is not allowed\n"
                    << "Face: " << face_
                    << " faceID:" << faceIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << abort(FatalError);
            }

            if (min(face_) < 0)
            {
                FatalErrorInFunction
                    << "This is not allowed.\n"
                    << " faceID:" << faceIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << abort(FatalError);
            }

            if (min(owner_, neighbour_) >= 0 && owner_ == neighbour_)
            {
                FatalErrorInFunction
                    << "This is not allowed.\n"
                    << "Face: " << face_
                    << " faceID:" << faceIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << abort(FatalError);
            }

            if (neighbour_ >= 0 && patchIndex_ >= 0)
            {
                FatalErrorInFunction
                    << "This is not allowed.\n"
                    << "Face: " << face_
                    << " faceID:" << faceIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << " patchID:" << patchIndex_
                    << abort(FatalError);
            }

            if (zoneIndex_ < 0 && zoneFlip )
            {
                FatalErrorInFunction
                    << "belong to zone.  This is not allowed.\n"
                    << "Face: " << face_
                    << " faceID:" << faceIndex_
                    << " owner:" << owner_
                    << " neighbour:" << neighbour_
                    << abort(FatalError);
            }
        }

        //- Construct and return a clone
        virtual autoPtr<topoAction> clone() const
        {
            return autoPtr<topoAction>(new polyModifyFace(*this));
        }


    // Default Destructor

    // Member Functions

        //- Return face
        const face& newFace() const
        {
            return face_;
        }

        //- Return master face ID
        label faceIndex() const
        {
            return faceIndex_;
        }

        //- Return owner cell ID
        label owner() const
        {
            return owner_;
        }

        //- Return owner cell ID
        label neighbour() const
        {
            return neighbour_;
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

        bool removeFromZone() const
        {
            return removeFromZone_;
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
