/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::ensightPartFaces

Description
    An implementation of ensightPart to hold volume mesh faces.

SourceFiles
    ensightPartFaces.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPartFaces_H
#define ensightPartFaces_H

#include "ensightPart.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ensightPartFaces Declaration
\*---------------------------------------------------------------------------*/

class ensightPartFaces
:
    public ensightPart
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ensightPartFaces&);

        //- Track points used
        virtual localPoints calcLocalPoints() const;

        //- Element connectivity
        virtual void writeConnectivity
        (
            ensightGeoFile&,
            const word& key,
            const labelUList& idList,
            const labelUList& pointMap
        ) const;


protected:

        //- Addressable ensight element types
        enum elemType
        {
            tria3Elements,
            quad4Elements,
            nsidedElements
        };


    // Static data members

        static const List<word> elemTypes_;


    // Protected data

        //- Faces referenced
        const faceList& faces_;

        //- Can skip local point renumbering when points are contiguous
        const bool contiguousPoints_;


    // Protected Member Functions

        //- Classify the face shapes, set elemLists.
        void classify(const faceList&);

        //- Helper: write connectivity
        void writeConnectivity
        (
            ensightGeoFile&,
            const word& key,
            const faceList&,
            const labelUList& idList,
            const labelUList& pointMap
        ) const;


public:

    //- Runtime type information
    TypeName("ensightFaces");

    // Constructors

        //- Construct empty part with number and description
        ensightPartFaces(label partNumber, const string& partDescription);

        //- Construct part with number, description, points and faces
        //  Can skip local point renumbering when points are contiguous
        ensightPartFaces
        (
            label partNumber,
            const string& partDescription,
            const pointField&,
            const faceList&,
            const bool contiguousPoints = false
        );

        //- Construct from polyMesh and polyPatch
        ensightPartFaces
        (
            label partNumber,
            const polyMesh&,
            const polyPatch&
        );

        //- Construct as copy
        ensightPartFaces(const ensightPartFaces&);

        //- Reconstruct part characteristics (eg, element types) from Istream
        //  A part reconstructed in this manner can be used when writing fields,
        //  but cannot be used to write a new geometry
        //  \sa Foam::ensightPart::reconstruct
        ensightPartFaces(Istream&);

        //- Reconstruct part characteristics on freestore from Istream
        static autoPtr<ensightPartFaces> New(Istream& is)
        {
            return autoPtr<ensightPartFaces>(new ensightPartFaces(is));
        }


    //- Destructor
    virtual ~ensightPartFaces();


    // Member Functions

        //- Write geometry
        virtual void writeGeometry(ensightGeoFile&) const;

        //- Static listing of the element types
        virtual const List<word>& elementTypes() const
        {
            return elemTypes_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
