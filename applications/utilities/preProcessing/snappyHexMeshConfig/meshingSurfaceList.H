/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::meshingSurfaceList

Description
    List of meshingSurfaces which stores the overall bounding box of all the
    meshingSurfaces.

SourceFiles
    meshingSurfaceList.C

\*---------------------------------------------------------------------------*/

#ifndef meshingSurfaceList_H
#define meshingSurfaceList_H

#include "meshingSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class meshingSurfaceList Declaration
\*---------------------------------------------------------------------------*/

class meshingSurfaceList
:
    public PtrList<meshingSurface>
{
    // Private Typedefs

        //- The surface type
        typedef meshingSurface::surfaceType surfaceType;


    // Private Data

        //- Bounding box for all the surfaces
        boundBox bb_;

        //- Bounding box for the rotatingZone surfaces
        boundBox rzbb_;


    // Private Member Functions

        //- Merge the second bounding box into the first
        void mergeBoundingBoxes
        (
            boundBox& bb1,
            const boundBox& bb2
        );

        //- Swap surfaces list to make external surface index 0
        void swapExternalIndexZero(const label index);

        //- Return true if specified regions > 0 and match surface regions
        bool regionsValid
        (
            const wordList& specifiedRegions,
            const wordList& regions,
            const word& opt
        );

        //- Set the specified surface types
        void setSurfaceTypes
        (
            const List<word>& surfaces,
            const surfaceType& type
        );

        //- Identify cell zones by analysing the surfaces
        void identifyCellZones();

        //- Write a word list as comma separated entries
        void reportWordList(const wordList& wl);

        //- Write a report on the surfaces to stdout
        void reportSurfaces();

        //- Set bounds specified in boundBox
        void setBounds(const boundBox& bb);

        //- Set bounds of rotatingZone surfaces
        void setRotatingZoneBounds();

public:

    // Constructors

        //- Construct from components
        meshingSurfaceList
        (
            const Time& time,
            const fileNameList& surfaces,
            const wordList& cellZones,
            const wordList& rotatingZones,
            const wordList& baffles,
            const boundBox& bb,
            const wordList& inletRegions,
            const wordList& outletRegions
        );

        //- Disallow default bitwise copy construction
        meshingSurfaceList(const meshingSurfaceList&) = delete;


    //- Destructor
    ~meshingSurfaceList();


    // Member Functions

        // Bounding box
        const boundBox& bb() const
        {
            return bb_;
        }

        // Bounding box of rotating zone surfaces
        const boundBox& rzbb() const
        {
            return rzbb_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshingSurfaceList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
