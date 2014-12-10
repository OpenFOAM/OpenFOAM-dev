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
    Foam::faceZoneSet

Description
    Like faceSet but updates faceZone when writing.

SourceFiles
    faceZone.C

\*---------------------------------------------------------------------------*/

#ifndef faceZoneSet_H
#define faceZoneSet_H

#include "faceSet.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceZoneSet Declaration
\*---------------------------------------------------------------------------*/

class faceZoneSet
:
    public faceSet
{
    // Private data

        const polyMesh& mesh_;

        labelList addressing_;

        boolList flipMap_;

   // Private Member Functions


public:

    //- Runtime type information
    TypeName("faceZoneSet");


    // Constructors

        //- Construct from objectRegistry and name
        faceZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            readOption r=MUST_READ,
            writeOption w=NO_WRITE
        );

        //- Construct from additional size of labelHashSet
        faceZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            const label,
            writeOption w=NO_WRITE
        );

        //- Construct from existing set
        faceZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            const topoSet&,
            writeOption w=NO_WRITE
        );



    //- Destructor
    virtual ~faceZoneSet();


    // Member functions

        const labelList& addressing() const
        {
            return addressing_;
        }

        labelList& addressing()
        {
            return addressing_;
        }


        const boolList& flipMap() const
        {
            return flipMap_;
        }

        boolList& flipMap()
        {
            return flipMap_;
        }


        //- Sort addressing and make faceSet part consistent with addressing
        void updateSet();

        //- Invert contents. (insert all members 0..maxLen-1 which were not in
        //  set)
        virtual void invert(const label maxLen);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const topoSet& set);

        //- Add elements present in set.
        virtual void addSet(const topoSet& set);

        //- Delete elements present in set.
        virtual void deleteSet(const topoSet& set);

        //- Sync faceZoneSet across coupled patches.
        virtual void sync(const polyMesh& mesh);

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh&,
            const label maxLen
        ) const;

        //- Write faceZone
        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType
        ) const;

        //- Update any stored data for new labels
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;




};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
