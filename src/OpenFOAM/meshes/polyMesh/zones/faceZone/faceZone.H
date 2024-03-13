/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::faceZone

Description
    A subset of mesh faces organised as a primitive patch.

    For quick check whether a face belongs to the zone use the lookup
    mechanism in meshFaceZones, where all the zoned faces are registered
    with their zone number.

SourceFiles
    faceZone.C
    faceZoneNew.C

\*---------------------------------------------------------------------------*/

#ifndef faceZone_H
#define faceZone_H

#include "zone.H"
#include "meshFaceZonesFwd.H"
#include "boolList.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class faceZone;
Ostream& operator<<(Ostream&, const faceZone&);


/*---------------------------------------------------------------------------*\
                          Class faceZone Declaration
\*---------------------------------------------------------------------------*/

class faceZone
:
    public zone
{
    // Private Data

        //- The name associated with the zone-labels dictionary entry
        static const word labelsName_;


protected:

    // Protected data

        //- Flip map for all faces in the zone.  Set to true if the
        //  face needs to be flipped to achieve the correct orientation.
        boolList flipMap_;

        //- Reference to zone list
        const meshFaceZones& meshZones_;


        // Demand-driven private data

            //- Primitive patch made out of correctly flipped faces
            mutable primitiveFacePatch* patchPtr_;

            //- Master cell layer
            mutable labelList* masterCellsPtr_;

            //- Slave cell layer
            mutable labelList* slaveCellsPtr_;

            //- Global edge addressing
            mutable labelList* mePtr_;


    // Protected Member Functions

        //- Build primitive patch
        void calcFaceZonePatch() const;

        //- Return map of local face indices
        const Map<label>& faceLookupMap() const;

        //- Calculate master and slave face layer
        void calcCellLayers() const;

        //- Check addressing
        void checkAddressing() const;


public:

    // Static Data Members

        //- The name associated with the zone-labels dictionary entry
        static const char* const labelsName;


    //- Runtime type information
    TypeName("faceZone");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            faceZone,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const meshFaceZones& mz
            ),
            (name, dict, mz)
        );


    // Constructors

        //- Construct from components
        faceZone
        (
            const word& name,
            const labelUList& addr,
            const boolList& fm,
            const meshFaceZones& mz
        );

        //- Construct from components, moving contents
        faceZone
        (
            const word& name,
            labelList&& addr,
            boolList&& fm,
            const meshFaceZones&
        );

        //- Construct from dictionary
        faceZone
        (
            const word& name,
            const dictionary&,
            const meshFaceZones&
        );

        //- Construct given the original zone and resetting the
        //  face list and mesh zones information
        faceZone
        (
            const faceZone&,
            const labelUList& addr,
            const boolList& fm,
            const meshFaceZones&
        );

        //- Construct given the original zone, resetting the
        //  face list and mesh zones information
        faceZone
        (
            const faceZone&,
            labelList&& addr,
            boolList&& fm,
            const meshFaceZones&
        );

        //- Disallow default bitwise copy construction
        faceZone(const faceZone&) = delete;


        //- Construct and return a clone, resetting the mesh zones
        virtual autoPtr<faceZone> clone(const meshFaceZones& mz) const
        {
            return autoPtr<faceZone>
            (
                new faceZone(*this, *this, flipMap(), mz)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and mesh zones
        virtual autoPtr<faceZone> clone
        (
            const labelUList& addr,
            const boolList& fm,
            const meshFaceZones& mz
        ) const
        {
            return autoPtr<faceZone>
            (
                new faceZone(*this, addr, fm, mz)
            );
        }


    // Selectors

        //- Return a pointer to a new face zone
        //  created on freestore from dictionary
        static autoPtr<faceZone> New
        (
            const word& name,
            const dictionary&,
            const meshFaceZones&
        );


    //- Destructor
    virtual ~faceZone();


    // Member Functions

        //- Return face flip map
        const boolList& flipMap() const
        {
            return flipMap_;
        }

        //- Convenient renaming of zone::localIndex(globalIndex)
        label whichFace(const label globalIndex) const;

        //- Return reference to primitive patch
        const primitiveFacePatch& operator()() const;

        //- Return meshZones reference
        const meshFaceZones& meshZones() const;


        // Addressing into mesh

            //- Return labels of master cells (cells next to the master face
            //  zone in the prescribed direction)
            const labelList& masterCells() const;

            //- Return labels of slave cells
            const labelList& slaveCells() const;

            //- Return global edge index for local edges
            const labelList& meshEdges() const;


        //- Clear addressing
        virtual void clearAddressing();

        //- Reset addressing and flip map (clearing demand-driven data)
        virtual void resetAddressing(const labelUList&, const boolList&);

        //- Check zone definition. Return true if in error.
        virtual bool checkDefinition(const bool report = false) const;

        //- Check whether all procs have faces synchronised. Return
        //  true if in error.
        virtual bool checkParallelSync(const bool report = false) const;

        //- Swap two faceZones
        virtual void swap(faceZone&);

        //- Correct patch after moving points
        virtual void movePoints(const pointField&);

        //- Update zone using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;


    // Member Operators

        //- Assignment to zone, clearing demand-driven data
        void operator=(const faceZone&);

        //- Move assignment to zone, clearing demand-driven data
        void operator=(faceZone&&);


    // I-O

        //- Ostream Operator
        friend Ostream& operator<<(Ostream&, const faceZone&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
