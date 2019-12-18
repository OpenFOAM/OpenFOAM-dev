/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    mechanism in faceZoneMesh, where all the zoned faces are registered
    with their zone number.

SourceFiles
    faceZone.C
    faceZoneNew.C

\*---------------------------------------------------------------------------*/

#ifndef faceZone_H
#define faceZone_H

#include "zone.H"
#include "faceZoneMeshFwd.H"
#include "boolList.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

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
        const faceZoneMesh& zoneMesh_;


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
        static const char * const labelsName;


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
                const label index,
                const faceZoneMesh& zm
            ),
            (name, dict, index, zm)
        );


    // Constructors

        //- Construct from components
        faceZone
        (
            const word& name,
            const labelUList& addr,
            const boolList& fm,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from components, moving contents
        faceZone
        (
            const word& name,
            labelList&& addr,
            boolList&& fm,
            const label index,
            const faceZoneMesh&
        );

        //- Construct from dictionary
        faceZone
        (
            const word& name,
            const dictionary&,
            const label index,
            const faceZoneMesh&
        );

        //- Construct given the original zone and resetting the
        //  face list and zone mesh information
        faceZone
        (
            const faceZone&,
            const labelUList& addr,
            const boolList& fm,
            const label index,
            const faceZoneMesh&
        );

        //- Construct given the original zone, resetting the
        //  face list and zone mesh information
        faceZone
        (
            const faceZone&,
            labelList&& addr,
            boolList&& fm,
            const label index,
            const faceZoneMesh&
        );

        //- Disallow default bitwise copy construction
        faceZone(const faceZone&) = delete;


        //- Construct and return a clone, resetting the zone mesh
        virtual autoPtr<faceZone> clone(const faceZoneMesh& zm) const
        {
            return autoPtr<faceZone>
            (
                new faceZone(*this, *this, flipMap(), index(), zm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and zone mesh
        virtual autoPtr<faceZone> clone
        (
            const labelUList& addr,
            const boolList& fm,
            const label index,
            const faceZoneMesh& zm
        ) const
        {
            return autoPtr<faceZone>
            (
                new faceZone(*this, addr, fm, index, zm)
            );
        }


    // Selectors

        //- Return a pointer to a new face zone
        //  created on freestore from dictionary
        static autoPtr<faceZone> New
        (
            const word& name,
            const dictionary&,
            const label index,
            const faceZoneMesh&
        );


    //- Destructor
    virtual ~faceZone();


    // Member Functions

        //- Return face flip map
        const boolList& flipMap() const
        {
            return flipMap_;
        }

        //- Helper function to re-direct to zone::localID(...)
        label whichFace(const label globalCellID) const;

        //- Return reference to primitive patch
        const primitiveFacePatch& operator()() const;

        //- Return zoneMesh reference
        const faceZoneMesh& zoneMesh() const;


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

        //- Correct patch after moving points
        virtual void movePoints(const pointField&);

        //- Update for changes in topology
        virtual void updateMesh(const mapPolyMesh&);

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
