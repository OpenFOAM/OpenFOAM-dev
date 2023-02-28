/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::polyCellSet

Description
    General run-time selected cell set selection class for polyMesh

    Currently supports cell selection from:
      - a set of points
      - a cellSet
      - a cellZone
      - all of the cells.

Usage
    Examples:
    \verbatim
        // Select all cells
        select   all;

        // Select the cells within the given cellSet
        select   cellSet;   // Optional
        cellSet  rotor;

        // Select the cells within the given cellZone
        select   cellZone;  // Optional
        cellZone rotor;

        // Select the cells containing a list of points
        select   points;    // Optional
        points
        (
            (2.25 0.5 0)
            (2.75 0.5 0)
        );
    \endverbatim

SourceFiles
    polyCellSet.C

\*---------------------------------------------------------------------------*/

#ifndef polyCellSet_H
#define polyCellSet_H

#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class polyMeshMap;
class polyDistributionMap;

/*---------------------------------------------------------------------------*\
                        Class polyCellSet Declaration
\*---------------------------------------------------------------------------*/

class polyCellSet
{
public:

    // Public data

        //- Enumeration for selection mode types
        enum class selectionTypes
        {
            points,
            cellSet,
            cellZone,
            all
        };

        //- Word list of selection type names
        static const NamedEnum<selectionTypes, 4> selectionTypeNames;


private:

    // Private data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Cell selection type
        selectionTypes selectionType_;

        //- Optional name of cell set or zone
        word cellSetName_;

        //- Optional list of points
        List<point> points_;

        //- Set of selected cells (not used for all or cellZone)
        mutable labelList cells_;


    // Private functions

        //- Set the cells
        void setCells();

        //- Return the identity map of length len
        labelUList identityMap(const label len) const;


public:

    // Constructors

        //- Construct from mesh. Will select all.
        polyCellSet(const polyMesh& mesh);

        //- Construct from mesh and dictionary
        polyCellSet(const polyMesh& mesh, const dictionary& dict);


    //- Destructor
    ~polyCellSet();


    // Member Functions

        // Access

            //- Return const access to the cell selection type
            inline const selectionTypes& selectionType() const;

            //- Return const access to the name of cellSet
            inline const word& cellSetName() const;

            //- Return const access to the cell set
            inline labelUList cells() const;

            //- Return the number of cells in the set
            inline label nCells() const;

            //- Return true if the set comprises all the cells
            inline bool all() const;

            //- Return the cell index corresponding to the cell set index
            inline label celli(const label i) const;


        // Mesh changes

            //- Update for mesh motion
            void movePoints();

            //- Update topology using the given map
            void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            void distribute(const polyDistributionMap&);


        // IO

            //- Read coefficients dictionary
            bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyCellSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
