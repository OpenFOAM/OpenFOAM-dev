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
    Foam::fvCellSet

Description
    General run-time selected cell set selection class for fvMesh

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

See also
    Foam::polyCellSet

SourceFiles
    fvCellSet.C

\*---------------------------------------------------------------------------*/

#ifndef fvCellSet_H
#define fvCellSet_H

#include "polyCellSet.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class polyMeshMap;
class polyDistributionMap;

/*---------------------------------------------------------------------------*\
                        Class fvCellSet Declaration
\*---------------------------------------------------------------------------*/

class fvCellSet
:
    public polyCellSet
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Sum of cell volumes
        mutable scalar V_;


    // Private functions

        //- Set the sum of scalar volumes
        void setV();


protected:

    // Protected Member Functions

        //- Output file header information
        void writeFileHeader
        (
            const functionObjects::writeFile& wf,
            Ostream& file
        );


public:

    // Constructors

        //- Construct from mesh. Will select all.
        fvCellSet(const fvMesh& mesh);

        //- Construct from mesh and dictionary
        fvCellSet(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    ~fvCellSet();


    // Member Functions

        // Access

            //- Return const access to the total cell volume
            inline scalar V() const;


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

#include "fvCellSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
