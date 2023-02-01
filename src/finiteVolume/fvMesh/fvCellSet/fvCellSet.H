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
    General cell set selection class for models that apply to sub-sets
    of the mesh.

    Currently supports cell selection from a set of points, a specified cellSet
    or cellZone or all of the cells.  The selection method can either be
    specified explicitly using the \c select entry or inferred from the presence
    of either a \c cellSet, \c cellZone or \c points entry.  The \c select entry
    is required to select \c all cells.

Usage
    Examples:
    \verbatim
        // Apply everywhere
        select          all;

        // Apply within a given cellSet
        select          cellSet; // Optional
        cellSet         rotor;

        // Apply within a given cellZone
        select          cellZone; // Optional
        cellSet         rotor;

        // Apply in cells containing a list of points
        select          points; // Optional
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

        const fvMesh& mesh_;

        //- Sum of cell volumes
        mutable scalar V_;


    // Private functions

        //- Set the sum of scalar volumes
        void setV();


public:

    // Constructors

        //- Construct from mesh and dictionary
        fvCellSet
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


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
