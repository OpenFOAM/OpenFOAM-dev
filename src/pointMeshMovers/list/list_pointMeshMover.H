/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2026 OpenFOAM Foundation
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
    Foam::pointMeshMovers::List

Description
    Motion of the mesh specified as a list of pointMeshMovers.

    The motion solvers are executed in order and the resulting displacements
    accumulated into an overall displacement and the displaced point positions
    returned.

SourceFiles
    list_pointMeshMover.C

\*---------------------------------------------------------------------------*/

#ifndef list_pointMeshMover_H
#define list_pointMeshMover_H

#include "pointMeshMover.H"
#include "PtrListDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
                            Class list Declaration
\*---------------------------------------------------------------------------*/

class list
:
    public pointMeshMover
{
    // Private Data

        PtrListDictionary<pointMeshMover> movers_;


public:

    //- Runtime type information
    TypeName("list");


    // Constructors

        //- Construct from mesh and dictionary
        list(const polyMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        list(const list&) = delete;


    //- Destructor
    virtual ~list();


    // Member Functions

        //- Provide current points for motion.  Uses current motion field
        virtual tmp<pointField> newPoints();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const list&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
