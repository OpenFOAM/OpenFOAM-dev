/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    Foam::externalDisplacementMeshMover

Description
    Virtual base class for mesh movers with externally provided displacement
    field giving the boundary conditions. Move the mesh from the current
    location to a new location (so modify the mesh; v.s. motionSolver that
    only returns the new location).

    All mesh movers are expected to read the dictionary settings at invocation
    of move(), i.e. not cache any settings.

SourceFiles
    externalDisplacementMeshMover.C

\*---------------------------------------------------------------------------*/

#ifndef externalDisplacementMeshMover_H
#define externalDisplacementMeshMover_H

#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                Class externalDisplacementMeshMover Declaration
\*---------------------------------------------------------------------------*/

class externalDisplacementMeshMover
{
protected:

    // Protected data

        //- Baffles in the mesh
        List<labelPair> baffles_;

        //- Reference to point motion field
        pointVectorField& pointDisplacement_;


public:

    //- Runtime type information
    TypeName("externalDisplacementMeshMover");


    // Declare run-time New selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            externalDisplacementMeshMover,
            dictionary,
            (
                const dictionary& dict,
                const List<labelPair>& baffles,
                pointVectorField& pointDisplacement
            ),
            (dict, baffles, pointDisplacement)
        );


    // Constructors

        //- Construct from dictionary and displacement field. Dictionary is
        //  allowed to go out of scope!
        externalDisplacementMeshMover
        (
            const dictionary& dict,
            const List<labelPair>& baffles,
            pointVectorField& pointDisplacement
        );

        //- Disallow default bitwise copy construction
        externalDisplacementMeshMover
        (
            const externalDisplacementMeshMover&
        ) = delete;


    // Selectors

        //- Return a reference to the selected meshMover model
        static autoPtr<externalDisplacementMeshMover> New
        (
            const word& type,
            const dictionary& dict,
            const List<labelPair>& baffles,
            pointVectorField& pointDisplacement
        );


    //- Destructor
    virtual ~externalDisplacementMeshMover();


    // Member Functions

        // Access

            //- Return reference to the point motion displacement field
            pointVectorField& pointDisplacement()
            {
                return pointDisplacement_;
            }

            //- Return const reference to the point motion displacement field
            const pointVectorField& pointDisplacement() const
            {
                return pointDisplacement_;
            }

            const pointMesh& pMesh() const
            {
                return pointDisplacement_.mesh();
            }

            const polyMesh& mesh() const
            {
                return pMesh()();
            }


        // Mesh mover

            //- Move mesh using current pointDisplacement boundary values
            //  and current dictionary settings. Return true if successful
            //  (errors  on checkFaces less than allowable). Updates
            //  pointDisplacement.
            virtual bool move
            (
                const dictionary&,
                const label nAllowableErrors,
                const labelList& checkFaces
            ) = 0;

            //- Update local data for geometry changes
            virtual void movePoints(const pointField&);

            //-  Update local data for topology changes
            virtual void topoChange(const polyTopoChangeMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const externalDisplacementMeshMover&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
