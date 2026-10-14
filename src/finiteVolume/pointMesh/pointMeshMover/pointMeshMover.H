/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::pointMeshMover

Description
    Abstract base class for pointMesh movers.

SourceFiles
    pointMeshMover.C

\*---------------------------------------------------------------------------*/

#ifndef pointMeshMover_H
#define pointMeshMover_H

#include "pointField.H"
#include "Time.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class polyMesh;
class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;

/*---------------------------------------------------------------------------*\
                        Class pointMeshMover Declaration
\*---------------------------------------------------------------------------*/

class pointMeshMover
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;


public:

    //- Runtime type information
    TypeName("pointMeshMover");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointMeshMover,
            dictionary,
            (
                const polyMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Selectors

        //- Select constructed from polyMesh and dictionary
        static autoPtr<pointMeshMover> New
        (
            const polyMesh&,
            const dictionary&
        );


    // Constructors

        //- Construct from polyMesh and dictionary and type.
        pointMeshMover(const polyMesh& mesh, const word& type);

        //- Clone function
        virtual autoPtr<pointMeshMover> clone() const;


    //- Destructor
    virtual ~pointMeshMover();


    // Member Functions

        //- Return reference to mesh
        const polyMesh& poly() const
        {
            return mesh_;
        }

        //- Provide new points for motion.  Solves for motion
        virtual tmp<pointField> newPoints() = 0;

        //- Correct point field for reduced-dimensionality cases
        virtual void twoDCorrectPoints(pointField&) const;

        //- Is the motion solid body? I.e., are the volumes and area magnitudes
        //  unchanged? Defaults to false. Set to true by
        //  pointMeshMovers::solidBodyMotion.
        virtual bool solidBodyMotion() const
        {
            return false;
        }

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&) = 0;

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&) = 0;

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&) = 0;

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&) = 0;

        //- Optionally write motion state information for restart
        virtual bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
