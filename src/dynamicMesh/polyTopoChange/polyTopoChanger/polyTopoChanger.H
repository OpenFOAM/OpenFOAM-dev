/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::polyTopoChanger

Description
    List of mesh modifiers defining the mesh dynamics.

SourceFiles
    polyTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef polyTopoChanger_H
#define polyTopoChanger_H

#include "regIOobject.H"
#include "PtrList.H"
#include "polyMeshModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMesh;
class mapPolyMesh;
class polyBoundaryMesh;


// Forward declaration of friend functions and operators

class polyTopoChanger;

Ostream& operator<<(Ostream&, const polyTopoChanger&);


/*---------------------------------------------------------------------------*\
                      Class polyTopoChanger Declaration
\*---------------------------------------------------------------------------*/

class polyTopoChanger
:
    public PtrList<polyMeshModifier>,
    public regIOobject
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        polyTopoChanger(const polyTopoChanger&);

        //- Disallow default bitwise assignment
        void operator=(const polyTopoChanger&);

        void readModifiers();


protected:

    // Protected data

        //- Reference to mesh
        polyMesh& mesh_;

public:

    //- Runtime type information
    TypeName("polyTopoChanger");


    // Constructors

        //- Read constructor given IOobject and a polyMesh
        polyTopoChanger(const IOobject&, polyMesh&);

        //- Read constructor for given polyMesh
        explicit polyTopoChanger(polyMesh&);


    //- Destructor
    virtual ~polyTopoChanger()
    {}


    // Member functions

        //- Return the mesh reference
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Return a list of patch types
        wordList types() const;

        //- Return a list of patch names
        wordList names() const;

        //- Is topology change required
        bool changeTopology() const;

        //- Return topology change request
        autoPtr<polyTopoChange> topoChangeRequest() const;

        //- Modify point motion
        void modifyMotionPoints(pointField&) const;

        autoPtr<mapPolyMesh> changeMesh
        (
            const bool inflate,
            const bool syncParallel = true,
            const bool orderCells = false,
            const bool orderPoints = false
        );

        //- Force recalculation of locally stored data on topological change
        void update(const mapPolyMesh& m);

        //- Add given set of topology modifiers to the topoChanger
        void addTopologyModifiers(const List<polyMeshModifier*>& tm);

        //- Find modifier given a name
        label findModifierID(const word& modName) const;


        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Member Operators

        bool operator!=(const polyTopoChanger&) const;
        bool operator==(const polyTopoChanger&) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const polyTopoChanger&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
