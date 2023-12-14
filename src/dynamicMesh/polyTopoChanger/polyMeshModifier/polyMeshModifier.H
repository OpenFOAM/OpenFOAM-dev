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
    Foam::polyMeshModifier

Description
    Virtual base class for mesh modifiers.

SourceFiles
    polyMeshModifier.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshModifier_H
#define polyMeshModifier_H

#include "typeInfo.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class polyTopoChangeMap;


/*---------------------------------------------------------------------------*\
                      Class polyMeshModifier Declaration
\*---------------------------------------------------------------------------*/

class polyMeshModifier
{
    // Private Data

        //- Name of modifier
        word name_;

        //- Reference to mesh
        const polyMesh& mesh_;


public:

    // Static Data Members

        //- Runtime type information
        TypeName("meshModifier");


    // Constructors

        //- Construct from components
        polyMeshModifier
        (
            const word& name,
            const polyMesh& mme
        );

        //- Disallow default bitwise copy construction
        polyMeshModifier(const polyMeshModifier&) = delete;


    //- Destructor
    virtual ~polyMeshModifier();


    // Member Functions

        //- Return name of this modifier
        const word& name() const
        {
            return name_;
        }

        //- Return reference to mesh
        const polyMesh& mesh() const;

        //- Check for topology change
        virtual bool changeTopology() const = 0;

        //- Insert the topological change instructions
        virtual void setRefinement(polyTopoChange&) const = 0;

        //- Modify motion points to comply with the topological change
        virtual void modifyMotionPoints(pointField& motionPoints) const = 0;

        //- Force recalculation of locally stored data on topological change
        virtual void topoChange(const polyTopoChangeMap&) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const polyMeshModifier&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
