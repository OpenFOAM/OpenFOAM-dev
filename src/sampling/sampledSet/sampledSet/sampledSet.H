/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::sampledSet

Description
    Holds list of sampling points which is filled at construction time.
    Various implementations of this base class to e.g. get sampling points
    at uniform distance along a line (lineUniformSet) or directly specified
    (pointsSet)

    Each 'sampledSet' has a name and a specifier of how the axis should be
    write (x/y/z component or all 3 components)

SourceFiles
    sampledSet.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSet_H
#define sampledSet_H

#include "coordSet.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class meshSearch;

/*---------------------------------------------------------------------------*\
                         Class sampledSet Declaration
\*---------------------------------------------------------------------------*/

class sampledSet
:
    public coordSet
{
    // Private Data

        //- Name of the set
        const word name_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Reference to mesh searching class
        const meshSearch& searchEngine_;


protected:

    // Protected data

        //- Cell numbers
        labelList cells_;

        //- Face numbers (-1 if not known)
        labelList faces_;


    // Protected Member Functions

        //- Sets sample data
        void setSamples
        (
            const List<point>& samplingPositions,
            const labelList& samplingSegments,
            const labelList& samplingCells,
            const labelList& samplingFaces
        );

        //- Sets sample data
        void setSamples
        (
            const List<point>& samplingPositions,
            const List<scalar>& samplingDistances,
            const labelList& samplingSegments,
            const labelList& samplingCells,
            const labelList& samplingFaces
        );


public:

    //- Runtime type information
    TypeName("sampledSet");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sampledSet,
            word,
            (
                const word& name,
                const polyMesh& mesh,
                const meshSearch& searchEngine,
                const dictionary& dict
            ),
            (name, mesh, searchEngine, dict)
        );


        //- Class used for the read-construction of
        //  PtrLists of sampledSet
        class iNew
        {
            const polyMesh& mesh_;
            const meshSearch& searchEngine_;

        public:

            iNew(const polyMesh& mesh, const meshSearch& searchEngine)
            :
                mesh_(mesh),
                searchEngine_(searchEngine)
            {}

            autoPtr<sampledSet> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);
                return sampledSet::New(name, mesh_, searchEngine_, dict);
            }
        };


    // Constructors

        //- Construct from components
        sampledSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const word& axis
        );

        //- Construct from dictionary
        sampledSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );

        //- Clone
        autoPtr<sampledSet> clone() const
        {
            NotImplemented;
            return autoPtr<sampledSet>(nullptr);
        }


    // Selectors

        //- Return a reference to the selected sampledSet
        static autoPtr<sampledSet> New
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledSet();


    // Member Functions

        //- Access the name
        const word& name() const
        {
            return name_;
        }

        //- Access the mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Access the search engine
        const meshSearch& searchEngine() const
        {
            return searchEngine_;
        }

        //- Access the positions
        const pointField& positions() const
        {
            return positions_();
        }

        //- Access the cells
        const labelList& cells() const
        {
            return cells_;
        }

        //- Access the faces
        const labelList& faces() const
        {
            return faces_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
