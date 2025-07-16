/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;

/*---------------------------------------------------------------------------*\
                         Class sampledSet Declaration
\*---------------------------------------------------------------------------*/

class sampledSet
{
    // Private Data

        //- Name of the set
        const word name_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Coordinate set
        mutable autoPtr<coordSet> coordsPtr_;

        //- Cell numbers
        mutable autoPtr<labelList> cellsPtr_;

        //- Face numbers (-1 if not known)
        mutable autoPtr<labelList> facesPtr_;

        //- Axis
        const coordSet::axisType axis_;


    // Private Member Functions

        //- Calc the samples
        virtual bool calcSamples
        (
            DynamicList<point>& samplingPositions,
            DynamicList<scalar>& samplingDistances,
            DynamicList<label>& samplingSegments,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces
        ) const = 0;

        //- Generate the samples
        void setSamples() const;


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
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


        //- Class used for the read-construction of PtrLists of sampledSet
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<sampledSet> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);
                return sampledSet::New(name, mesh_, dict);
            }
        };


    // Constructors

        //- Construct from components
        sampledSet
        (
            const word& name,
            const polyMesh& mesh,
            const word& axis
        );

        //- Construct from dictionary
        sampledSet
        (
            const word& name,
            const polyMesh& mesh,
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
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledSet();


    // Member Functions

        //- Access the name
        inline const word& name() const;

        //- Access the mesh
        inline const polyMesh& mesh() const;

        //- Access the coordinate set
        inline const coordSet& coords() const;

        //- Return the size
        inline label size() const;

        //- Return the size
        inline const labelList& segments() const;

        //- Access the positions
        inline const pointField& positions() const;

        //- Access the cells
        inline const labelList& cells() const;

        //- Access the faces
        inline const labelList& faces() const;

        //- Update for mesh point-motion
        virtual void movePoints();

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Cast to the coordinate set
        inline operator const coordSet&() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sampledSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
