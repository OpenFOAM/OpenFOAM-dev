/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::sampledSets

Description
    Set of sets to sample.
    Call sampledSets.write() to sample&write files.

SourceFiles
    sampledSets.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSets_H
#define sampledSets_H

#include "fvMeshFunctionObject.H"
#include "sampledSet.H"
#include "volFieldsFwd.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "coordSet.H"
#include "setWriter.H"
#include "wordReList.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class objectRegistry;
class dictionary;
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class sampledSets Declaration
\*---------------------------------------------------------------------------*/

class sampledSets
:
    public fvMeshFunctionObject,
    public PtrList<sampledSet>
{
    // Private Data

        //- Keep the dictionary to recreate sets for moving mesh cases
        dictionary dict_;

        //- Output path
        fileName outputPath_;

        //- Mesh search engine
        meshSearch searchEngine_;


        // Read from dictionary

            //- Names of fields to sample
            wordList fields_;

            //- Interpolation scheme to use
            word interpolationScheme_;


        // Merging structures

            //- Global sets gathered on the master
            PtrList<coordSet> masterSets_;

            //- Mapping from local to global sets
            labelListList masterSetOrders_;


        // Calculated

            //- Set formatter
            autoPtr<setWriter> formatter_;


    // Private Member Functions

        //- Combine sets on all processors. Sort by curveDist and produce
        //  index list. Valid result only on master processor.
        void combineSampledSets();

        //- Correct for mesh changes
        void correct();

        //- Sample all fields of a type on a given set
        template<class Type>
        PtrList<Field<Type>> sampleLocalType
        (
            const label seti,
            const wordList& fieldNames,
            HashPtrTable<interpolation<Type>>& interpolations
        );

        //- Sample all fields of a type on a given set and combine on the
        //  master for writing
        template<class Type>
        PtrList<Field<Type>> sampleType
        (
            const label seti,
            const wordList& fieldNames,
            HashPtrTable<interpolation<Type>>& interpolations
        );


public:

    //- Runtime type information
    TypeName("sets");


    // Constructors

        //- Construct from Time and dictionary
        sampledSets
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        sampledSets(const sampledSets&) = delete;


    //- Destructor
    virtual ~sampledSets();


    // Member Functions

        //- Read the sampledSets
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sampledSets&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSetsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
