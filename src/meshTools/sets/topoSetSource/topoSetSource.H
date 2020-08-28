/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::topoSetSource

Description
    Base class of a source for a topoSet.

    Implementer has to modify the given set (see applyToSet) according to
    its function and the setAction (one of add/delete/new)

SourceFiles
    topoSetSource.C

\*---------------------------------------------------------------------------*/

#ifndef topoSetSource_H
#define topoSetSource_H

#include "pointField.H"
#include "word.H"
#include "labelList.H"
#include "faceList.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "NamedEnum.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class topoSet;

/*---------------------------------------------------------------------------*\
                        Class topoSetSource Declaration
\*---------------------------------------------------------------------------*/

class topoSetSource
{
public:

    // Public data types

        //- Enumeration defining the types of sources
        enum sourceType
        {
            CELLSETSOURCE,
            FACESETSOURCE,
            POINTSETSOURCE,

            CELLZONESOURCE,
            FACEZONESOURCE,
            POINTZONESOURCE
        };

        //- Enumeration defining the valid actions
        enum setAction
        {
            CLEAR,
            NEW,
            INVERT,
            ADD,
            DELETE,
            SUBSET,
            LIST,
            REMOVE
        };

protected:

        //- A table of usage strings
        static HashTable<string>* usageTablePtr_;

        //- Class with constructor to add usage string to table
        class addToUsageTable
        {
        public:

            addToUsageTable(const word& name, const string& msg)
            {
                if (!usageTablePtr_)
                {
                    usageTablePtr_ = new HashTable<string>();
                }
                usageTablePtr_->insert(name, msg);
            }

            ~addToUsageTable()
            {
                if (usageTablePtr_)
                {
                    delete usageTablePtr_;
                    usageTablePtr_ = nullptr;
                }
            }
        };


    // Protected data

        const polyMesh& mesh_;

        //- Add (if bool) celli to set or delete celli from set.
        void addOrDelete(topoSet& set, const label celli, const bool) const;


private:

        static const NamedEnum<setAction, 8> actionNames_;

        static const string illegalSource_;


public:

    //- Runtime type information
    TypeName("topoSetSource");


    // Static Functions

        //- Convert string to action
        static setAction toAction(const word& actionName)
        {
            return actionNames_[actionName];
        }

        //- Check state of stream.
        static Istream& checkIs(Istream& is);

    // Declare run-time constructor selection table

        // For the dictionary constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSetSource,
            word,
            (
                const polyMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

        // For the Istream constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSetSource,
            istream,
            (
                const polyMesh& mesh,
                Istream& is
            ),
            (mesh, is)
        );


        //- Class used for the read-construction of
        //  PtrLists of topoSetSource
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<topoSetSource> operator()(Istream& is) const
            {
                word topoSetSourceType(is);
                dictionary dict(is);
                return topoSetSource::New(topoSetSourceType, mesh_, dict);
            }
        };


        static const string& usage(const word& name)
        {
            if (!usageTablePtr_)
            {
                usageTablePtr_ = new HashTable<string>();
            }

            const HashTable<string>& usageTable = *usageTablePtr_;

            if (usageTable.found(name))
            {
                return usageTable[name];
            }
            else
            {
                return illegalSource_;
            }
        }


    // Constructors

        //- Construct from components
        topoSetSource(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        topoSetSource(const topoSetSource&) = delete;

        //- Clone
        autoPtr<topoSetSource> clone() const
        {
            NotImplemented;
            return autoPtr<topoSetSource>(nullptr);
        }


    // Selectors

        //- Return a reference to the selected topoSetSource
        static autoPtr<topoSetSource> New
        (
            const word& topoSetSourceType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Return a reference to the selected topoSetSource
        static autoPtr<topoSetSource> New
        (
            const word& topoSetSourceType,
            const polyMesh& mesh,
            Istream& is
        );


    //- Destructor
    virtual ~topoSetSource();


    // Member Functions

        const polyMesh& mesh() const
        {
            return mesh_;
        }


    // Member Functions

        virtual sourceType setType() const = 0;

        virtual void applyToSet(const setAction action, topoSet&) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const topoSetSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
