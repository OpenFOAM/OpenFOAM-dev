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
    Foam::coordinateSystem

Description
    Base class for other coordinate system specifications.

    All systems are defined by an origin point and a co-ordinate rotation.

    \verbatim
    coordinateSystem
    {
        type    cartesian;
        origin  (0 0 0);
        coordinateRotation
        {
            type        cylindrical;
            e3          (0 0 1);
        }
    }
    \endverbatim

    Types of coordinateRotation:
      -# axesRotation
      -# \link STARCDCoordinateRotation STARCDRotation \endlink
      -# cylindricalCS cylindrical
      -# EulerCoordinateRotation

    Type of co-ordinates:
      -# \link cartesianCS cartesian \endlink


SourceFiles
    coordinateSystem.C
    coordinateSystemNew.C

\*---------------------------------------------------------------------------*/

#ifndef coordinateSystem_H
#define coordinateSystem_H

#include "vector.H"
#include "point.H"
#include "tensor.H"
#include "vectorField.H"
#include "pointField.H"
#include "tmp.H"
#include "coordinateRotation.H"
#include "objectRegistry.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class coordinateSystem;

bool operator!=(const coordinateSystem&, const coordinateSystem&);
Ostream& operator<<(Ostream&, const coordinateSystem&);


/*---------------------------------------------------------------------------*\
                     Class coordinateSystem Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystem
{
    // Private data

        //- Name of coordinate system
        word name_;

        //- Optional note
        string note_;

        //- Origin
        point origin_;

        //- Local-to-Global transformation tensor
        autoPtr<coordinateRotation> R_;


protected:

    // Protected Member Functions

        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual vector localToGlobal(const vector&, bool translate) const;

        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual tmp<vectorField> localToGlobal
        (
            const vectorField&,
            bool translate
        ) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual vector globalToLocal(const vector&, bool translate) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual tmp<vectorField> globalToLocal
        (
            const vectorField&,
            bool translate
        ) const;

        //- Init from dict and obr
        void init(const dictionary&);

        //- Init from dictionary
        void init(const dictionary&, const objectRegistry&);


public:

    //- Runtime type information
    TypeName("coordinateSystem");


    // Constructors

        //- Construct null. This is equivalent to an identity coordinateSystem
        coordinateSystem();

        //- Construct copy with a different name
        coordinateSystem
        (
            const word& name,
            const coordinateSystem&
        );

        //- Construct from origin and rotation
        coordinateSystem
        (
            const word& name,
            const point& origin,
            const coordinateRotation&
        );

        //- Construct from origin and 2 axes
        coordinateSystem
        (
            const word& name,
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct from dictionary with a given name
        coordinateSystem(const word& name, const dictionary&);

        //- Construct from dictionary with default name
        coordinateSystem(const dictionary&);

        //- Construct from dictionary (default name)
        //  With the ability to reference global coordinateSystems
        coordinateSystem(const objectRegistry&, const dictionary&);

        //- Construct from Istream
        //  The Istream contains a word followed by a dictionary
        coordinateSystem(Istream&);


    //- Return clone
    autoPtr<coordinateSystem> clone() const
    {
        return autoPtr<coordinateSystem>(new coordinateSystem(*this));
    }


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        coordinateSystem,
        dictionary,
        (
            const objectRegistry& obr,
            const dictionary& dict
        ),
        (obr, dict)
    );


    // Selectors

        //- Select constructed from dictionary and objectRegistry
        static autoPtr<coordinateSystem> New
        (
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Select constructed from dictionary
        static autoPtr<coordinateSystem> New
        (
            const dictionary& dict
        );

        //- Select constructed from Istream
        static autoPtr<coordinateSystem> New(Istream& is);


    //- Destructor
    virtual ~coordinateSystem();


    // Member Functions

        // Access

            //- Return name
            const word& name() const
            {
                return name_;
            }

            //- Return non-constant access to the optional note
            string& note()
            {
                return note_;
            }

            //- Return the optional note
            const string& note() const
            {
                return note_;
            }

            //- Return origin
            const point& origin() const
            {
                return origin_;
            }

            //- Return const reference to co-ordinate rotation
            const coordinateRotation& R() const
            {
                return R_();
            }

            //- Return non const reference to co-ordinate rotation
            coordinateRotation& R()
            {
                return R_();
            }

            //- Update and return the co-ordinate rotation for a list of cells
            const coordinateRotation& R
            (
                const polyMesh& mesh,
                const labelList& cells
            )
            {
                R_->updateCells(mesh, cells);
                return R_();
            }

            //- Return as dictionary of entries
            //  \param[in] ignoreType drop type (cartesian, cylindrical, etc)
            //  when generating the dictionary
            virtual dictionary dict(bool ignoreType=false) const;


        // Edit

            //- Rename
            void rename(const word& newName)
            {
                name_ = newName;
            }

            //- Edit access to origin
            point& origin()
            {
                return origin_;
            }

            //- Reset origin and rotation to an identity coordinateSystem
            //  Also resets the note
            virtual void clear();


        // Write

            //- Write
            virtual void write(Ostream&) const;

            //- Write dictionary
            void writeDict(Ostream&, bool subDict=true) const;


        // Transformations

            //- Convert from position in local coordinate system to global
            //  Cartesian position
            point globalPosition(const point& local) const
            {
                return localToGlobal(local, true);
            }

            //- Convert from position in local coordinate system to global
            //  Cartesian position
            tmp<pointField> globalPosition(const pointField& local) const
            {
                return localToGlobal(local, true);
            }

            //- Convert from vector components in local coordinate system to
            //  global Cartesian vector
            vector globalVector(const vector& local) const
            {
                return localToGlobal(local, false);
            }

            //- Convert from vector components in local coordinate system to
            //  global Cartesian vector
            tmp<vectorField> globalVector(const vectorField& local) const
            {
                return localToGlobal(local, false);
            }

            //- Convert from global Cartesian position to position in local
            //  coordinate system
            point localPosition(const point& global) const
            {
                return globalToLocal(global, true);
            }

            //- Convert from global Cartesian position to position in local
            //  coordinate system
            tmp<pointField> localPosition(const pointField& global) const
            {
                return globalToLocal(global, true);
            }

            //- Convert from global Cartesian vector to components in local
            //  coordinate system
            vector localVector(const vector& global) const
            {
                return globalToLocal(global, false);
            }

            //- Convert from global Cartesian vector to components in local
            //  coordinate system
            tmp<vectorField> localVector(const vectorField& global) const
            {
                return globalToLocal(global, false);
            }


    // Member Operators

        // friend Operators

            friend bool operator!=
            (
                const coordinateSystem&,
                const coordinateSystem&
            );


        // IOstream Operators

            friend Ostream& operator<<(Ostream&, const coordinateSystem&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
