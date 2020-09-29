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

#include "coordinateRotation.H"
#include "objectRegistry.H"


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
    // Private Data

        //- Name of coordinate system
        word name_;

        //- Origin
        point origin_;

        //- Local-to-Global transformation tensor
        mutable autoPtr<coordinateRotation> R_;


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


public:

    //- Runtime type information
    TypeName("coordinateSystem");


    // Constructors

        //- Construct from origin
        coordinateSystem
        (
            const word& name,
            const point& origin
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

        //- Construct and return a clone
        virtual autoPtr<coordinateSystem> clone() const
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
            const word& name,
            const dictionary& dict
        ),
        (name, dict)
    );


    // Selectors

        //- Select constructed from dictionary and objectRegistry
        static autoPtr<coordinateSystem> New
        (
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Select constructed from name and dictionary
        static autoPtr<coordinateSystem> New
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    virtual ~coordinateSystem();


    // Member Functions

        // Access

            //- Return name
            const word& name() const
            {
                return name_;
            }

            //- Return keyword
            const word& keyword() const
            {
                return name_;
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

            //- Update and return the co-ordinate rotation for a list of points
            const coordinateRotation& R(const UList<vector>& points) const
            {
                R_->updatePoints(points);
                return R_();
            }


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


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const coordinateSystem&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
