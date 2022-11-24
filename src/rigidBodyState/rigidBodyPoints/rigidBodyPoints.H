/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::functionObjects::rigidBodyPoints

Description
    Writes the position, linear and angular velocities and accelerations of a
    list of points on a body specified in the body local coordinate system.

Usage
    \table
        Property     | Description                  | Required | Default value
        type         | type name: rigidBodyPoints   | yes      |
        angleUnits   | degrees or radians           | no       | radians
        body         | name of the body             | yes      |
        points       | list of points on the body   | yes      |
    \endtable

    Example of function object specification:
    \verbatim
    rigidBodyPoints
    {
        type           rigidBodyPoints;
        libs           ("librigidBodyState.so");

        angleUnits     degrees;

        body           floatingObject;

        points
        (
            point1     (0 0 0)
            point2     (0.1 0.1 0.25)
        );
    }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    rigidBodyPoints.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyPoints_H
#define rigidBodyPoints_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "rigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class rigidBodyPoints Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyPoints
:
    public fvMeshFunctionObject,
    public logFiles
{
    // Private Data

        //- Angle units, radians (default) or degrees
        word angleUnits_;

        //- Name of the body
        word body_;

        //- List of points on the body
        List<point> points_;

        //- Names of the body point files
        wordList names_;


    // Private Member Functions

        const RBD::rigidBodyMotion& motion() const;


protected:

    // Protected Member Functions

        //- overloaded writeFileHeader from writeFile
        virtual void writeFileHeader(const label i = 0);


public:

    //- Runtime type information
    TypeName("rigidBodyPoints");


    // Constructors

        //- Construct from Time and dictionary
        rigidBodyPoints
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rigidBodyPoints(const rigidBodyPoints&) = delete;


    //- Destructor
    virtual ~rigidBodyPoints();


    // Member Functions

        //- Read the rigidBodyPoints data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the rigidBodyPoints
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rigidBodyPoints&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
