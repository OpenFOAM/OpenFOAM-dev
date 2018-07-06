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
    Foam::engineTime

Description
    An abstract class for the time description of the piston motion

SourceFiles
    engineTime.C

\*---------------------------------------------------------------------------*/

#ifndef engineTime_H
#define engineTime_H

#include "Time.H"
#include "IOdictionary.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class engineTime Declaration
\*---------------------------------------------------------------------------*/

class engineTime
:
    public Time
{

protected:

    const IOdictionary dict_;


public:

    //- Runtime type information
    TypeName("engineTime");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        engineTime,
        dictionary,
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const fileName& systemName,
            const fileName& constantName,
            const fileName& dictName
        ),
        (name, rootPath, caseName, systemName, constantName, dictName)
    );


    // Constructors

        //- Construct from objectRegistry arguments
        engineTime
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const fileName& systemName = "system",
            const fileName& constantName = "constant",
            const fileName& dictName = "engineGeometry"
        );


    // Selector

        static autoPtr<engineTime> New
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const fileName& systemName = "system",
            const fileName& constantName = "constant",
            const fileName& dictName = "engineGeometry"
        );


    //- Destructor
    virtual ~engineTime()
    {}


    // Member Functions
        // Conversion

            //- Calculate the piston position from the engine geometry
            //  and given timr (CA or s)
            virtual scalar pistonPosition(const scalar theta) const = 0;


        // Access

            //- Return the engine geometry dictionary
            inline const IOdictionary& engineDict() const
            {
                return dict_;
            }

            //- Return current engine time
            //  (value might be expressed in CA or s depending on the model)
            virtual scalar theta() const = 0;

            //- Return time unit
            virtual word unit() const = 0;

            //- Return engine time increment
            //  (value might be expressed in CA or s depending on the model)
            virtual scalar deltaTheta() const = 0;

            //- Return current piston position
            dimensionedScalar pistonPosition() const;

            //- Return piston displacement for current time step
            dimensionedScalar pistonDisplacement() const;

            //- Return piston speed for current time step
            dimensionedScalar pistonSpeed() const;


        // Member functions overriding the virtual functions in time

            //- Read the control dictionary and set the write controls etc.
            virtual void readDict();


        // Edit

            //- Read the controlDict and set all the parameters
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
