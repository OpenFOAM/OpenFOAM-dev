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
    Foam::solidBodyMotionFunctions::multiMotion

Description
    Combination of SolidBodyMotionFvMesh 6DoF motion functions.

SourceFiles
    multiMotion.C

\*---------------------------------------------------------------------------*/

#ifndef multiMotion_H
#define multiMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                         Class multiMotion Declaration
\*---------------------------------------------------------------------------*/

class multiMotion
:
    public solidBodyMotionFunction
{
    // Private Data

        //- Motions to combine
        PtrList<solidBodyMotionFunction> SBMFs_;


public:

    //- Runtime type information
    TypeName("multiMotion");


    // Constructors

        //- Construct from components
        multiMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        multiMotion(const multiMotion&);

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new multiMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~multiMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
