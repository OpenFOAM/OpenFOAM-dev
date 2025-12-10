/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::solidBodyMotionFunctions::multiValveEnginePistonMotion

Description
    Motion function for solid regions representing a piston. Looks up the
    multiValveEngine mover from the associated fluid region and uses the motion
    specified for the piston.

    Example specification, in constant/<solidRegion>/dynamicMeshDict:
    \verbatim
    mover
    {
        type            motionSolver;
        libs            ("libfvMotionSolvers.so");
        motionSolver    solidBody;
        cellZone        all;
        solidBodyMotionFunction multiValveEnginePistonMotion;
        fluidRegion     fluid;
    }
    \endverbatim

SourceFiles
    multiValveEnginePistonMotion.C

\*---------------------------------------------------------------------------*/

#ifndef multiValveEnginePistonMotion_H
#define multiValveEnginePistonMotion_H

#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                 Class multiValveEnginePistonMotion Declaration
\*---------------------------------------------------------------------------*/

class multiValveEnginePistonMotion
:
    public solidBodyMotionFunction
{
private:

    // Private Data

        //- Name of the region in which the multiValveEngine mover is used
        word fluidRegionName_;


public:

    //- Runtime type information
    TypeName("multiValveEnginePistonMotion");


    // Constructors

        //- Construct from components
        multiValveEnginePistonMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        multiValveEnginePistonMotion
        (
            const multiValveEnginePistonMotion&
        ) = delete;

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new multiValveEnginePistonMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~multiValveEnginePistonMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiValveEnginePistonMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
