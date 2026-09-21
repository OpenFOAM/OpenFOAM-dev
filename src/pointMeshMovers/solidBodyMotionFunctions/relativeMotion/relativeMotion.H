/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::solidBodyMotionFunctions::relativeMotion

Description
    Solid-body motion relative to a reference motion specified in the same list.

Usage
    Example of multiSolidBody motion specified in dynamicMeshDict:
    \verbatim
    mover
    {
        type        multiSolidBody;
        libs        ("libpointMeshMovers.so");

        outer
        {
            type        rotatingMotion;

            origin      (5 5 0);
            axis        (0 0 1);
            omega       0.1;
        }

        inner
        {
            type        relativeMotion;

            relativeMotion  rotatingMotion;
            relativeTo  outer;

            origin      (4 5 0);
            axis        (0 0 1);
            omega       -0.2;
        }
    }
    \endverbatim
    The inner region is specified to rotate within and relative to the
    rotation of the outer region.

See also
    Foam::solidBodyMotionFunctions::multiMotion

SourceFiles
    relativeMotion.C

\*---------------------------------------------------------------------------*/

#ifndef relativeMotion_H
#define relativeMotion_H

#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                       Class relativeMotion Declaration
\*---------------------------------------------------------------------------*/

class relativeMotion
:
    public solidBodyMotionFunction
{
    // Private Data

        //- List of solidBodyMotionFunctions
        const PtrListDictionary<solidBodyMotionFunction>& SBMFs_;

        //- Name of reference motion function
        const word referenceName_;

        //- Reference motion
        autoPtr<solidBodyMotionFunction> SBMF_;


public:

    //- Runtime type information
    TypeName("relativeMotion");


    // Constructors

        //- Construct from solidBodyMotionFunction list
        //  and coefficients dictionary
        relativeMotion
        (
            const word& name,
            const PtrListDictionary<solidBodyMotionFunction>& SBMFs,
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        relativeMotion(const relativeMotion&) = delete;

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            NotImplemented;
            return autoPtr<solidBodyMotionFunction>(nullptr);
        }


    //- Destructor
    virtual ~relativeMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const relativeMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
