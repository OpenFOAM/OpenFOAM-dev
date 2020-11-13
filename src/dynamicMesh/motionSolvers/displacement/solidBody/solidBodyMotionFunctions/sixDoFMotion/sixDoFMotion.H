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
    Foam::solidBodyMotionFunctions::sixDoFMotion

Description
    Tabulated 6DoF motion function.

    Obtained by interpolating tabulated data for surge (x-translation),
    sway (y-translation), heave (z-translation), roll (rotation about x),
    pitch (rotation about y) and yaw (rotation about z).

SourceFiles
    sixDoFMotion.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFMotion_H
#define sixDoFMotion_H

#include "solidBodyMotionFunction.H"
#include "Function1.H"
#include "Vector2D.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                     Class sixDoFMotion Declaration
\*---------------------------------------------------------------------------*/

class sixDoFMotion
:
    public solidBodyMotionFunction
{
public:

        //- Type used to read in the translation and rotation "vectors"
        typedef Vector2D<vector> translationRotationVectors;


private:

    // Private Data

        //- Center of gravity read from dictionary
        vector CofG_;

        //- Translation and rotation function
        autoPtr<Function1<translationRotationVectors>> translationRotation_;


public:

    //- Runtime type information
    TypeName("sixDoFMotion");


    // Constructors

        //- Construct from components
        sixDoFMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        sixDoFMotion(const sixDoFMotion&);


        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new sixDoFMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~sixDoFMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sixDoFMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
