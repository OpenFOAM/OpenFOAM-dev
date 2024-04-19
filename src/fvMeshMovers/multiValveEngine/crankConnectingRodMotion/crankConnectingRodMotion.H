/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::Function1s::crankConnectingRodMotion

Description
    Crank and connecting-rod motion function for piston engines etc.

    \c value returns the current piston position as a function of
    the crank-angle in degrees.

    Usage:
    \verbatim
        <name>
        {
            type            crankConnectingRodMotion;

            conRodLength    0.2;
            stroke          0.12;
        }
    \endverbatim

See also
    Foam::Function1

SourceFiles
    crankConnectingRodMotionI.H
    crankConnectingRodMotion.C

\*---------------------------------------------------------------------------*/

#ifndef crankConnectingRodMotion_H
#define crankConnectingRodMotion_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                   Class crankConnectingRodMotion Declaration
\*---------------------------------------------------------------------------*/

class crankConnectingRodMotion
:
    public Function1<scalar>
{
    // Private data

        //- Connecting-rod length
        const scalar conRodLength_;

        //- Stroke
        const scalar stroke_;


public:

    // Runtime type information
    TypeName("crankConnectingRodMotion");


    // Constructors

        //- Construct from name and dictionary
        crankConnectingRodMotion
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const
        {
            return tmp<Function1<scalar>>
            (
                new crankConnectingRodMotion(*this)
            );
        }


    //- Destructor
    virtual ~crankConnectingRodMotion();


    // Member Functions

        //- Return position for crank-angle theta in deg
        virtual inline scalar value(const scalar theta) const;

        //- Not implemented
        virtual inline tmp<Field<scalar>> value(const scalarField&) const;

        //- Not implemented
        virtual inline scalar integral
        (
            const scalar,
            const scalar
        ) const;

        //- Not implemented
        virtual inline tmp<Field<scalar>> integral
        (
            const scalarField&,
            const scalarField&
        ) const;


        //- Write data to dictionary stream
        virtual void write(Ostream& os, const unitConversions& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const crankConnectingRodMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "crankConnectingRodMotionI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
