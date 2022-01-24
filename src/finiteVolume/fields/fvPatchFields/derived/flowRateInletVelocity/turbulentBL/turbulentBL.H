/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::Function1s::turbulentBL

Description
    Turbulent boundary-layer function pow(yHat, exponent)

    where yHat is the normalised distance to the wall and the optional exponent
    defaults to 1/7.

    Usage:
    \verbatim
        <name> turbulentBL;
    \endverbatim
    or
    \verbatim
        <name>
        {
            type      turbulentBL;
            exponent  0.143;
        }
    \endverbatim

See also
    Foam::Function1
    Foam::Function1s::laminarBL

SourceFiles
    turbulentBL.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentBL_H
#define turbulentBL_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class turbulentBL Declaration
\*---------------------------------------------------------------------------*/

class turbulentBL
:
    public FieldFunction1<scalar, turbulentBL>
{
    // Private data

        //- Exponent of the pow function
        scalar exponent_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& dict);


public:

    // Runtime type information
    TypeName("turbulentBL");


    // Constructors

        //- Construct from name and dictionary
        turbulentBL
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    virtual ~turbulentBL();


    // Member Functions

        //- Return value for time t
        virtual inline scalar value(const scalar t) const;

        //- Return the integral between times t1 and t2
        virtual inline scalar integral
        (
            const scalar t1,
            const scalar t2
        ) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const turbulentBL&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentBLI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
