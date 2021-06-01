/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::cubicEqn

Description
    Cubic equation of the form a*x^3 + b*x^2 + c*x + d = 0

SourceFiles
    cubicEqnI.H
    cubicEqn.C

\*---------------------------------------------------------------------------*/

#ifndef cubicEqn_H
#define cubicEqn_H

#include "Roots.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class cubicEqn Declaration
\*---------------------------------------------------------------------------*/

class cubicEqn
:
    public VectorSpace<cubicEqn, scalar, 4>
{
public:

    //- Coefficient indexing enumeration
    struct coefficient
    {
        enum {a, b, c, d};
    };


    // Constructors

        //- Construct null
        inline cubicEqn();

        //- Construct initialised to zero
        inline cubicEqn(const Foam::zero);

        //- Construct from components
        inline cubicEqn
        (
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar d
        );


    // Member Functions

        // Access

            inline scalar a() const;
            inline scalar b() const;
            inline scalar c() const;
            inline scalar d() const;

            inline scalar& a();
            inline scalar& b();
            inline scalar& c();
            inline scalar& d();

        //- Evaluate at x
        inline scalar value(const scalar x) const;

        //- Evaluate the derivative at x
        inline scalar derivative(const scalar x) const;

        //- Estimate the error of evaluation at x
        inline scalar error(const scalar x) const;

        //- Get the roots
        Roots<3> roots() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cubicEqnI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
