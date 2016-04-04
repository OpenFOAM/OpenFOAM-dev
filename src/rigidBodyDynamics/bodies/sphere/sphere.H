/*---------------------------------------------------------------------------*\3
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::sphere

Description

SourceFiles
    sphereI.H
    sphere.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_sphere_H
#define RBD_sphere_H

#include "rigidBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of Foam classes
class Istream;
class Ostream;

namespace RBD
{

// Forward declaration of friend functions and operators
class sphere;
Ostream& operator<<(Ostream&, const sphere&);


/*---------------------------------------------------------------------------*\
                         Class sphere Declaration
\*---------------------------------------------------------------------------*/

class sphere
:
    public rigidBody
{

public:

    // Constructors

        //- Construct from name, mass and radius
        inline sphere(const word& name, const scalar m, const scalar r);


    //- Destructor
    virtual ~sphere();


    // Member Functions

    // IOstream Operators

        //friend Ostream& operator<<(Ostream&, const sphere&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sphereI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
