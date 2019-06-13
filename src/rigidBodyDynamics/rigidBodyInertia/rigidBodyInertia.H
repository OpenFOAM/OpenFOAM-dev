/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::RBD::rigidBodyInertia

Description
    This class represents the linear and angular inertia of a rigid body
    by the mass, centre of mass and moment of inertia tensor about the
    centre of mass.

    Reference:
    \verbatim
        Featherstone, R. (2008).
        Rigid body dynamics algorithms.
        Springer.
    \endverbatim

SourceFiles
    rigidBodyInertiaI.H
    rigidBodyInertia.C
    rigidBodyInertiaIO.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_rigidBodyInertia_H
#define RBD_rigidBodyInertia_H

#include "vector.H"
#include "symmTensor.H"
#include "spatialVector.H"
#include "spatialTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

namespace RBD
{

// Forward declaration of friend functions and operators
class rigidBodyInertia;
Istream& operator>>(Istream&, rigidBodyInertia&);
Ostream& operator<<(Ostream&, const rigidBodyInertia&);


/*---------------------------------------------------------------------------*\
                         Class rigidBodyInertia Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyInertia
{
    // Private Data

        //- Mass of the rigid-body
        scalar m_;

        //- Centre of mass of the rigid-body
        vector c_;

        //- Inertia tensor about the centre of mass
        symmTensor Ic_;


public:

    // Static member functions

        //- Return the difference between the inertia tensor of the rigid-body
        //  about the origin - about the centre of mass
        //  for the given mass and centre of mass
        inline static symmTensor Ioc(const scalar m, const vector& c);


    // Constructors

        //- Null constructor, initializes to zero
        inline rigidBodyInertia();

        //- Construct from mass, centre of mass and moment of inertia tensor
        //  about the centre of mass
        inline rigidBodyInertia
        (
            const scalar m,
            const vector& c,
            const symmTensor& Ic
        );

        //- Construct from dictionary
        inline rigidBodyInertia(const dictionary& dict);

        //- Construct from the components of a spatial tensor
        inline explicit rigidBodyInertia(const spatialTensor& st);

        //- Construct from Istream
        inline explicit rigidBodyInertia(Istream& is);


    // Member Functions

        //- Return the mass of the rigid-body
        inline scalar m() const;

        //- Return the centre of mass of the rigid-body
        inline const vector& c() const;

        //- Return the inertia tensor of the rigid-body about the centre of mass
        inline const symmTensor& Ic() const;

        //- Return the difference between the inertia tensor of the rigid-body
        //  about the origin - about the centre of mass
        inline symmTensor Ioc() const;

        //- Return the difference between the inertia tensor of the rigid-body
        //  about the a new centre of mass - about the current centre of mass
        inline symmTensor Icc(const vector& c) const;

        //- Return the inertia tensor of the rigid-body about the origin
        inline symmTensor Io() const;

        //- Return the kinetic energy of the body with the given velocity
        inline scalar kineticEnergy(const spatialVector& v);


    // Member Operators

        //- Conversion to spatial tensor
        inline operator spatialTensor() const;

        inline void operator+=(const rigidBodyInertia&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, rigidBodyInertia&);
        friend Ostream& operator<<(Ostream&, const rigidBodyInertia&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rigidBodyInertiaI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
