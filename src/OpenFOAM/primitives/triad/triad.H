/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::triad

Description
    Representation of a 3D Cartesian coordinate system as a Vector of vectors.

See also
    Foam::quaternion

SourceFiles
    triadI.H
    triad.C

\*---------------------------------------------------------------------------*/

#ifndef triad_H
#define triad_H

#include "vector.H"
#include "tensor.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class triad;
Istream& operator>>(Istream&, triad&);
Ostream& operator<<(Ostream&, const triad&);

class quaternion;

/*---------------------------------------------------------------------------*\
                            Class triad Declaration
\*---------------------------------------------------------------------------*/

class triad
:
    public Vector<vector>
{

public:

    // Constructors

        //- Construct null
        inline triad();

        //- Construct from components
        inline triad(const Vector<vector>& vv);

        //- Construct from coordinate axes
        inline triad(const vector& x, const vector& y, const vector& z);

        //- Construct from a primary axis with the other two unset
        inline triad(const vector& pa);

        //- Construct from a quaternion
        triad(const quaternion& q);

        //- Construct from a tensor
        triad(const tensor& t);

        //- Construct from Istream
        inline triad(Istream&);


    // Static Data Members

        static const triad I;
        static const triad unset;


    // Member Functions

        //- Is the vector in the direction d set
        inline bool set(const direction d) const;

        //- Are all the vector set
        inline bool set() const;

        //- Return the primary direction of the vector v
        static inline direction primaryDirection(const vector& v);

        //- Return the vector orthogonal to the two provided
        static inline vector orthogonal(const vector& v1, const vector& v2);

        //- Orthogonalize this triad so that it is ortho-normal
        void orthogonalize();

        //- Normalize each set axis vector to have a unit magnitude
        void normalize();

        //- Align this triad with the given vector v
        //  by rotating the most aligned axis to be coincident with v
        void align(const vector& v);

        //- Sort the axes such that they are closest to the x, y and z axes
        triad sortxyz() const;

        //- Convert to a quaternion
        operator quaternion() const;


    // Member Operators

        inline void operator=(const Vector<vector>&);

        void operator=(const tensor& t);

        //- Add the triad t2 to this triad
        //  without normalizing or orthogonalizing
        void operator+=(const triad& t2);


    // IOstream Operators

        inline friend Istream& operator>>(Istream&, triad&);
        inline friend Ostream& operator<<(Ostream&, const triad&);
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return a quantity of the difference between two triads
scalar diff(const triad& A, const triad& B);

//- Data associated with quaternion type are contiguous
template<>
inline bool contiguous<triad>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triadI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
