/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::spatialTransform

Description
    Compact representation of the Plücker spatial transformation tensor
    in terms of the rotation tensor \c E and translation vector \c r .

    See Chapter 2 and Appendix A in reference:
    \verbatim
        Featherstone, R. (2008).
        Rigid body dynamics algorithms.
        Springer.
    \endverbatim

SourceFiles
    spatialTransformI.H

\*---------------------------------------------------------------------------*/

#ifndef spatialTransform_H
#define spatialTransform_H

#include "tensor.H"
#include "spatialVector.H"
#include "spatialTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class spatialTransform;
Istream& operator>>(Istream&, spatialTransform&);
Ostream& operator<<(Ostream&, const spatialTransform&);


/*---------------------------------------------------------------------------*\
                      Class spatialTransform Declaration
\*---------------------------------------------------------------------------*/

class spatialTransform
{
    // Private Data

        //- Rotation tensor
        tensor E_;

        //- Translation vector
        vector r_;


    // Private Member Functions

        //- Return E . *r
        inline tensor Erx() const;


public:

        //- Wrapper-class to provide transpose functions and operators
        class transpose
        {
            const spatialTransform& X_;

        public:

            //- Construct from a spatialTransform
            inline transpose(const spatialTransform& X);

            //- Return the transpose transformation tensor ^A{X^*}_B
            //  X^T
            inline operator spatialTensor() const;

            //- Transpose transform dual f: ^A{X^*}_B & f
            //  X^T . f = (E^T . fl + r ^ E^T . fw, E^T . fl)
            inline spatialVector operator&(const spatialVector& f) const;
        };


        //- Wrapper-class to provide dual functions and operators
        class dual
        {
            const spatialTransform& X_;

        public:

            //- Construct from a spatialTransform
            inline dual(const spatialTransform& X);

            //- Return dual transformation tensor ^B{X^*}_A
            inline operator spatialTensor() const;

            //- Transform dual f: ^B{X^*}_A & f
            //  X^* . f = (E . fw - r ^ fl, E . fl)
            inline spatialVector operator&(const spatialVector& f) const;
        };


    // Constructors

        //- Construct null
        inline spatialTransform();

        //- Construct from components
        inline spatialTransform(const tensor& E, const vector& r);

        //- Construct from Istream
        inline spatialTransform(Istream&);


    // Member Functions

        //- Return the rotation tensor
        inline const tensor& E() const;

        //- Return non-const access to the rotation tensor
        inline tensor& E();

        //- Return the translation vector
        inline const vector& r() const;

        //- Return non-const access to the translation vector
        inline vector& r();

        //- Return the transpose transformation tensor ^A{X^*}_B
        //  X^T
        inline transpose T() const;

        //- Return the inverse transformation tensor: X^-1
        //  X^-1 = (E^T, −E.r)
        inline spatialTransform inv() const;


    // Member Operators

        //- Return the dual transformation tensor ^B{X^*}_A
        inline dual operator*() const;

        //- Return transformation tensor ^BX_A
        //  X
        inline operator spatialTensor() const;

        //- Inner-product multiply with a transformation tensor
        inline void operator&=(const spatialTransform& X);

        //- Return the inner-product of two transformation tensors
        inline spatialTransform operator&(const spatialTransform& X) const;

        //- Transform v: ^BX_A . v
        //  X.v = (E . vw, E . (vl - r^vw))
        inline spatialVector operator&(const spatialVector& v) const;

        //- Transform position p
        //  X:p = E . (pl - r)
        inline vector transformPoint(const vector& p) const;

        //- Transform position p
        //  X:p = (E . pw, E . (vl - r))
        inline spatialVector operator&&(const spatialVector& v) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, spatialTransform&);
        friend Ostream& operator<<(Ostream&, const spatialTransform&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "spatialTransformI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
