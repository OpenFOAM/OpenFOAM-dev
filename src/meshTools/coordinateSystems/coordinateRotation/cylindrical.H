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
    Foam::cylindrical

Description
    A local coordinate rotation.

    The rotational field can be created in two ways:
      -# Each rotational tensor is defined with two vectors (\c dir and \c axis)
         where <tt>dir = point - origin</tt> and \c axis is the rotation axis.
          Per each point an axesRotation type of rotation is created
          (cylindrical coordinates). For example:
          \verbatim
          cylindrical
          {
              type        localAxes;
              axis          (0 0 1);
          }
          \endverbatim

      -# The rotational tensor field is provided at construction.

SourceFiles
    cylindrical.C

\*---------------------------------------------------------------------------*/

#ifndef cylindrical_H
#define cylindrical_H

#include "point.H"
#include "vector.H"
#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cylindrical Declaration
\*---------------------------------------------------------------------------*/

class cylindrical
:
    public coordinateRotation
{
    // Private Data

        //- AutoPtr to transformation tensor
        autoPtr<tensorField> Rptr_;

        //- Origin of the coordinate system
        point origin_;

        //- Rotation axis
        vector axis_;


    // Private members

        //- Return the local transformation tensor
        //  corresponding to the given point
        tensor R(const vector& p) const;

        //- Initialise transformation tensor field for given points
        void init(const UList<vector>& points);


public:

    //- Runtime type information
    TypeName("cylindrical");

    // Constructors

        //- Construct from components for list of points
        cylindrical
        (
            const vector& axis,
            const point& origin,
            const UList<vector>& points
        );

        //- Construct from dictionary
        cylindrical(const dictionary&);

        //- Construct from dictionary
        cylindrical(const dictionary&, const UList<vector>& points);

        //- Return clone
        virtual autoPtr<coordinateRotation> clone() const
        {
            return autoPtr<coordinateRotation>(new cylindrical(*this));
        }


    //- Destructor
    virtual ~cylindrical()
    {}


    // Member Functions

        //- Update the rotation for a list of points
        virtual void updatePoints(const UList<vector>& points);

        //- Return local-to-global transformation tensor
        virtual const tensor& R() const
        {
            NotImplemented;
            return tensor::zero;
        }

        //- Return local Cartesian x-axis in global coordinates
        virtual const vector e1() const
        {
            NotImplemented;
            return vector::zero;
        }

        //- Return local Cartesian y-axis in global coordinates
        virtual const vector e2() const
        {
            NotImplemented;
            return vector::zero;
        }

        //- Return local Cartesian z-axis in global coordinates
        virtual const vector e3() const
        {
            return axis_;
        }

        //- Return local Cartesian z-axis in global coordinates
        virtual const vector axis() const
        {
            return axis_;
        }

        //- Return if the rotation is uniform
        virtual bool uniform() const
        {
            return false;
        }

        //- Transform vector using transformation tensor
        virtual vector transform(const vector& v) const;

        //- Transform vectorField using transformation tensor field
        virtual tmp<vectorField> transform(const vectorField& tf) const;

        //- Transform vector using transformation tensor for component
        virtual vector transform(const vector& v, const label cmptI) const;

        //- Inverse transform vector using transformation tensor
        virtual vector invTransform(const vector& v) const;

        //- Inverse transform vectorField using transformation tensor field
        virtual tmp<vectorField> invTransform(const vectorField& vf) const;

        //- Inverse transform vector using transformation tensor for component
        virtual vector invTransform(const vector& v, const label cmptI) const;

        //- Transform tensor using transformation tensorField
        virtual tensor transformTensor(const tensor& t) const;

        //- Transform tensor field using transformation tensorField
        virtual tmp<tensorField> transformTensor(const tensorField& tf) const;

        //- Transform vector using transformation tensor and return
        //  symmetrical tensor
        virtual symmTensor transformVector(const vector& v) const;

        //- Transform vectorField using transformation tensorField and return
        //  symmetrical tensorField
        virtual tmp<symmTensorField> transformVector
        (
            const vectorField& vf
        ) const;


    // Write

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
