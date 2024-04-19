/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::EulerCoordinateRotation

Description
    A coordinateRotation defined in the z-x-y Euler convention.

    The 3 rotations are defined in the Euler convention
    (around Z, around X' and around Z').
    For reference and illustration, see
    http://mathworld.wolfram.com/EulerAngles.html
    Note, however, that it is the reverse transformation
    (local->global) that is defined here.

    - the rotation angles are in degrees, unless otherwise explicitly specified:

    \verbatim
    coordinateRotation
    {
        type        EulerRotation;
        degrees     false;
        rotation    (0 0 3.141592654);
    }
    \endverbatim

SourceFiles
    EulerCoordinateRotation.C

\*---------------------------------------------------------------------------*/

#ifndef EulerCoordinateRotation_H
#define EulerCoordinateRotation_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class EulerCoordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class EulerCoordinateRotation
:
    public coordinateRotation
{

    // Private Member Data

        //- Local-to-global transformation tensor
        tensor R_;

        //- Global-to-Local transformation tensor
        tensor Rtr_;


    // Private Member Functions

        //- Calculate transformation tensor
        void calcTransform
        (
            const scalar phiAngle,
            const scalar thetaAngle,
            const scalar psiAngle
        );


public:

    //- Runtime type information
    TypeName("EulerRotation");


    // Constructors

        //- Construct from dictionary
        EulerCoordinateRotation(const dictionary&);

        //- Construct from dictionary and list of points
        EulerCoordinateRotation(const dictionary&, const UList<vector>& points);

        //- Construct and return a clone
        virtual autoPtr<coordinateRotation> clone() const
        {
            return autoPtr<coordinateRotation>
            (
                new EulerCoordinateRotation(*this)
            );
        }


    // Member Functions

        //- Update the rotation for a list of points
        virtual void updatePoints(const UList<vector>& points)
        {}

        //- Return local-to-global transformation tensor
        virtual const tensor& R() const
        {
            return R_;
        }

        //- Return local Cartesian x-axis in global coordinates
        virtual const vector e1() const
        {
            return Rtr_.x();
        }

        //- Return local Cartesian y-axis in global coordinates
        virtual const vector e2() const
        {
            return Rtr_.y();
        }

        //- Return local Cartesian z-axis in global coordinates
        virtual const vector e3() const
        {
            return Rtr_.z();
        }

        //- Transform vectorField using transformation tensor field
        virtual tmp<vectorField> transform(const vectorField& vf) const;

        //- Transform vector using transformation tensor
        virtual vector transform(const vector& v) const;

        //- Inverse transform vectorField using transformation tensor field
        virtual tmp<vectorField> invTransform(const vectorField& vf) const;

        //- Inverse transform vector using transformation tensor
        virtual vector invTransform(const vector& v) const;

        //- Transform tensor field using transformation tensorField
        virtual tmp<tensorField> transform(const tensorField& tf) const;

        //- Transform tensor using transformation tensorField
        virtual tensor transform(const vector& p, const tensor& t) const;

        //- Transform diagTensor masquerading as a vector using transformation
        //  tensor and return symmTensor
        virtual symmTensor transformDiagTensor
        (
            const vector& p,
            const vector& v
        ) const;

        //- Transform diagTensorField masquerading as a vectorField
        //  using transformation tensorField and return symmTensorField
        virtual tmp<symmTensorField> transformDiagTensor
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
