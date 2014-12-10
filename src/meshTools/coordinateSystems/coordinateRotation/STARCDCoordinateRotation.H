/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::STARCDCoordinateRotation

Description
    A coordinateRotation defined by the STAR-CD convention.

    The 3 rotations are defined in the STAR-CD convention
    (around Z, around X' and around Y'').
    The order of the parameter arguments matches this rotation order.

    - the rotation angles are in degrees, unless otherwise explictly specified:

    \verbatim
        coordinateRotation
        {
            type        STARCDRotation;
            degrees     false;
            rotation    (0 0 3.141592654);
        }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef STARCDCoordinateRotation_H
#define STARCDCoordinateRotation_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class STARCDCoordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class STARCDCoordinateRotation
:
    public coordinateRotation
{

    // Private Member Data

        //- Local-to-Global transformation tensor
        tensor R_;

        //- Global-to-Local transformation tensor
        tensor Rtr_;


    // Private Member Functions

        //- Calculate transformation tensor
        void calcTransform
        (
            const scalar rotZ,
            const scalar rotX,
            const scalar rotY,
            const bool inDegrees=true
        );


public:

    //- Runtime type information
    TypeName("STARCDRotation");


    // Constructors

        //- Construct null
        STARCDCoordinateRotation();

        //- Construct from rotation vector
        STARCDCoordinateRotation
        (
            const vector& rotZrotXrotY,
            const bool inDegrees=true
        );

        //- Construct from components of rotation vector
        STARCDCoordinateRotation
        (
            const scalar rotZ,
            const scalar rotX,
            const scalar rotY,
            const bool inDegrees=true
        );

        //- Construct from dictionary
        STARCDCoordinateRotation(const dictionary&);

        //- Construct from dictionary and mesh
        STARCDCoordinateRotation(const dictionary&, const objectRegistry&);

         // Member Functions

        //- Reset rotation to an identity rotation
        virtual void clear()
        {
            R_ = sphericalTensor::I;
            Rtr_ = sphericalTensor::I;
        }

        //- Update the rotation for a list of cells
        virtual void updateCells(const polyMesh&, const labelList&)
        {
            // do nothing
        }

        //- Return local-to-global transformation tensor
        virtual const tensor& R() const
        {
            return R_;
        }

        //- Return global-to-local transformation tensor
        virtual const tensor& Rtr() const
        {
            return Rtr_;
        };

        //- Return local Cartesian x-axis
        virtual const vector e1() const
        {
            return R_.x();
        }

        //- Return local Cartesian y-axis
        virtual const vector e2() const
        {
            return R_.y();
        }

        //- Return local Cartesian z-axis
        virtual const vector e3() const
        {
            return R_.z();
        }

        //- Return transformation tensor field
        virtual const tensorField& Tr() const;

        //- Transform vectorField using transformation tensor field
        virtual tmp<vectorField> transform(const vectorField& st) const;

        //- Transform vector using transformation tensor
        virtual vector transform(const vector& st) const;

        //- Inverse transform vectorField using transformation tensor field
        virtual tmp<vectorField> invTransform(const vectorField& st) const;

        //- Inverse transform vector using transformation tensor
        virtual vector invTransform(const vector& st) const;

        //- Transform tensor field using transformation tensorField
        virtual tmp<tensorField> transformTensor(const tensorField& st) const;

        //- Transform tensor using transformation tensorField
        virtual tensor transformTensor(const tensor& st) const;

        //- Transform tensor sub-field using transformation tensorField
        virtual tmp<tensorField> transformTensor
        (
            const tensorField& st,
            const labelList& cellMap
        ) const;

        //- Transform vectorField using transformation tensorField and return
        // symmetrical tensorField
        virtual tmp<symmTensorField> transformVector
        (
            const vectorField& st
        ) const;

        //- Transform vector using transformation tensor and return
        // symmetrical tensor
        virtual symmTensor transformVector(const vector& st) const;


    // Write

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
