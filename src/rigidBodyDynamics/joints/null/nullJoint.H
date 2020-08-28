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
    Foam::RBD::joints::null

Description
    Null joint for the root-body.

    Reference:
    \verbatim
        Featherstone, R. (2008).
        Rigid body dynamics algorithms.
        Springer.
        Chapter 4.
    \endverbatim

SourceFiles
    nullJoint.C

\*---------------------------------------------------------------------------*/

#ifndef nullJoint_H
#define nullJoint_H

#include "joint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

class rigidBodyModel;

namespace joints
{

/*---------------------------------------------------------------------------*\
                            Class null Declaration
\*---------------------------------------------------------------------------*/

class null
:
    public joint
{

public:

    //- Runtime type information
    TypeName("null");


    // Constructors

        //- Construct for given model
        null(const rigidBodyModel& model);

        //- Construct for given model from dictionary
        null(const rigidBodyModel& model, const dictionary& dict);

        //- Clone this joint
        virtual autoPtr<joint> clone() const;


    //- Destructor
    virtual ~null();


    // Member Functions

        //- Update the model state for this joint
        virtual void jcalc
        (
            joint::XSvc& J,
            const rigidBodyModelState& state
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace joints
} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
