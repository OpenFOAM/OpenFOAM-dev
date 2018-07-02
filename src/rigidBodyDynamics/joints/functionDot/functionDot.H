/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::RBD::joints::functionDot

Description
    Joint in which the position is a function of the parent joint's velocity

    Reference:
    \verbatim
        Featherstone, R. (2008).
        Rigid body dynamics algorithms.
        Springer.
        Chapter 4.
    \endverbatim

SourceFiles
    functionDot.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_joints_functionDot_H
#define RBD_joints_functionDot_H

#include "joint.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace joints
{

/*---------------------------------------------------------------------------*\
                         Class functionDot Declaration
\*---------------------------------------------------------------------------*/

class functionDot
:
    public joint
{
private:

    // Private data

        //- Function
        autoPtr<Function1<scalar>> f_;

        //- Difference used to calculate the gradient of the function
        const scalar delta_;


public:

    //- Runtime type information
    TypeName("functionDot");


    // Constructors

        //- Construct for given model from dictionary
        functionDot(const rigidBodyModel& model, const dictionary& dict);

        //- Clone this joint
        virtual autoPtr<joint> clone() const;


    //- Destructor
    virtual ~functionDot();


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
