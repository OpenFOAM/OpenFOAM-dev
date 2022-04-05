/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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
    Foam::RBD::joints::rotating

Description
    Joint with a specified rotational speed

    Reference:
    \verbatim
        Featherstone, R. (2008).
        Rigid body dynamics algorithms.
        Springer.
        Chapter 4.
    \endverbatim

SourceFiles
    rotating.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_joints_rotating_H
#define RBD_joints_rotating_H

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
                         Class rotating Declaration
\*---------------------------------------------------------------------------*/

class rotating
:
    public joint
{
private:

    // Private data

        //- Rotation vector function
        autoPtr<Function1<vector>> omega_;


public:

    //- Runtime type information
    TypeName("rotating");


    // Constructors

        //- Construct for given model from dictionary
        rotating(const rigidBodyModel& model, const dictionary& dict);

        //- Clone this joint
        virtual autoPtr<joint> clone() const;


    //- Destructor
    virtual ~rotating();


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
