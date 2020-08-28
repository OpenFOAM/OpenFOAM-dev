/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::RBD::restraints::externalForce

Description
    Time-dependent external force restraint using Function1.

Usage
    Example applying a constant force to the floatingObject:
    \verbatim
    restraints
    {
        force
        {
            type        externalForce;
            body        floatingObject;
            location    (0 0 0);
            force       (100 0 0);
        }
    }
    \endverbatim

SourceFiles
    externalForce.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_restraints_externalForce_H
#define RBD_restraints_externalForce_H

#include "rigidBodyRestraint.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{

/*---------------------------------------------------------------------------*\
                        Class externalForce Declaration
\*---------------------------------------------------------------------------*/

class externalForce
:
    public restraint
{
    // Private data

        //- External force (N)
        autoPtr<Function1<vector>> externalForce_;

        //- Point of application of the force
        vector location_;


public:

    //- Runtime type information
    TypeName("externalForce");


    // Constructors

        //- Construct from components
        externalForce
        (
            const word& name,
            const dictionary& dict,
            const rigidBodyModel& model
        );

        //- Construct and return a clone
        virtual autoPtr<restraint> clone() const
        {
            return autoPtr<restraint>
            (
                new externalForce(*this)
            );
        }


    //- Destructor
    virtual ~externalForce();


    // Member Functions

        //- Accumulate the retraint internal joint forces into the tau field and
        //  external forces into the fx field
        virtual void restrain
        (
            scalarField& tau,
            Field<spatialVector>& fx,
            const rigidBodyModelState& state
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dict);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace restraints
} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
