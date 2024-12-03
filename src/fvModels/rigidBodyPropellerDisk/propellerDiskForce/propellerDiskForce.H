/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::RBD::restraints::propellerDiskForce

Description
     rigidBodyPropellerDisk force restraint.

     Applies the opposite force and moment of the propellerDisk onto the body
     the propeller is attached to.

Usage
    Example applying the propeller force and moment to the floatingObject:
    \verbatim
    restraints
    {
        force
        {
            type        propellerDiskForce;
            body        floatingObject;
        }
    }
    \endverbatim

See also
    Foam::fv::rigidBodyPropellerDisk
    Foam::fv::propellerDisk

SourceFiles
    propellerDiskForce.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_restraints_propellerDiskForce_H
#define RBD_restraints_propellerDiskForce_H

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
                        Class propellerDiskForce Declaration
\*---------------------------------------------------------------------------*/

class propellerDiskForce
:
    public restraint
{

public:

    //- Runtime type information
    TypeName("propellerDiskForce");


    // Constructors

        //- Construct from components
        propellerDiskForce
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
                new propellerDiskForce(*this)
            );
        }


    //- Destructor
    virtual ~propellerDiskForce();


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
