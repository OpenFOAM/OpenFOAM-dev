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
    Foam::fv::rigidBodyPropellerDisk

Description
    Disk momentum source derived from Foam::fv::propellerDisk with support
    for motion of the body the propeller is attached to.

Usage
    Example usage:
    \verbatim
    diskSource
    {
        type            rigidBodyPropellerDisk;

        libs            ("librigidBodyPropellerDisk.so");

        selectionMode   cellZone;
        cellZone        propeller;

        body            hull;

        centre          (0 0 0);    // Centre of the propeller
        normal          (1 0 0);    // Normal direction of the propeller

        n               26.03;      // Rotation speed [1/s]

        dPropeller      0.203;      // Propeller diameter
        dHub            0.039179;   // Hub diameter

        propellerCurve
        {
            type table;

            //   J     Kt     Kq
            values
            (
                (0.10 (0.3267 0.03748))
                (0.15 (0.3112 0.03629))
                (0.20 (0.2949 0.03500))
                (0.25 (0.2777 0.03361))
                (0.30 (0.2598 0.03210))
                (0.35 (0.2410 0.03047))
                (0.40 (0.2214 0.02871))
                (0.45 (0.2010 0.02682))
                (0.50 (0.1798 0.02479))
                (0.55 (0.1577 0.02261))
                (0.60 (0.1349 0.02027))
                (0.65 (0.1112 0.01777))
                (0.70 (0.0867 0.01509))
                (0.75 (0.0614 0.01224))
                (0.80 (0.0353 0.00921))
            );
        }
    }
    \endverbatim

See also
    Foam::fv::propellerDisk

SourceFiles
    rigidBodyPropellerDisk.C
    rigidBodyPropellerDiskTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyPropellerDisk_H
#define rigidBodyPropellerDisk_H

#include "propellerDisk.H"
#include "rigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class rigidBodyPropellerDisk Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyPropellerDisk
:
    public propellerDisk
{
    // Private Member Data

        const RBD::rigidBodyMotion& motion_;

        //- Name of the body
        word body_;

        //- ID of the body
        label bodyID_;


protected:

    // Protected Member Functions

        //- Read the model coefficients
        virtual void readCoeffs(const dictionary& dict);

        //- Return the current propeller normal
        virtual vector normal() const;


public:

    //- Runtime type information
    TypeName("rigidBodyPropellerDisk");


    // Constructors

        //- Construct from components
        rigidBodyPropellerDisk
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rigidBodyPropellerDisk(const rigidBodyPropellerDisk&) = delete;


    //- Destructor
    virtual ~rigidBodyPropellerDisk()
    {}


    // Member Functions

        //- Return the current propeller centre
        virtual vector centre() const;

        label bodyID() const
        {
            return bodyID_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
