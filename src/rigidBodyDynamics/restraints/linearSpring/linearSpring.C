/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "linearSpring.H"
#include "rigidBodyModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{
    defineTypeNameAndDebug(linearSpring, 0);

    addToRunTimeSelectionTable
    (
        restraint,
        linearSpring,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraints::linearSpring::linearSpring
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    restraint(name, dict, model)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::restraints::linearSpring::~linearSpring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::restraints::linearSpring::restrain
(
    scalarField& tau,
    Field<spatialVector>& fx
) const
{
    point attachmentPt = bodyPoint(refAttachmentPt_);

    // Current axis of the spring
    vector r = attachmentPt - anchor_;
    scalar magR = mag(r);
    r /= (magR + vSmall);

    // Velocity of the attached end of the spring
    vector v = bodyPointVelocity(refAttachmentPt_).l();

    // Force and moment on the master body including optional damping
    vector force
    (
        (-stiffness_*(magR - restLength_) - damping_*(r & v))*r
    );

    vector moment(attachmentPt ^ force);

    if (model_.debug)
    {
        Info<< " attachmentPt " << attachmentPt
            << " attachmentPt - anchor " << r*magR
            << " spring length " << magR
            << " force " << force
            << " moment " << moment
            << endl;
    }

    // Accumulate the force for the restrained body
    fx[bodyIndex_] += spatialVector(moment, force);
}


bool Foam::RBD::restraints::linearSpring::read
(
    const dictionary& dict
)
{
    restraint::read(dict);

    coeffs_.lookup("anchor") >> anchor_;
    coeffs_.lookup("refAttachmentPt") >> refAttachmentPt_;
    coeffs_.lookup("stiffness") >> stiffness_;
    coeffs_.lookup("damping") >> damping_;
    coeffs_.lookup("restLength") >> restLength_;

    return true;
}


void Foam::RBD::restraints::linearSpring::write
(
    Ostream& os
) const
{
    restraint::write(os);

    writeEntry(os, "anchor", anchor_);

    writeEntry(os, "refAttachmentPt", refAttachmentPt_);

    writeEntry(os, "stiffness", stiffness_);

    writeEntry(os, "damping", damping_);

    writeEntry(os, "restLength", restLength_);
}


// ************************************************************************* //
