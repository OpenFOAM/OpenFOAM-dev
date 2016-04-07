/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "joint.H"
#include "rigidBodyModel.H"
#include "Rs.H"
#include "Rzyx.H"
#include "Pxyz.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
    defineTypeNameAndDebug(joint, 0);
    defineRunTimeSelectionTable(joint, dictionary);
}
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joint::New(joint* jointPtr)
{
    return autoPtr<joint>(jointPtr);
}


Foam::autoPtr<Foam::RBD::joint> Foam::RBD::joint::New
(
    const dictionary& dict
)
{
    const word bodyType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(bodyType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown joint type "
            << bodyType << nl << nl
            << "Valid joint types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<joint>(cstrIter()(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::joint::~joint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::PtrList<Foam::RBD::joint> Foam::RBD::joint::floating()
{
    PtrList<joint> cj(2);
    cj.set(0, new joints::Pxyz());
    //cj.set(1, new joints::Rs());
    cj.set(1, new joints::Rzyx());
    return cj;
}


void Foam::RBD::joint::write(Ostream& os) const
{}


// ************************************************************************* //
