/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "blockMeshConfigurationBase.H"
#include "polyPatch.H"
#include "wallPolyPatch.H"
#include "blockMeshFunctions.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::blockMeshConfigurationBase::roundBoundingBox
(
    boundBox& bb,
    const scalar s
)
{
    Info<< "Rounding bounding box to multiples of " << s << endl;

    bb.min() = roundDown(bb.min(), s);
    bb.max() = roundUp(bb.max(), s);
}


Foam::Pair<Foam::word> Foam::blockMeshConfigurationBase::readPatchOption
(
    const word& option
) const
{
    const Pair<word> patchOpt(patchOpts_.find(option)());

    if
    (
        !(
            polyPatch::constraintType(patchOpt.second())
         || patchOpt.second() == wallPolyPatch::typeName
         || patchOpt.second() == polyPatch::typeName
        )
    )
    {
        FatalErrorInFunction<< "Argument to '-"
            << option << " should be of the form "
            << "'(<name> <type>)'" << nl
            << "where <type> must be a generic \"patch\", \"wall\" "
            << "or a constraint condition:" << nl << nl
            << polyPatch::constraintTypes()
            << exit(FatalError);
    }

    return patchOpt;
}


void Foam::blockMeshConfigurationBase::writeVertex
(
    const word& x,
    const word& y,
    const word& z
)
{
    const word bgm("$!backgroundMesh/");
    os_ << indent << "("
        << bgm << x << " "
        << bgm << y << " "
        << bgm << z << ")"
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMeshConfigurationBase::blockMeshConfigurationBase
(
    const fileName& name,
    const fileName& dir,
    const Time& time,
    const meshingSurfaceList& surfaces,
    const HashTable<Pair<word>>& patchOpts
)
:
    caseFileConfiguration(name, dir, time),
    bb_(surfaces.bb()),
    patchOpts_(patchOpts)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockMeshConfigurationBase::~blockMeshConfigurationBase()
{}


// ************************************************************************* //
