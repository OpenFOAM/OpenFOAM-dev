/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "surfaceInterpolate.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(surfaceInterpolate, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        surfaceInterpolate,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::surfaceInterpolate::surfaceInterpolate
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::surfaceInterpolate::~surfaceInterpolate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::surfaceInterpolate::read
(
    const dictionary& dict
)
{
    dict.lookup("fields") >> fieldSet_;

    return true;
}


bool Foam::functionObjects::surfaceInterpolate::execute()
{
    Info<< type() << " " << name() << " write:" << nl;

    // Clear out any previously loaded fields
    ssf_.clear();
    svf_.clear();
    sSpheretf_.clear();
    sSymmtf_.clear();
    stf_.clear();

    interpolateFields<scalar>(ssf_);
    interpolateFields<vector>(svf_);
    interpolateFields<sphericalTensor>(sSpheretf_);
    interpolateFields<symmTensor>(sSymmtf_);
    interpolateFields<tensor>(stf_);

    Info<< endl;

    return true;
}


bool Foam::functionObjects::surfaceInterpolate::write()
{
    Info<< type() << " " << name() << " write:" << nl;

    Info<< "    Writing interpolated surface fields to "
        << obr_.time().timeName() << endl;

    forAll(ssf_, i)
    {
        ssf_[i].write();
    }
    forAll(svf_, i)
    {
        svf_[i].write();
    }
    forAll(sSpheretf_, i)
    {
        sSpheretf_[i].write();
    }
    forAll(sSymmtf_, i)
    {
        sSymmtf_[i].write();
    }
    forAll(stf_, i)
    {
        stf_[i].write();
    }

    return true;
}


// ************************************************************************* //
