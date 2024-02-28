/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2024 OpenFOAM Foundation
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

#include "cylindricalFunctionObject.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cylindrical, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        cylindrical,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tensor Foam::functionObjects::cylindrical::R(const vector& p) const
{
    const vector dir = normalised(p - origin_);

    const vector axis = normalised(axis_);
    const vector r = dir - (dir & axis)*axis;

    return tensor(normalised(r), normalised(axis^r), axis);
}


void Foam::functionObjects::cylindrical::transform
(
    vectorField& vf,
    const vectorField& points
) const
{
    if (toCartesian_)
    {
        forAll(points, i)
        {
            vf[i] = this->R(points[i]).T() & vf[i];
        }
    }
    else
    {
        forAll(points, i)
        {
            vf[i] = this->R(points[i]) & vf[i];
        }
    }
}


bool Foam::functionObjects::cylindrical::calc()
{
    if (foundObject<volVectorField>(fieldName_))
    {
        const volVectorField& vf = lookupObject<volVectorField>(fieldName_);
        const volVectorField& C = mesh_.C();

        tmp<volVectorField> tcvf(volVectorField::New(resultName_, vf));

        volVectorField& cvf = tcvf.ref();
        transform(cvf.primitiveFieldRef(), C);

        forAll(vf.boundaryField(), patchi)
        {
            transform
            (
                cvf.boundaryFieldRef()[patchi],
                C.boundaryField()[patchi]
            );
        }

        return store(resultName_, tcvf);
    }
    else
    {
        cannotFindObject<volVectorField>(fieldName_);

        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cylindrical::cylindrical
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression
    (
        name,
        runTime,
        dict,
        dict.lookupOrDefault<Switch>("toCartesian", false)
      ? "cartesian"
      : typeName
    ),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    toCartesian_(dict.lookupOrDefault<Switch>("toCartesian", false))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cylindrical::~cylindrical()
{}


// ************************************************************************* //
