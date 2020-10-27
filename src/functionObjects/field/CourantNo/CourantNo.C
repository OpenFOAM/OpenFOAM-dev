/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "CourantNo.H"
#include "surfaceFields.H"
#include "fvcSurfaceIntegrate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(CourantNo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        CourantNo,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::functionObjects::CourantNo::byRho
(
    const tmp<volScalarField::Internal>& Co
) const
{
    if (Co().dimensions() == dimDensity)
    {
        return Co/obr_.lookupObject<volScalarField>(rhoName_);
    }
    else
    {
        return Co;
    }
}


bool Foam::functionObjects::CourantNo::calc()
{
    if (foundObject<surfaceScalarField>(fieldName_))
    {
        const surfaceScalarField& phi =
            lookupObject<surfaceScalarField>(fieldName_);

        tmp<volScalarField> tCo
        (
            volScalarField::New
            (
                resultName_,
                mesh_,
                dimensionedScalar(dimless, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        tCo->ref() =
            byRho
            (
                (0.5*mesh_.time().deltaT())
               *fvc::surfaceSum(mag(phi))()()
               /mesh_.V()
            );

        tCo->correctBoundaryConditions();

        return store(resultName_, tCo);
    }
    else
    {
        cannotFindObject<surfaceScalarField>(fieldName_);

        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::CourantNo::CourantNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "Co", "phi")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::CourantNo::~CourantNo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CourantNo::read(const dictionary& dict)
{
    fieldExpression::read(dict);

    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

    return true;
}


// ************************************************************************* //
