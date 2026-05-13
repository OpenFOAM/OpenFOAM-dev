/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2026 OpenFOAM Foundation
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

#include "PecletNo.H"
#include "incompressibleMomentumTransportModel.H"
#include "compressibleMomentumTransportModel.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(PecletNo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        PecletNo,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::PecletNo::calc()
{
    if (foundObject<surfaceScalarField>(fieldName_))
    {
        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(fieldName_);

        tmp<volScalarField> nuOrMuEff;
        if (phi.dimensions() == dimVolumetricFlux)
        {
            const incompressibleMomentumTransportModel& imtm =
                mesh_.lookupType<incompressibleMomentumTransportModel>();

            nuOrMuEff = imtm.nuEff();
        }
        else if (phi.dimensions() == dimMassFlux)
        {
            const compressibleMomentumTransportModel& cmtm =
                mesh_.lookupType<compressibleMomentumTransportModel>();

            nuOrMuEff = cmtm.rho()*cmtm.nuEff();
        }
        else
        {
            FatalErrorInFunction
                << "dimensions of flux " << phi.name() << " are "
                << phi.dimensions() << ", but they must be either "
                << dimMassFlux << " (mass flux) or " << dimVolumetricFlux
                << " (volume flux) " << exit(FatalError);
        }

        store
        (
            resultName_,
            mag(phi)
           /(
                mesh_.magSf()
               *mesh_.surfaceInterpolation::deltaCoeffs()
               *fvc::interpolate(nuOrMuEff)
            )
        );

        return true;
    }
    else
    {
        cannotFindObject<surfaceScalarField>(fieldName_);

        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::PecletNo::PecletNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "Pe", "phi")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::PecletNo::~PecletNo()
{}


// ************************************************************************* //
