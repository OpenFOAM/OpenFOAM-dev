/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "adjustTimeStepToNucleation.H"
#include "fvModels.H"
#include "nucleation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(adjustTimeStepToNucleation, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        adjustTimeStepToNucleation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::adjustTimeStepToNucleation::adjustTimeStepToNucleation
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    modelName_(word::null),
    maxCo_(NaN)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::adjustTimeStepToNucleation::~adjustTimeStepToNucleation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::adjustTimeStepToNucleation::read
(
    const dictionary& dict
)
{
    modelName_ = dict.lookupOrDefault<word>("model", word::null);
    maxCo_ = dict.lookupOrDefault<scalar>("maxCo", 1);

    return true;
}


bool Foam::functionObjects::adjustTimeStepToNucleation::execute()
{
    return true;
}


bool Foam::functionObjects::adjustTimeStepToNucleation::write()
{
    return true;
}


Foam::scalar
Foam::functionObjects::adjustTimeStepToNucleation::maxDeltaT() const
{
    if (!time_.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        return vGreat;
    }

    const Foam::fvModels& fvModels = Foam::fvModels::New(mesh_);

    if (modelName_ == word::null)
    {
        bool found = false;

        tmp<volScalarField::Internal> tTau =
            volScalarField::Internal::New
            (
                typedName("tau"),
                mesh(),
                dimensionedScalar(dimTime, vGreat)
            );

        forAll(fvModels, fvModeli)
        {
            if (isA<fv::nucleation>(fvModels[fvModeli]))
            {
                found = true;

                const fv::nucleation& nucleationModel =
                    refCast<const fv::nucleation>(fvModels[fvModeli]);

                tTau = min(tTau, nucleationModel.tau());
            }
        }

        if (!found)
        {
            WarningInFunction
                << "No nucleation models found"
                << exit(FatalError);
        }

        return gMin(tTau());
    }
    else
    {
        const fv::nucleation& nucleationModel =
            refCast<const fv::nucleation>(fvModels[modelName_]);

        return gMin(nucleationModel.tau()());
    }
}


// ************************************************************************* //
