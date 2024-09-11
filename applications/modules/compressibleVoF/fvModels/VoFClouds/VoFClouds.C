/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "VoFClouds.H"
#include "compressibleTwoPhaseVoFMixture.H"
#include "fvmSup.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(VoFClouds, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            VoFClouds,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::VoFClouds::VoFClouds
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    phaseName_(dict.lookup("phase")),
    carrierPhaseName_(dict.lookup("carrierPhase")),
    thermo_
    (
        mesh.lookupObject<fluidThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        )
    ),
    carrierThermo_
    (
        mesh.lookupObject<fluidThermo>
        (
            IOobject::groupName(physicalProperties::typeName, carrierPhaseName_)
        )
    ),
    clouds_
    (
        dict.lookupOrDefault<wordList>
        (
            "clouds",
            parcelCloudList::defaultCloudNames
        ),
        carrierThermo_.rho(),
        mesh.lookupObject<volVectorField>("U"),
        mesh.lookupObject<uniformDimensionedVectorField>("g"),
        carrierThermo_
    ),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::VoFClouds::addSupFields() const
{
    return wordList({thermo_.rho()().name(), thermo_.he().name(), "U"});
}


void Foam::fv::VoFClouds::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    clouds_.evolve();

    curTimeIndex_ = mesh().time().timeIndex();
}


void Foam::fv::VoFClouds::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&rho == &thermo_.rho()())
    {
        eqn += clouds_.Srho();
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << rho.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFClouds::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&he == &thermo_.he())
    {
        eqn += clouds_.Sh(eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << he.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFClouds::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (U.name() == "U")
    {
        eqn += clouds_.SU(eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << U.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFClouds::preUpdateMesh()
{
    // Store the particle positions
    clouds_.storeGlobalPositions();
}


void Foam::fv::VoFClouds::topoChange(const polyTopoChangeMap& map)
{
    clouds_.topoChange(map);
}


void Foam::fv::VoFClouds::mapMesh(const polyMeshMap& map)
{
    clouds_.mapMesh(map);
}


void Foam::fv::VoFClouds::distribute(const polyDistributionMap& map)
{
    clouds_.distribute(map);
}


bool Foam::fv::VoFClouds::movePoints()
{
    return true;
}


// ************************************************************************* //
