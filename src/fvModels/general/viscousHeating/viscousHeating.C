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

\*---------------------------------------------------------------------------*/

#include "viscousHeating.H"
#include "basicThermo.H"
#include "compressibleMomentumTransportModel.H"
#include "fvMatrices.H"
#include "fvcDiv.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(viscousHeating, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        viscousHeating,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::viscousHeating::readCoeffs(const dictionary& dict)
{
    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::viscousHeating::viscousHeating
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    phaseName_(word::null)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::viscousHeating::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        );

    return wordList(1, thermo.he().name());
}


void Foam::fv::viscousHeating::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const compressible::momentumTransportModel& momentumTransport =
        mesh().lookupType<compressible::momentumTransportModel>();

    volVectorField& U = const_cast<volVectorField&>(momentumTransport.U());

    mesh().schemes().setFluxRequired(U.name());

    eqn -= fvc::div
    (
        fvc::dotInterpolate(momentumTransport.divDevTau(U)->flux(), U)
    );
}


void Foam::fv::viscousHeating::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const compressible::momentumTransportModel& momentumTransport =
        mesh().lookupType<compressible::momentumTransportModel>(phaseName_);

    volVectorField& U = const_cast<volVectorField&>(momentumTransport.U());

    mesh().schemes().setFluxRequired(U.name());

    eqn -= fvc::div
    (
        fvc::dotInterpolate(momentumTransport.divDevTau(U)->flux(), U)
    );
}


bool Foam::fv::viscousHeating::movePoints()
{
    return true;
}


void Foam::fv::viscousHeating::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::viscousHeating::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::viscousHeating::distribute(const polyDistributionMap&)
{}


bool Foam::fv::viscousHeating::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
