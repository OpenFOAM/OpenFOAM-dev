/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "clouds.H"
#include "basicSpecieMixture.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(clouds, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            clouds,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::clouds::clouds
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(sourceName, modelType, dict, mesh),
    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector(dimAcceleration, Zero)
    ),
    carrierHasThermo_
    (
        mesh.foundObject<fluidThermo>(physicalProperties::typeName)
    ),
    tCarrierThermo_
    (
        carrierHasThermo_
      ? tmpNrc<fluidThermo>
        (
            mesh.lookupObject<fluidThermo>(physicalProperties::typeName)
        )
      : tmpNrc<fluidThermo>(nullptr)
    ),
    tCarrierViscosity_
    (
        carrierHasThermo_
      ? tmpNrc<viscosityModel>(nullptr)
      : tmpNrc<viscosityModel>
        (
            mesh.lookupObject<viscosityModel>(physicalProperties::typeName)
        )
    ),
    tRho_
    (
        carrierHasThermo_
      ? tmp<volScalarField>(nullptr)
      : tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("rho", dimDensity, tCarrierViscosity_())
            )
        )
    ),
    tMu_
    (
        carrierHasThermo_
      ? tmp<volScalarField>(nullptr)
      : tmp<volScalarField>
        (
            new volScalarField("mu", tRho_()*tCarrierViscosity_().nu())
        )
    ),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    cloudsPtr_
    (
        carrierHasThermo_
      ? new parcelCloudList
        (
            mesh.lookupObject<volScalarField>(rhoName_),
            mesh.lookupObject<volVectorField>(UName_),
            g_,
            tCarrierThermo_()
        )
      : new parcelCloudList
        (
            tRho_(),
            mesh.lookupObject<volVectorField>(UName_),
            tMu_(),
            g_
        )
    ),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::clouds::addSupFields() const
{
    wordList fieldNames(1, UName_);

    if (carrierHasThermo_)
    {
        const fluidThermo& carrierThermo = tCarrierThermo_();

        fieldNames.append(rhoName_);

        fieldNames.append(carrierThermo.he().name());

        if (isA<basicSpecieMixture>(carrierThermo))
        {
            const basicSpecieMixture& composition =
                refCast<const basicSpecieMixture>(carrierThermo);

            const PtrList<volScalarField>& Y = composition.Y();

            forAll(Y, i)
            {
                if (composition.solve(i))
                {
                    fieldNames.append(Y[i].name());
                }
            }
        }
    }

    return fieldNames;
}


void Foam::fv::clouds::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    if (!carrierHasThermo_)
    {
        tMu_.ref() = tRho_()*tCarrierViscosity_().nu();
    }

    cloudsPtr_().evolve();

    curTimeIndex_ = mesh().time().timeIndex();
}


void Foam::fv::clouds::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (!carrierHasThermo_)
    {
        FatalErrorInFunction
            << "Applying source to compressible equation when carrier thermo "
            << "is not available"
            << exit(FatalError);
    }

    if (fieldName == rhoName_)
    {
        eqn += cloudsPtr_().Srho(eqn.psi());
    }
}


void Foam::fv::clouds::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (!carrierHasThermo_)
    {
        FatalErrorInFunction
            << "Applying source to compressible equation when carrier thermo "
            << "is not available"
            << exit(FatalError);
    }

    const fluidThermo& carrierThermo = tCarrierThermo_();

    if (fieldName == rhoName_)
    {
        eqn += cloudsPtr_().Srho(eqn.psi());
    }
    else if (fieldName == carrierThermo.he().name())
    {
        eqn += cloudsPtr_().Sh(eqn.psi());
    }
    else if
    (
        isA<basicSpecieMixture>(carrierThermo)
     && refCast<const basicSpecieMixture>(carrierThermo).contains
        (
            eqn.psi().name()
        )
    )
    {
        eqn += cloudsPtr_().SYi
        (
            refCast<const basicSpecieMixture>(carrierThermo).index(eqn.psi()),
            eqn.psi()
        );
    }
}


void Foam::fv::clouds::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (carrierHasThermo_)
    {
        FatalErrorInFunction
            << "Applying source to incompressible equation when carrier thermo "
            << "is available"
            << exit(FatalError);
    }

    if (fieldName == UName_)
    {
        eqn += cloudsPtr_().SU(eqn.psi())/tRho_();
    }
}


void Foam::fv::clouds::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (!carrierHasThermo_)
    {
        FatalErrorInFunction
            << "Applying source to compressible equation when carrier thermo "
            << "is not available"
            << exit(FatalError);
    }

    if (fieldName == UName_)
    {
        eqn += cloudsPtr_().SU(eqn.psi());
    }
}


void Foam::fv::clouds::preUpdateMesh()
{
    // Store the particle positions
    cloudsPtr_().storeGlobalPositions();
}


void Foam::fv::clouds::updateMesh(const mapPolyMesh&)
{}


void Foam::fv::clouds::distribute(const mapDistributePolyMesh& map)
{
    cloudsPtr_().distribute(map);
}


bool Foam::fv::clouds::movePoints()
{
    return true;
}


// ************************************************************************* //
