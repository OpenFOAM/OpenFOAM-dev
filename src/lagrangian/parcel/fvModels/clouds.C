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
    carrierThermo_
    (
        mesh.lookupObject<fluidThermo>(basicThermo::dictName)
    ),
    clouds_
    (
        mesh.lookupObject<volScalarField>("rho"),
        mesh.lookupObject<volVectorField>("U"),
        mesh.lookupObject<uniformDimensionedVectorField>("g"),
        carrierThermo_
    ),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::clouds::addSupFields() const
{
    wordList fieldNames({"rho", "U", carrierThermo_.he().name()});

    if (isA<basicSpecieMixture>(carrierThermo_))
    {
        const basicSpecieMixture& composition =
            refCast<const basicSpecieMixture>(carrierThermo_);

        const PtrList<volScalarField>& Y = composition.Y();

        forAll(Y, i)
        {
            if (composition.solve(i))
            {
                fieldNames.append(Y[i].name());
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

    clouds_.evolve();

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

    if (fieldName == "rho")
    {
        eqn += clouds_.Srho(eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
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

    if (fieldName == "rho")
    {
        eqn += clouds_.Srho(eqn.psi());
    }
    else if (fieldName == carrierThermo_.he().name())
    {
        eqn += clouds_.Sh(eqn.psi());
    }
    else if
    (
        isA<basicSpecieMixture>(carrierThermo_)
     && refCast<const basicSpecieMixture>(carrierThermo_).contains
        (
            eqn.psi().name()
        )
    )
    {
        eqn += clouds_.SYi
        (
            refCast<const basicSpecieMixture>(carrierThermo_).index(eqn.psi()),
            eqn.psi()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
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

    if (fieldName == "U")
    {
        eqn += clouds_.SU(eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::clouds::preUpdateMesh()
{
    // Store the particle positions
    clouds_.storeGlobalPositions();
}


// ************************************************************************* //
