/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2021 OpenFOAM Foundation
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

#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicSpecieMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSpecieMixture::basicSpecieMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicMixture(thermoDict, mesh, phaseName),
    phaseName_(phaseName),
    species_(specieNames),
    defaultSpecie_
    (
        species_.size()
      ? thermoDict.lookupBackwardsCompatible<word>
        (
            {"defaultSpecie", "inertSpecie"}
        )
      : word("undefined")
    ),
    defaultSpecieIndex_(-1),
    active_(species_.size(), true),
    Y_(species_.size())
{
    if (species_.size())
    {
        if (species_.found(defaultSpecie_))
        {
            defaultSpecieIndex_ = species_[defaultSpecie_];
        }
        else
        {
            FatalIOErrorInFunction(thermoDict)
                << "default specie " << defaultSpecie_
                << " not found in available species " << species_
                << exit(FatalIOError);
        }
    }

    tmp<volScalarField> tYdefault;

    forAll(species_, i)
    {
        typeIOobject<volScalarField> header
        (
            IOobject::groupName(species_[i], phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.headerOk())
        {
            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], phaseName),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            // Read Ydefault if not already read
            if (!tYdefault.valid())
            {
                const word YdefaultName
                (
                    IOobject::groupName("Ydefault", phaseName)
                );

                typeIOobject<volScalarField> timeIO
                (
                    YdefaultName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                typeIOobject<volScalarField> constantIO
                (
                    YdefaultName,
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                typeIOobject<volScalarField> time0IO
                (
                    YdefaultName,
                    Time::timeName(0),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (timeIO.headerOk())
                {
                    tYdefault = new volScalarField(timeIO, mesh);
                }
                else if (constantIO.headerOk())
                {
                    tYdefault = new volScalarField(constantIO, mesh);
                }
                else
                {
                    tYdefault = new volScalarField(time0IO, mesh);
                }
            }

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], phaseName),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tYdefault()
                )
            );
        }
    }

    correctMassFractions();
}


void Foam::basicSpecieMixture::correctMassFractions()
{
    if (species_.size())
    {
        tmp<volScalarField> tYt
        (
            volScalarField::New
            (
                IOobject::groupName("Yt", phaseName_),
                Y_[0],
                calculatedFvPatchScalarField::typeName
            )
        );
        volScalarField& Yt = tYt.ref();

        for (label i=1; i<Y_.size(); i++)
        {
            Yt += Y_[i];
        }

        if (mag(max(Yt).value()) < rootVSmall)
        {
            FatalErrorInFunction
                << "Sum of mass fractions is zero for species " << species()
                << exit(FatalError);
        }

        forAll(Y_, i)
        {
            Y_[i] /= Yt;
        }
    }
}


void Foam::basicSpecieMixture::normalise()
{
    if (species_.size())
    {
        tmp<volScalarField> tYt
        (
            volScalarField::New
            (
                IOobject::groupName("Yt", phaseName_),
                Y_[0].mesh(),
                dimensionedScalar(dimless, 0)
            )
        );
        volScalarField& Yt = tYt.ref();

        forAll(Y_, i)
        {
            if (solve(i))
            {
                Y_[i].max(scalar(0));
                Yt += Y_[i];
            }
        }

        Y_[defaultSpecieIndex_] = scalar(1) - Yt;
        Y_[defaultSpecieIndex_].max(scalar(0));
    }
}


// ************************************************************************* //
