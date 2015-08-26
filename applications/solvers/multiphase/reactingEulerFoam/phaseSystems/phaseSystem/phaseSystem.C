/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "phaseSystem.H"
#include "surfaceTensionModel.H"
#include "aspectRatioModel.H"
#include "surfaceInterpolate.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
}

const Foam::word Foam::phaseSystem::propertiesName("phaseProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::phaseSystem::phaseModelTable
Foam::phaseSystem::generatePhaseModels(const wordList& phaseNames) const
{
    phaseModelTable phaseModels;

    forAllConstIter(wordList, phaseNames, phaseNameIter)
    {
        phaseModels.insert
        (
            *phaseNameIter,
            phaseModel::New
            (
                *this,
                *phaseNameIter
            )
        );
    }

    // normalise ?

    return phaseModels;
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::generatePhi
(
    const phaseModelTable& phaseModels
) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels.begin();

    tmp<surfaceScalarField> tmpPhi
    (
        new surfaceScalarField
        (
            "phi",
            fvc::interpolate(phaseModelIter()())*phaseModelIter()->phi()
        )
    );

    ++phaseModelIter;

    for (; phaseModelIter != phaseModels.end(); ++ phaseModelIter)
    {
        tmpPhi() +=
            fvc::interpolate(phaseModelIter()())
           *phaseModelIter()->phi();
    }

    return tmpPhi;
}


void Foam::phaseSystem::generatePairs
(
    const dictTable& modelDicts
)
{
    forAllConstIter(dictTable, modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        // pair already exists
        if (phasePairs_.found(key))
        {
            // do nothing ...
        }

        // new ordered pair
        else if (key.ordered())
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new orderedPhasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }

        // new unordered pair
        else
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    phaseNames_(lookup("phases")),

    phaseModels_(generatePhaseModels(phaseNames_)),

    phi_(generatePhi(phaseModels_)),

    dpdt_
    (
        IOobject
        (
            "dpdt",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dpdt", dimPressure/dimTime, 0)
    ),

    MRF_(mesh_),
    fvOptions_(mesh_)
{
    // Blending methods
    forAllConstIter(dictionary, subDict("blending"), iter)
    {
        blendingMethods_.insert
        (
            iter().dict().dictName(),
            blendingMethod::New
            (
                iter().dict(),
                wordList(lookup("phases"))
            )
        );
    }

    // Sub-models
    generatePairsAndSubModels("surfaceTension", surfaceTensionModels_);
    generatePairsAndSubModels("aspectRatio", aspectRatioModels_);

    correctKinematics();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseSystem::rho() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpRho
    (
        phaseModelIter()()*phaseModelIter()->rho()
    );

    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpRho() += phaseModelIter()()*phaseModelIter()->rho();
    }

    return tmpRho;
}


Foam::tmp<Foam::volVectorField> Foam::phaseSystem::U() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volVectorField> tmpU
    (
        phaseModelIter()()*phaseModelIter()->U()
    );

    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpU() += phaseModelIter()()*phaseModelIter()->U();
    }

    return tmpU;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::sigma(const phasePairKey& key) const
{
    if (surfaceTensionModels_.found(key))
    {
        return surfaceTensionModels_[key]->sigma();
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    surfaceTensionModel::typeName + ":sigma",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar("zero", surfaceTensionModel::dimSigma, 0)
            )
        );
    }
}


void Foam::phaseSystem::solve()
{}


void Foam::phaseSystem::correct()
{
    forAllIter(phaseModelTable, phaseModels_, phaseModelIter)
    {
        phaseModelIter()->correct();
    }
}


void Foam::phaseSystem::correctKinematics()
{
    bool updateDpdt = false;

    forAllIter(phaseModelTable, phaseModels_, phaseModelIter)
    {
        phaseModelIter()->correctKinematics();

        updateDpdt = updateDpdt || phaseModelIter()->thermo().dpdt();
    }

    // Update the pressure time-derivative if required
    if (updateDpdt)
    {
        dpdt_ = fvc::ddt(phaseModels_.begin()()().thermo().p());
    }
}


void Foam::phaseSystem::correctThermo()
{
    forAllIter(phaseModelTable, phaseModels_, phaseModelIter)
    {
        phaseModelIter()->correctThermo();
    }
}


void Foam::phaseSystem::correctTurbulence()
{
    forAllIter(phaseModelTable, phaseModels_, phaseModelIter)
    {
        phaseModelIter()->correctTurbulence();
    }
}


bool Foam::phaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        forAllIter(phaseModelTable, phaseModels_, phaseModelIter)
        {
            readOK &= phaseModelIter()->read();
        }

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
