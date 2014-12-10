/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "MultiInteraction.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
bool Foam::MultiInteraction<CloudType>::read(const dictionary& dict)
{
    // Count dictionaries

    Info<< "Patch interaction model " << typeName << nl
        << "Executing in turn " << endl;

    label nModels = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            Info<< "    " << iter().name() << endl;

            nModels++;
        }
    }

    models_.setSize(nModels);
    nModels = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            models_.set
            (
                nModels++,
                PatchInteractionModel<CloudType>::New
                (
                    iter().dict(),
                    this->owner()
                )
            );
        }
    }

    oneInteractionOnly_ = Switch(dict.lookup("oneInteractionOnly"));

    if (oneInteractionOnly_)
    {
        Info<< "Stopping upon first model that interacts with particle."
            << nl << endl;
    }
    else
    {
        Info<< "Allowing multiple models to interact."
            << nl << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MultiInteraction<CloudType>::MultiInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName)
{
    read(this->coeffDict());
}


template<class CloudType>
Foam::MultiInteraction<CloudType>::MultiInteraction
(
    const MultiInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    oneInteractionOnly_(pim.oneInteractionOnly_),
    models_(pim.models_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MultiInteraction<CloudType>::~MultiInteraction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::MultiInteraction<CloudType>::active() const
{
    forAll(models_, i)
    {
        if (models_[i].active())
        {
            return true;
        }
    }
    return false;
}


template<class CloudType>
bool Foam::MultiInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    label origFacei = p.face();
    label patchi = pp.index();

    bool interacted = false;

    forAll(models_, i)
    {
        bool myInteracted = models_[i].correct
        (
            p,
            this->owner().pMesh().boundaryMesh()[patchi],
            keepParticle,
            trackFraction,
            tetIs
        );

        if (myInteracted && oneInteractionOnly_)
        {
            break;
        }

        interacted = (interacted || myInteracted);


        // Check if perhaps the interaction model has changed patches
        // (CoincidentBaffleInteraction can do this)

        if (p.face() != origFacei)
        {
            origFacei = p.face();
            patchi = p.patch(p.face());

            // Interaction model has moved particle off wall?
            if (patchi == -1)
            {
                break;
            }
        }
    }

    return interacted;
}


// ************************************************************************* //
