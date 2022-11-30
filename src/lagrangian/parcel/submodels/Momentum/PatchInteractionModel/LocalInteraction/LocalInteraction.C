/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "LocalInteraction.H"
#include "wordAndDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct wordReAndDictionary
:
    public Tuple2<wordRe, dictionary>
{
    using Tuple2<wordRe, dictionary>::Tuple2;

    wordReAndDictionary();

    wordReAndDictionary(Istream& is);
};


inline Istream& operator>>(Istream& is, wordReAndDictionary& wd)
{
    wd.first() = wordRe(is);
    dictionary d(is);
    wd.second().transfer(d);
    return is;
}


inline Ostream& operator<<(Ostream& os, const wordReAndDictionary& wd)
{
    return os << wd.first() << token::SPACE << wd.second();
}


inline wordReAndDictionary::wordReAndDictionary()
{}


inline wordReAndDictionary::wordReAndDictionary(Istream& is)
{
    is >> *this;
}

}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LocalInteraction<CloudType>::LocalInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    patchInteractionTypes_
    (
        this->owner().mesh().boundaryMesh().size(),
        PatchInteractionModel<CloudType>::itOther
    ),
    patchEs_(this->owner().mesh().boundaryMesh().size(), NaN),
    patchMus_(this->owner().mesh().boundaryMesh().size(), NaN),
    nEscape_(this->owner().mesh().boundaryMesh().size(), 0),
    massEscape_(this->owner().mesh().boundaryMesh().size(), scalar(0)),
    nStick_(this->owner().mesh().boundaryMesh().size(), 0),
    massStick_(this->owner().mesh().boundaryMesh().size(), scalar(0)),
    writeFields_(this->coeffDict().lookupOrDefault("writeFields", false)),
    massEscapePtr_(nullptr),
    massStickPtr_(nullptr)
{
    const polyBoundaryMesh& patches = this->owner().mesh().boundaryMesh();

    if (writeFields_)
    {
        const word massEscapeName(this->owner().name() + ":massEscape");
        const word massStickName(this->owner().name() + ":massStick");

        Info<< "    Interaction fields will be written to " << massEscapeName
            << " and " << massStickName << endl;

        (void)massEscape();
        (void)massStick();
    }
    else
    {
        Info<< "    Interaction fields will not be written" << endl;
    }

    // Get the patch-settings dictionaries
    dictionary patchesDict;
    if (this->coeffDict().isDict("patches"))
    {
        patchesDict = this->coeffDict().subDict("patches");
    }
    else
    {
        const List<wordReAndDictionary> patchNameAndDicts
        (
            this->coeffDict().lookup("patches")
        );

        forAll(patchNameAndDicts, dicti)
        {
            patchesDict.set
            (
                keyType(string(patchNameAndDicts[dicti].first())),
                patchNameAndDicts[dicti].second()
            );
        }
    }

    // Read the patch settings
    wordList unspecifiedNonConstraintPatches;
    forAll(patches, patchi)
    {
        const word& patchName = patches[patchi].name();

        const bool havePatchDict = patchesDict.found(patchName);

        const bool patchIsConstraint =
            polyPatch::constraintType(patches[patchi].type());

        // No settings for constrained patch. No model.
        if (!havePatchDict && patchIsConstraint)
        {
            patchInteractionTypes_[patchi] =
                PatchInteractionModel<CloudType>::itNone;
            continue;
        }

        // No settings for non-constrained patch. Error.
        if (!havePatchDict && !patchIsConstraint)
        {
            unspecifiedNonConstraintPatches.append(patches[patchi].name());
            continue;
        }

        const dictionary& patchDict = patchesDict.subDict(patchName);

        // Settings for constrained patch. Ignored unless "patchType" is
        // correctly specified.
        if (havePatchDict && patchIsConstraint)
        {
            if (!patchDict.found("patchType"))
            {
                patchInteractionTypes_[patchi] =
                    PatchInteractionModel<CloudType>::itNone;
                continue;
            }

            const word patchType = patchDict.lookup<word>("patchType");

            if (patchType != patches[patchi].type())
            {
                FatalErrorInFunction
                    << "Type " << patchType
                    << " specified for patch " << patchName
                    << " does not match the patch type "
                    << patches[patchi].type() << exit(FatalError);
            }
        }

        // Read and set the interaction model
        const word itName = patchDict.lookup<word>("type");
        const interactionType it = this->wordToInteractionType(itName);

        if (it == PatchInteractionModel<CloudType>::itOther)
        {
            FatalErrorInFunction
                << "Unknown patch interaction type "
                << itName << " for patch " << patchName
                << ". Valid types are:"
                << PatchInteractionModel<CloudType>::interactionTypeNames_
                << nl << exit(FatalError);
        }

        patchInteractionTypes_[patchi] = it;

        if (it == PatchInteractionModel<CloudType>::itRebound)
        {
            patchEs_[patchi] = patchDict.lookupOrDefault<scalar>("e", 1);
            patchMus_[patchi] = patchDict.lookupOrDefault<scalar>("mu", 0);
        }
    }

    // Error if interactions are unspecified for non-constraint patches
    if (!unspecifiedNonConstraintPatches.empty())
    {
        FatalErrorInFunction
            << "No interaction type was specified for non-constraint patches: "
            << unspecifiedNonConstraintPatches
            << exit(FatalError);
    }
}


template<class CloudType>
Foam::LocalInteraction<CloudType>::LocalInteraction
(
    const LocalInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    patchInteractionTypes_(pim.patchInteractionTypes_),
    patchEs_(pim.patchEs_),
    patchMus_(pim.patchMus_),
    nEscape_(pim.nEscape_),
    massEscape_(pim.massEscape_),
    nStick_(pim.nStick_),
    massStick_(pim.massStick_),
    writeFields_(pim.writeFields_),
    massEscapePtr_(nullptr),
    massStickPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LocalInteraction<CloudType>::~LocalInteraction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::volScalarField& Foam::LocalInteraction<CloudType>::massEscape()
{
    if (!massEscapePtr_.valid())
    {
        const fvMesh& mesh = this->owner().mesh();

        massEscapePtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + ":massEscape",
                    mesh.time().name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimMass, 0)
            )
        );
    }

    return massEscapePtr_();
}


template<class CloudType>
Foam::volScalarField& Foam::LocalInteraction<CloudType>::massStick()
{
    if (!massStickPtr_.valid())
    {
        const fvMesh& mesh = this->owner().mesh();

        massStickPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + ":massStick",
                    mesh.time().name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimMass, 0)
            )
        );
    }

    return massStickPtr_();
}


template<class CloudType>
bool Foam::LocalInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    const label patchi = pp.index();

    switch (patchInteractionTypes_[patchi])
    {
        case PatchInteractionModel<CloudType>::itNone:
        {
            return false;
        }

        case PatchInteractionModel<CloudType>::itEscape:
        {
            const scalar dm = p.mass()*p.nParticle();

            keepParticle = false;
            p.moving() = false;
            p.U() = Zero;
            nEscape_[patchi] ++;
            massEscape_[patchi] += dm;

            if (writeFields_)
            {
                const label patchFacei = pp.whichFace(p.face());
                massEscape().boundaryFieldRef()[patchi][patchFacei] += dm;
            }

            return true;
        }

        case PatchInteractionModel<CloudType>::itStick:
        {
            const scalar dm = p.mass()*p.nParticle();

            keepParticle = true;
            p.moving() = false;
            p.U() = Zero;
            nStick_[patchi] ++;
            massStick_[patchi] += dm;

            if (writeFields_)
            {
                const label patchFacei = pp.whichFace(p.face());
                massStick().boundaryFieldRef()[patchi][patchFacei] += dm;
            }

            return true;
        }

        case PatchInteractionModel<CloudType>::itRebound:
        {
            keepParticle = true;
            p.moving() = true;

            vector nw, Up;
            this->owner().patchData(p, pp, nw, Up);

            // Make motion relative to patch velocity
            p.U() -= Up;

            const scalar Un = p.U() & nw;
            const vector Ut = p.U() - Un*nw;

            if (Un > 0)
            {
                p.U() -= (1 + patchEs_[patchi])*Un*nw;
            }

            p.U() -= patchMus_[patchi]*Ut;

            // Return velocity to global space
            p.U() += Up;

            return true;
        }

        default:
        {
            return false;
        }
    }

    return false;
}


template<class CloudType>
void Foam::LocalInteraction<CloudType>::info(Ostream& os)
{
    const polyBoundaryMesh& patches = this->owner().mesh().boundaryMesh();

    // Determine the number of non-processor patches
    label nPatches = patches.size();
    for (; isA<processorPolyPatch>(patches[nPatches - 1]); nPatches --);

    // Retrieve any stored data
    labelList npe0(nPatches, 0);
    this->getModelProperty("nEscape", npe0);

    scalarList mpe0(nPatches, scalar(0));
    this->getModelProperty("massEscape", mpe0);

    labelList nps0(nPatches, 0);
    this->getModelProperty("nStick", nps0);

    scalarList mps0(nPatches, scalar(0));
    this->getModelProperty("massStick", mps0);

    // Accumulate current data
    labelList npe(SubList<label>(nEscape_, nPatches));
    Pstream::listCombineGather(npe, plusEqOp<label>());
    npe = npe + npe0;

    scalarList mpe(SubList<scalar>(massEscape_, nPatches));
    Pstream::listCombineGather(mpe, plusEqOp<scalar>());
    mpe = mpe + mpe0;

    labelList nps(SubList<label>(nStick_, nPatches));
    Pstream::listCombineGather(nps, plusEqOp<label>());
    nps = nps + nps0;

    scalarList mps(SubList<scalar>(massStick_, nPatches));
    Pstream::listCombineGather(mps, plusEqOp<scalar>());
    mps = mps + mps0;

    for (label patchi = 0; patchi < nPatches; ++ patchi)
    {
        if
        (
            patchInteractionTypes_[patchi]
         != PatchInteractionModel<CloudType>::itNone
        )
        {
            os  << "    Parcel fate (number, mass)      : patch "
                << this->owner().mesh().boundaryMesh()[patchi].name() << nl
                << "      - escape                      = " << npe[patchi]
                << ", " << mpe[patchi] << nl
                << "      - stick                       = " << nps[patchi]
                << ", " << mps[patchi] << nl;
        }
    }

    if (this->writeTime())
    {
        this->setModelProperty("nEscape", npe);
        nEscape_ = 0;

        this->setModelProperty("massEscape", mpe);
        massEscape_ = scalar(0);

        this->setModelProperty("nStick", nps);
        nStick_ = 0;

        this->setModelProperty("massStick", mps);
        massStick_ = scalar(0);
    }
}


// ************************************************************************* //
