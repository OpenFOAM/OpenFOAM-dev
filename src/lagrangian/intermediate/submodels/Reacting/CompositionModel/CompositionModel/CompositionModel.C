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

#include "CompositionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    thermo_(owner.thermo()),
    phaseProps_()
{}


template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    thermo_(owner.thermo()),
    phaseProps_
    (
        this->coeffDict().lookup("phases"),
        thermo_.carrier().species(),
        thermo_.liquids().components(),
        thermo_.solids().components()
    )
{}


template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel
(
    const CompositionModel<CloudType>& cm
)
:
    CloudSubModelBase<CloudType>(cm),
    thermo_(cm.thermo_),
    phaseProps_(cm.phaseProps_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CompositionModel<CloudType>::~CompositionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const Foam::SLGThermo& Foam::CompositionModel<CloudType>::thermo() const
{
    return thermo_;
}


template<class CloudType>
const Foam::basicMultiComponentMixture&
Foam::CompositionModel<CloudType>::carrier() const
{
    return thermo_.carrier();
}


template<class CloudType>
const Foam::liquidMixtureProperties&
Foam::CompositionModel<CloudType>::liquids() const
{
    return thermo_.liquids();
}


template<class CloudType>
const Foam::solidMixtureProperties&
Foam::CompositionModel<CloudType>::solids() const
{
    return thermo_.solids();
}


template<class CloudType>
const Foam::phasePropertiesList&
Foam::CompositionModel<CloudType>::phaseProps() const
{
    return phaseProps_;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::nPhase() const
{
    return phaseProps_.size();
}


template<class CloudType>
const Foam::wordList& Foam::CompositionModel<CloudType>::phaseTypes() const
{
    // if only 1 phase, return the constituent component names
    if (phaseProps_.size() == 1)
    {
        return phaseProps_[0].names();
    }
    else
    {
        return phaseProps_.phaseTypes();
    }
}


template<class CloudType>
const Foam::wordList& Foam::CompositionModel<CloudType>::stateLabels() const
{
    return phaseProps_.stateLabels();
}


template<class CloudType>
const Foam::wordList&
Foam::CompositionModel<CloudType>::componentNames(const label phaseI) const
{
    return phaseProps_[phaseI].names();
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::globalCarrierId
(
    const word& cmptName,
    const bool allowNotFound
) const
{
    label id = thermo_.carrierId(cmptName);

    if (id < 0 && !allowNotFound)
    {
        FatalErrorIn
        (
            "Foam::label Foam::CompositionModel<CloudType>::globalCarrierId"
            "("
                "const word&, "
                "const bool"
            ") const"
        )   << "Unable to determine global id for requested component "
            << cmptName << ". Available components are " << nl
            << thermo_.carrier().species() << abort(FatalError);
    }

    return id;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::globalId
(
    const label phaseI,
    const word& cmptName,
    const bool allowNotFound
) const
{
    label id = phaseProps_[phaseI].globalId(cmptName);

    if (id < 0 && !allowNotFound)
    {
        FatalErrorIn
        (
            "Foam::label Foam::CompositionModel<CloudType>::globalId"
            "("
                "const label, "
                "const word&, "
                "const bool"
            ") const"
        )   << "Unable to determine global id for requested component "
            << cmptName << abort(FatalError);
    }

    return id;
}


template<class CloudType>
const Foam::labelList& Foam::CompositionModel<CloudType>::globalIds
(
    const label phaseI
) const
{
    return phaseProps_[phaseI].globalIds();
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::localId
(
    const label phaseI,
    const word& cmptName,
    const bool allowNotFound
) const
{
    label id = phaseProps_[phaseI].id(cmptName);

    if (id < 0 && !allowNotFound)
    {
        FatalErrorIn
        (
            "Foam::label Foam::CompositionModel<CloudType>::localId"
            "("
                "const label, "
                "const word&, "
                "const bool"
            ") const"
        )   << "Unable to determine local id for component " << cmptName
            << abort(FatalError);
    }

    return id;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::localToGlobalCarrierId
(
    const label phaseI,
    const label id,
    const bool allowNotFound
) const
{
    label gid = phaseProps_[phaseI].globalCarrierIds()[id];

    if (gid < 0 && !allowNotFound)
    {
        FatalErrorIn
        (
            "Foam::label "
            "Foam::CompositionModel<CloudType>::localToGlobalCarrierId"
            "("
                "const label, "
                "const label, "
                "const bool"
            ") const"
        )   << "Unable to determine global carrier id for phase "
            << phaseI << " with local id " << id
            << abort(FatalError);
    }

    return gid;
}


template<class CloudType>
const Foam::scalarField& Foam::CompositionModel<CloudType>::Y0
(
    const label phaseI
) const
{
    return phaseProps_[phaseI].Y();
}


template<class CloudType>
Foam::scalarField Foam::CompositionModel<CloudType>::X
(
    const label phaseI,
    const scalarField& Y
) const
{
    const phaseProperties& props = phaseProps_[phaseI];
    scalarField X(Y.size(), 0.0);
    scalar WInv = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                WInv += Y[i]/thermo_.carrier().W(gid);
                X[i] = Y[i]/thermo_.carrier().W(gid);
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                WInv += Y[i]/thermo_.liquids().properties()[gid].W();
                X[i] += Y[i]/thermo_.liquids().properties()[gid].W();
            }
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalarField Foam::CompositionModel<CloudType>::X"
                "("
                    "const label, "
                    "const scalarField&"
                ") const"
            )   << "Only possible to convert gas and liquid mass fractions"
                << abort(FatalError);
        }
    }

    tmp<scalarField> tfld = X/WInv;
    return tfld();
}


template<class CloudType>
const Foam::scalarField& Foam::CompositionModel<CloudType>::YMixture0() const
{
    notImplemented
    (
        "const scalarField& Foam::CompositionModel<CloudType>::YMixture0() "
        "const"
    );

    return scalarField::null();
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::idGas() const
{
    notImplemented
    (
        "Foam::label Foam::CompositionModel<CloudType>::idGas() const"
    );

    return -1;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::idLiquid() const
{
    notImplemented
    (
        "Foam::label Foam::CompositionModel<CloudType>::idLiquid() const"
    );

    return -1;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::idSolid() const
{
    notImplemented
    (
        "Foam::label Foam::CompositionModel<CloudType>::idSolid() const"
    );

    return -1;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::H
(
    const label phaseI,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phaseI];
    scalar HMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HMixture += Y[i]*thermo_.carrier().Ha(gid, p, T);
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HMixture += Y[i]*thermo_.liquids().properties()[gid].h(p, T);
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HMixture +=
                     Y[i]
                    *(
                        thermo_.solids().properties()[gid].Hf()
                      + thermo_.solids().properties()[gid].Cp()*T
                     );
            }
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::CompositionModel<CloudType>::H"
                "("
                "    const label, "
                "    const scalarField&, "
                "    const scalar, "
                "    const scalar"
                ") const"
            )   << "Unknown phase enumeration" << abort(FatalError);
        }
    }

    return HMixture;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::Hs
(
    const label phaseI,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phaseI];
    scalar HsMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HsMixture += Y[i]*thermo_.carrier().Hs(gid, p, T);
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HsMixture +=
                    Y[i]
                   *(
                       thermo_.liquids().properties()[gid].h(p, T)
                     - thermo_.liquids().properties()[gid].h(p, 298.15)
                    );
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HsMixture += Y[i]*thermo_.solids().properties()[gid].Cp()*T;
            }
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::CompositionModel<CloudType>::Hs"
                "("
                "    const label, "
                "    const scalarField&, "
                "    const scalar, "
                "    const scalar"
                ") const"
            )   << "Unknown phase enumeration" << abort(FatalError);
        }
    }

    return HsMixture;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::Hc
(
    const label phaseI,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phaseI];
    scalar HcMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HcMixture += Y[i]*thermo_.carrier().Hc(gid);
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HcMixture +=
                    Y[i]*thermo_.liquids().properties()[gid].h(p, 298.15);
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HcMixture += Y[i]*thermo_.solids().properties()[gid].Hf();
            }
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::CompositionModel<CloudType>::Hc"
                "("
                "    const label, "
                "    const scalarField&, "
                "    const scalar, "
                "    const scalar"
                ") const"
            )   << "Unknown phase enumeration" << abort(FatalError);
        }
    }

    return HcMixture;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::Cp
(
    const label phaseI,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phaseI];
    scalar CpMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                CpMixture += Y[i]*thermo_.carrier().Cp(gid, p, T);
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                CpMixture += Y[i]*thermo_.liquids().properties()[gid].Cp(p, T);
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                CpMixture += Y[i]*thermo_.solids().properties()[gid].Cp();
            }
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::CompositionModel<CloudType>::Cp"
                "("
                    "const label, "
                    "const scalarField&, "
                    "const scalar, "
                    "const scalar"
                ") const"
            )   << "Unknown phase enumeration" << abort(FatalError);
        }
    }

    return CpMixture;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::L
(
    const label phaseI,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phaseI];
    scalar LMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            if (debug)
            {
                WarningIn
                (
                    "Foam::scalar Foam::CompositionModel<CloudType>::L"
                    "("
                        "const label, "
                        "const scalarField&, "
                        "const scalar, "
                        "const scalar"
                    ") const\n"
                )   << "No support for gaseous components" << endl;
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                LMixture += Y[i]*thermo_.liquids().properties()[gid].hl(p, T);
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            if (debug)
            {
                WarningIn
                (
                    "Foam::scalar Foam::CompositionModel<CloudType>::L"
                    "("
                        "const label, "
                        "const scalarField&, "
                        "const scalar, "
                        "const scalar"
                    ") const\n"
                )   << "No support for solid components" << endl;
            }
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::CompositionModel<CloudType>::L"
                "("
                    "const label, "
                    "const scalarField&, "
                    "const scalar, "
                    "const scalar"
                ") const"
            )   << "Unknown phase enumeration" << abort(FatalError);
        }
    }

    return LMixture;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CompositionModelNew.C"

// ************************************************************************* //

