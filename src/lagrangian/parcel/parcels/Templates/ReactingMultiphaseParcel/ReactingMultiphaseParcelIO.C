/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "ReactingMultiphaseParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ReactingMultiphaseParcel<ParcelType>::propertyList_ =
    Foam::ReactingMultiphaseParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::ReactingMultiphaseParcel<ParcelType>::sizeofFields_
(
    sizeof(scalar)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    Istream& is,
    bool readFields
)
:
    ParcelType(is, readFields),
    mass0_(0.0),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0)
{
    if (readFields)
    {
        DynamicList<scalar> Yg;
        DynamicList<scalar> Yl;
        DynamicList<scalar> Ys;

        if (is.format() == IOstream::ASCII)
        {
            is >> mass0_ >> Yg >> Yl >> Ys;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&mass0_), sizeofFields_);
            is >> Yg >> Yl >> Ys;
        }

        YGas_.transfer(Yg);
        YLiquid_.transfer(Yl);
        YSolid_.transfer(Ys);
    }

    // Check state of Istream
    is.check
    (
        "ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel"
        "("
            "const polyMesh&, "
            "Istream&, "
            "bool"
        ")"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::readFields(CloudType& c)
{
    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingMultiphaseParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    bool valid = c.size();

    ParcelType::readFields(c, compModel);

    IOField<scalar> mass0
    (
        c.fieldIOobject("mass0", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, mass0);

    label i = 0;
    forAllIter(typename CloudType, c, iter)
    {
        ReactingMultiphaseParcel<ParcelType>& p = iter();
        p.mass0_ = mass0[i++];
    }

    // Get names and sizes for each Y...
    const label idGas = compModel.idGas();
    const wordList& gasNames = compModel.componentNames(idGas);
    const label idLiquid = compModel.idLiquid();
    const wordList& liquidNames = compModel.componentNames(idLiquid);
    const label idSolid = compModel.idSolid();
    const wordList& solidNames = compModel.componentNames(idSolid);
    const wordList& stateLabels = compModel.stateLabels();

    // Set storage for each Y... for each parcel
    forAllIter(typename CloudType, c, iter)
    {
        ReactingMultiphaseParcel<ParcelType>& p = iter();
        p.YGas_.setSize(gasNames.size(), 0.0);
        p.YLiquid_.setSize(liquidNames.size(), 0.0);
        p.YSolid_.setSize(solidNames.size(), 0.0);
    }

    // Populate YGas for each parcel
    forAll(gasNames, j)
    {
        IOField<scalar> YGas
        (
            c.fieldIOobject
            (
                "Y" + gasNames[j] + stateLabels[idGas],
                IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        forAllIter(typename CloudType, c, iter)
        {
            ReactingMultiphaseParcel<ParcelType>& p = iter();
            p.YGas_[j] = YGas[i++]/(p.Y()[idGas] + rootVSmall);
        }
    }
    // Populate YLiquid for each parcel
    forAll(liquidNames, j)
    {
        IOField<scalar> YLiquid
        (
            c.fieldIOobject
            (
                "Y" + liquidNames[j] + stateLabels[idLiquid],
                 IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        forAllIter(typename CloudType, c, iter)
        {
            ReactingMultiphaseParcel<ParcelType>& p = iter();
            p.YLiquid_[j] = YLiquid[i++]/(p.Y()[idLiquid] + rootVSmall);
        }
    }
    // Populate YSolid for each parcel
    forAll(solidNames, j)
    {
        IOField<scalar> YSolid
        (
            c.fieldIOobject
            (
                "Y" + solidNames[j] + stateLabels[idSolid],
                IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        forAllIter(typename CloudType, c, iter)
        {
            ReactingMultiphaseParcel<ParcelType>& p = iter();
            p.YSolid_[j] = YSolid[i++]/(p.Y()[idSolid] + rootVSmall);
        }
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingMultiphaseParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c, compModel);

    label np = c.size();

    {
        IOField<scalar> mass0(c.fieldIOobject("mass0", IOobject::NO_READ), np);

        label i = 0;
        forAllConstIter(typename CloudType, c, iter)
        {
            const ReactingMultiphaseParcel<ParcelType>& p = iter();
            mass0[i++] = p.mass0_;
        }
        mass0.write(np > 0);

        // Write the composition fractions
        const wordList& stateLabels = compModel.stateLabels();

        const label idGas = compModel.idGas();
        const wordList& gasNames = compModel.componentNames(idGas);
        forAll(gasNames, j)
        {
            IOField<scalar> YGas
            (
                c.fieldIOobject
                (
                    "Y" + gasNames[j] + stateLabels[idGas],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter(typename CloudType, c, iter)
            {
                const ReactingMultiphaseParcel<ParcelType>& p0 = iter();
                YGas[i++] = p0.YGas()[j]*p0.Y()[idGas];
            }

            YGas.write(np > 0);
        }

        const label idLiquid = compModel.idLiquid();
        const wordList& liquidNames = compModel.componentNames(idLiquid);
        forAll(liquidNames, j)
        {
            IOField<scalar> YLiquid
            (
                c.fieldIOobject
                (
                    "Y" + liquidNames[j] + stateLabels[idLiquid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter(typename CloudType, c, iter)
            {
                const ReactingMultiphaseParcel<ParcelType>& p0 = iter();
                YLiquid[i++] = p0.YLiquid()[j]*p0.Y()[idLiquid];
            }

            YLiquid.write(np > 0);
        }

        const label idSolid = compModel.idSolid();
        const wordList& solidNames = compModel.componentNames(idSolid);
        forAll(solidNames, j)
        {
            IOField<scalar> YSolid
            (
                c.fieldIOobject
                (
                    "Y" + solidNames[j] + stateLabels[idSolid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter(typename CloudType, c, iter)
            {
                const ReactingMultiphaseParcel<ParcelType>& p0 = iter();
                YSolid[i++] = p0.YSolid()[j]*p0.Y()[idSolid];
            }

            YSolid.write(np > 0);
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingMultiphaseParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.mass0()
            << token::SPACE << p.YGas()
            << token::SPACE << p.YLiquid()
            << token::SPACE << p.YSolid();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.mass0_),
            ReactingMultiphaseParcel<ParcelType>::sizeofFields_
        );
        os  << p.YGas() << p.YLiquid() << p.YSolid();
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const ReactingMultiphaseParcel<ParcelType>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
