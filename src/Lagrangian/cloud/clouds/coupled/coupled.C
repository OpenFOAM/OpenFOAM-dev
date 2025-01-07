/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "autoPtr.H"
#include "coupled.H"
#include "fvcCurl.H"
#include "fvcDdt.H"
#include "LagrangianmDdt.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(coupled, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#define ACCESS_CARRIER_FIELDS(Type, nullArg)                                   \
namespace Foam                                                                 \
{                                                                              \
    namespace clouds                                                           \
    {                                                                          \
        template<>                                                             \
        PtrDictionary<CarrierField<Type>>& coupled::carrierFields() const      \
        {                                                                      \
            return CAT3(carrier, CAPITALIZE(Type), Fields_);                   \
        }                                                                      \
    }                                                                          \
}
FOR_ALL_FIELD_TYPES(ACCESS_CARRIER_FIELDS)
#undef ACCESS_CARRIER_FIELDS


void Foam::clouds::coupled::clearCarrierFields()
{
    #define CLEAR_TYPE_CARRIER_FIELDS(Type, nullArg)                           \
        forAllIter                                                             \
        (                                                                      \
            PtrDictionary<CarrierField<Type>>,                                 \
            carrierFields<Type>(),                                             \
            iter                                                               \
        )                                                                      \
        {                                                                      \
            iter().clear(true);                                                \
       }
    FOR_ALL_FIELD_TYPES(CLEAR_TYPE_CARRIER_FIELDS);
    #undef CLEAR_TYPE_CARRIER_FIELDS
}


void Foam::clouds::coupled::resetCarrierFields(const bool predict)
{
    #define RESET_TYPE_CARRIER_FIELDS(Type, nullArg)                           \
        forAllIter                                                             \
        (                                                                      \
            PtrDictionary<CarrierField<Type>>,                                 \
            carrierFields<Type>(),                                             \
            iter                                                               \
        )                                                                      \
        {                                                                      \
            iter().reset(predict);                                             \
       }
    FOR_ALL_FIELD_TYPES(RESET_TYPE_CARRIER_FIELDS);
    #undef RESET_TYPE_CARRIER_FIELDS
}


#define ACCESS_CARRIER_EQNS(Type, nullArg)                                     \
namespace Foam                                                                 \
{                                                                              \
    namespace clouds                                                           \
    {                                                                          \
        template<>                                                             \
        PtrDictionary<CarrierEqn<Type>>& coupled::carrierEqns() const          \
        {                                                                      \
            return CAT3(carrier, CAPITALIZE(Type), Eqns_);                     \
        }                                                                      \
    }                                                                          \
}
FOR_ALL_FIELD_TYPES(ACCESS_CARRIER_EQNS)
#undef ACCESS_CARRIER_EQNS


Foam::autoPtr<Foam::volVectorField> Foam::clouds::coupled::readDUdtc
(
    const cloud& c
) const
{
    const volVectorField& U =
        c.mesh().mesh().lookupObject<volVectorField>("U");

    typeIOobject<volVectorField> io
    (
        "ddt(" + U.name() + ")",
        U.mesh().time().name(),
        U.mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    return
        autoPtr<volVectorField>
        (
            io.headerOk()
          ? new volVectorField(io, U.mesh())
          : nullptr
        );
}


const Foam::volVectorField& Foam::clouds::coupled::dUdtc() const
{
    if (Uc.psi().hasStoredOldTimes())
    {
        dUdtcPtr_.reset(fvc::ddt(Uc.psi()).ptr());
    }
    else if (!dUdtcPtr_.valid())
    {
        dUdtcPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "ddt(" + Uc.psi().name() + ")",
                    Uc.psi().mesh().time().name(),
                    Uc.psi().mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                Uc.psi().mesh(),
                dimensionedVector(dimVelocity/dimTime, Zero)
            )
        );
    }

    return dUdtcPtr_();
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::clouds::coupled::clearCarrierEqns()
{
    #define CLEAR_TYPE_CARRIER_EQNS(Type, nullArg)                             \
        forAllIter                                                             \
        (                                                                      \
            PtrDictionary<CarrierEqn<Type>>,                                   \
            carrierEqns<Type>(),                                               \
            iter                                                               \
        )                                                                      \
        {                                                                      \
            iter().clear();                                                    \
        }
    FOR_ALL_FIELD_TYPES(CLEAR_TYPE_CARRIER_EQNS);
    #undef CLEAR_TYPE_CARRIER_EQNS
}


void Foam::clouds::coupled::initialise(const bool predict)
{
    resetCarrierFields(predict);
}


void Foam::clouds::coupled::partition()
{
    clearCarrierFields();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupled::coupled(const cloud& c)
:
    dUdtcPtr_(readDUdtc(c)),
    Uc
    (
        carrierField<vector>
        (
            c.mesh().mesh().lookupObject<volVectorField>("U")
        )
    ),
    curlUc
    (
        carrierField<vector>
        (
            "curlUc",
            [&]()
            {
                return fvc::curl(Uc.psi());
            }
        )
    ),
    DUDtc
    (
        carrierField<vector>
        (
            "DUDtc",
            [&]()
            {
                return dUdtc() + (Uc.psi() & fvc::grad(Uc.psi()));
            }
        )
    ),
    nuc(c.derivedField<scalar>(*this, &coupled::calcNuc))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupled::~coupled()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::clouds::coupled::carrierName(const word& name)
{
    return
        IOobject::groupName
        (
            IOobject::member(name) + 'c',
            IOobject::group(name)
        );
}


// ************************************************************************* //
