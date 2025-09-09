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

#include "carried.H"
#include "fvcCurl.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(carried, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#define ACCESS_CARRIER_FIELDS(Type, nullArg)                                   \
namespace Foam                                                                 \
{                                                                              \
    namespace clouds                                                           \
    {                                                                          \
        template<>                                                             \
        HashPtrTable<CarrierField<Type>>& carried::carrierFields() const       \
        {                                                                      \
            return CAT3(carrier, CAPITALIZE(Type), Fields_);                   \
        }                                                                      \
    }                                                                          \
}
FOR_ALL_FIELD_TYPES(ACCESS_CARRIER_FIELDS)
#undef ACCESS_CARRIER_FIELDS


Foam::autoPtr<Foam::volVectorField> Foam::clouds::carried::readDUdtc
(
    const cloud& c
) const
{
    typeIOobject<volVectorField> io
    (
        "ddt(" + Uc.psi().name() + ")",
        Uc.psi().mesh().time().name(),
        Uc.psi().mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    return
        autoPtr<volVectorField>
        (
            io.headerOk()
          ? new volVectorField(io, Uc.psi().mesh())
          : nullptr
        );
}


const Foam::volVectorField& Foam::clouds::carried::dUdtc() const
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

void Foam::clouds::carried::clearCarrierFields()
{
    #define CLEAR_TYPE_CARRIER_FIELDS(Type, nullArg)                           \
        forAllIter                                                             \
        (                                                                      \
            HashPtrTable<CarrierField<Type>>,                                  \
            carrierFields<Type>(),                                             \
            iter                                                               \
        )                                                                      \
        {                                                                      \
            iter()->clear(true);                                               \
        }
    FOR_ALL_FIELD_TYPES(CLEAR_TYPE_CARRIER_FIELDS);
    #undef CLEAR_TYPE_CARRIER_FIELDS
}


void Foam::clouds::carried::resetCarrierFields(const bool initial)
{
    #define RESET_TYPE_CARRIER_FIELDS(Type, nullArg)                           \
        forAllIter                                                             \
        (                                                                      \
            HashPtrTable<CarrierField<Type>>,                                  \
            carrierFields<Type>(),                                             \
            iter                                                               \
        )                                                                      \
        {                                                                      \
            iter()->reset(initial);                                            \
       }
    FOR_ALL_FIELD_TYPES(RESET_TYPE_CARRIER_FIELDS);
    #undef RESET_TYPE_CARRIER_FIELDS
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::carried::carried(const cloud& c, const dictionary& dict)
:
    carrierPhaseName_
    (
        dict.found("phase") || dict.found("carrierPhase")
      ? dict.lookup<word>("carrierPhase")
      : word::null
    ),
    phaseName_
    (
        dict.lookupOrDefault<word>("phase", word::null)
    ),
    Uc
    (
        carrierField<vector>
        (
            c.mesh().mesh().lookupObject<volVectorField>
            (
                c.mesh().mesh().foundObject<volVectorField>
                (
                    IOobject::groupName("U", carrierPhaseName())
                )
              ? IOobject::groupName("U", carrierPhaseName())
              : c.mesh().mesh().foundObject<volVectorField>
                (
                    "U"
                )
              ? "U"
              : IOobject::groupName("U", carrierPhaseName())
            )
        )
    ),
    curlUc
    (
        carrierField<vector>
        (
            IOobject::groupName("curlUc", carrierPhaseName()),
            [&]()
            {
                return fvc::curl(Uc.psi());
            }
        )
    ),
    dUdtcPtr_(readDUdtc(c)),
    DUDtc
    (
        carrierField<vector>
        (
            IOobject::groupName("DUDtc", carrierPhaseName()),
            [&]()
            {
                return dUdtc() + (Uc.psi() & fvc::grad(Uc.psi()));
            }
        )
    ),
    UcPhase
    (
        hasPhase()
      ? Uc.psi().group() == word::null
      ? Uc
      : carrierField<vector>
        (
            c.mesh().mesh().lookupObject<volVectorField>
            (
                IOobject::groupName("U", phaseName())
            )
        )
      : carrierField<vector>
        (
            "UcPhase",
            [&]()
            {
                FatalErrorInFunction
                    << "Cloud " << c.name() << " does not have a corresponding "
                    << "Eulerian phase velocity" << exit(FatalError);
                return tmp<volVectorField>(nullptr);
            }
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::carried::~carried()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::word& Foam::clouds::carried::carrierPhaseName() const
{
    return carrierPhaseName_;
}


const Foam::word& Foam::clouds::carried::phaseName() const
{
    return phaseName_;
}


bool Foam::clouds::carried::hasPhase() const
{
    return phaseName_ != word::null;
}


Foam::word Foam::clouds::carried::nameToCarrierName(const word& name)
{
    return
        IOobject::groupName
        (
            IOobject::member(name) + 'c',
            IOobject::group(name)
        );
}


Foam::word Foam::clouds::carried::carrierNameToName(const word& namec)
{
    const word memberc = IOobject::member(namec);

    if (memberc[memberc.size() - 1] != 'c')
    {
        FatalErrorInFunction
            << "Name " << namec << " is not a carrier name"
            << exit(FatalError);
    }

    return
        IOobject::groupName
        (
            memberc(memberc.size() - 1),
            IOobject::group(namec)
        );
}


// ************************************************************************* //
