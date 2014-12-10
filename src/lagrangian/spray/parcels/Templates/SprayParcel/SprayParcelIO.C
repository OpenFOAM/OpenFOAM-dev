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

#include "SprayParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::SprayParcel<ParcelType>::propHeader =
    ParcelType::propHeader
  + " d0"
  + " position0"
  + " sigma"
  + " mu"
  + " liquidCore"
  + " KHindex"
  + " y"
  + " yDot"
  + " tc"
  + " ms"
  + " injector"
  + " tMom"
  + " user";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SprayParcel<ParcelType>::SprayParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    d0_(0.0),
    position0_(vector::zero),
    sigma_(0.0),
    mu_(0.0),
    liquidCore_(0.0),
    KHindex_(0.0),
    y_(0.0),
    yDot_(0.0),
    tc_(0.0),
    ms_(0.0),
    injector_(1.0),
    tMom_(GREAT),
    user_(0.0)
{
    if (readFields)
    {

        if (is.format() == IOstream::ASCII)
        {
            d0_ = readScalar(is);
            is >> position0_;
            sigma_ = readScalar(is);
            mu_ = readScalar(is);
            liquidCore_ = readScalar(is);
            KHindex_ = readScalar(is);
            y_ = readScalar(is);
            yDot_ = readScalar(is);
            tc_ = readScalar(is);
            ms_ = readScalar(is);
            injector_ = readScalar(is);
            tMom_ = readScalar(is);
            user_ = readScalar(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&d0_),
                sizeof(d0_)
              + sizeof(position0_)
              + sizeof(sigma_)
              + sizeof(mu_)
              + sizeof(liquidCore_)
              + sizeof(KHindex_)
              + sizeof(y_)
              + sizeof(yDot_)
              + sizeof(tc_)
              + sizeof(ms_)
              + sizeof(injector_)
              + sizeof(tMom_)
              + sizeof(user_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "SprayParcel<ParcelType>::SprayParcel"
        "("
            "const polyMesh, "
            "Istream&, "
            "bool"
        ")"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::SprayParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::SprayParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c, compModel);

    IOField<scalar> d0(c.fieldIOobject("d0", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d0);

    IOField<vector> position0
    (
        c.fieldIOobject("position0", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, position0);

    IOField<scalar> sigma(c.fieldIOobject("sigma", IOobject::MUST_READ));
    c.checkFieldIOobject(c, sigma);

    IOField<scalar> mu(c.fieldIOobject("mu", IOobject::MUST_READ));
    c.checkFieldIOobject(c, mu);

    IOField<scalar> liquidCore(c.fieldIOobject
    (
        "liquidCore", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, liquidCore);

    IOField<scalar> KHindex(c.fieldIOobject("KHindex", IOobject::MUST_READ));
    c.checkFieldIOobject(c, KHindex);

    IOField<scalar> y(c.fieldIOobject("y", IOobject::MUST_READ));
    c.checkFieldIOobject(c, y);

    IOField<scalar> yDot(c.fieldIOobject("yDot", IOobject::MUST_READ));
    c.checkFieldIOobject(c, yDot);

    IOField<scalar> tc(c.fieldIOobject("tc", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tc);

    IOField<scalar> ms(c.fieldIOobject("ms", IOobject::MUST_READ));
    c.checkFieldIOobject(c, ms);

    IOField<scalar> injector(c.fieldIOobject("injector", IOobject::MUST_READ));
    c.checkFieldIOobject(c, injector);

    IOField<scalar> tMom(c.fieldIOobject("tMom", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tMom);

    IOField<scalar> user(c.fieldIOobject("user", IOobject::MUST_READ));
    c.checkFieldIOobject(c, user);

    label i = 0;
    forAllIter(typename Cloud<SprayParcel<ParcelType> >, c, iter)
    {
        SprayParcel<ParcelType>& p = iter();
        p.d0_ = d0[i];
        p.position0_ = position0[i];
        p.sigma_ = sigma[i];
        p.mu_ = mu[i];
        p.liquidCore_ = liquidCore[i];
        p.KHindex_ = KHindex[i];
        p.y_ = y[i];
        p.yDot_ = yDot[i];
        p.tc_ = tc[i];
        p.ms_ = ms[i];
        p.injector_ = injector[i];
        p.tMom_ = tMom[i];
        p.user_ = user[i];
        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::SprayParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::SprayParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c, compModel);

    label np = c.size();

    IOField<scalar> d0(c.fieldIOobject("d0", IOobject::NO_READ), np);
    IOField<vector> position0
    (
        c.fieldIOobject("position0", IOobject::NO_READ),
        np
    );
    IOField<scalar> sigma(c.fieldIOobject("sigma", IOobject::NO_READ), np);
    IOField<scalar> mu(c.fieldIOobject("mu", IOobject::NO_READ), np);
    IOField<scalar> liquidCore
    (
        c.fieldIOobject("liquidCore", IOobject::NO_READ),
        np
    );
    IOField<scalar> KHindex(c.fieldIOobject("KHindex", IOobject::NO_READ), np);
    IOField<scalar> y(c.fieldIOobject("y", IOobject::NO_READ), np);
    IOField<scalar> yDot(c.fieldIOobject("yDot", IOobject::NO_READ), np);
    IOField<scalar> tc(c.fieldIOobject("tc", IOobject::NO_READ), np);
    IOField<scalar> ms(c.fieldIOobject("ms", IOobject::NO_READ), np);
    IOField<scalar> injector
    (
        c.fieldIOobject("injector", IOobject::NO_READ),
        np
    );
    IOField<scalar> tMom(c.fieldIOobject("tMom", IOobject::NO_READ), np);
    IOField<scalar> user(c.fieldIOobject("user", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<SprayParcel<ParcelType> >, c, iter)
    {
        const SprayParcel<ParcelType>& p = iter();
        d0[i] = p.d0_;
        position0[i] = p.position0_;
        sigma[i] = p.sigma_;
        mu[i] = p.mu_;
        liquidCore[i] = p.liquidCore_;
        KHindex[i] = p.KHindex_;
        y[i] = p.y_;
        yDot[i] = p.yDot_;
        tc[i] = p.tc_;
        ms[i] = p.ms_;
        injector[i] = p.injector_;
        tMom[i] = p.tMom_;
        user[i] = p.user_;
        i++;
    }

    d0.write();
    position0.write();
    sigma.write();
    mu.write();
    liquidCore.write();
    KHindex.write();
    y.write();
    yDot.write();
    tc.write();
    ms.write();
    injector.write();
    tMom.write();
    user.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const SprayParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
        << token::SPACE << p.d0()
        << token::SPACE << p.position0()
        << token::SPACE << p.sigma()
        << token::SPACE << p.mu()
        << token::SPACE << p.liquidCore()
        << token::SPACE << p.KHindex()
        << token::SPACE << p.y()
        << token::SPACE << p.yDot()
        << token::SPACE << p.tc()
        << token::SPACE << p.ms()
        << token::SPACE << p.injector()
        << token::SPACE << p.tMom()
        << token::SPACE << p.user();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.d0_),
            sizeof(p.d0())
          + sizeof(p.position0())
          + sizeof(p.sigma())
          + sizeof(p.mu())
          + sizeof(p.liquidCore())
          + sizeof(p.KHindex())
          + sizeof(p.y())
          + sizeof(p.yDot())
          + sizeof(p.tc())
          + sizeof(p.ms())
          + sizeof(p.injector())
          + sizeof(p.tMom())
          + sizeof(p.user())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const SprayParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
