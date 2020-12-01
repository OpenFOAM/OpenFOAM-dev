/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "MomentumParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::MomentumParcel<ParcelType>::propertyList_ =
    Foam::MomentumParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::MomentumParcel<ParcelType>::sizeofFields_
(
    sizeof(MomentumParcel<ParcelType>)
  - offsetof(MomentumParcel<ParcelType>, moving_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::MomentumParcel<ParcelType>::MomentumParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    moving_(false),
    typeId_(0),
    nParticle_(0.0),
    d_(0.0),
    dTarget_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            moving_ = readBool(is);
            typeId_ = readLabel(is);
            nParticle_ = readScalar(is);
            d_ = readScalar(is);
            dTarget_ = readScalar(is);
            is >> U_;
            rho_ = readScalar(is);
            age_ = readScalar(is);
            tTurb_ = readScalar(is);
            is >> UTurb_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&moving_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check
    (
        "MomentumParcel<ParcelType>::MomentumParcel"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::MomentumParcel<ParcelType>::readFields(CloudType& c)
{
    bool write = c.size();

    ParcelType::readFields(c);

    IOField<label> moving
    (
        c.fieldIOobject("active", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, moving);

    IOField<label> typeId
    (
        c.fieldIOobject("typeId", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, typeId);

    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> d
    (
        c.fieldIOobject("d", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, d);

    IOField<scalar> dTarget
    (
        c.fieldIOobject("dTarget", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, dTarget);

    IOField<vector> U
    (
        c.fieldIOobject("U", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, U);

    IOField<scalar> rho
    (
        c.fieldIOobject("rho", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, rho);

    IOField<scalar> age
    (
        c.fieldIOobject("age", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, age);

    IOField<scalar> tTurb
    (
        c.fieldIOobject("tTurb", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb
    (
        c.fieldIOobject("UTurb", IOobject::MUST_READ),
        write
    );
    c.checkFieldIOobject(c, UTurb);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        MomentumParcel<ParcelType>& p = iter();

        p.moving_ = moving[i];
        p.typeId_ = typeId[i];
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        p.dTarget_ = dTarget[i];
        p.U_ = U[i];
        p.rho_ = rho[i];
        p.age_ = age[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::MomentumParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<label> moving(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> age(c.fieldIOobject("age", IOobject::NO_READ), np);
    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::NO_READ), np);
    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const MomentumParcel<ParcelType>& p = iter();

        moving[i] = p.moving();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        dTarget[i] = p.dTarget();
        U[i] = p.U();
        rho[i] = p.rho();
        age[i] = p.age();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();

        i++;
    }

    const bool write = np > 0;

    moving.write(write);
    typeId.write(write);
    nParticle.write(write);
    d.write(write);
    dTarget.write(write);
    U.write(write);
    rho.write(write);
    age.write(write);
    tTurb.write(write);
    UTurb.write(write);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MomentumParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.moving()
            << token::SPACE << p.typeId()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.d()
            << token::SPACE << p.dTarget()
            << token::SPACE << p.U()
            << token::SPACE << p.rho()
            << token::SPACE << p.age()
            << token::SPACE << p.tTurb()
            << token::SPACE << p.UTurb();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.moving_),
            MomentumParcel<ParcelType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const MomentumParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
