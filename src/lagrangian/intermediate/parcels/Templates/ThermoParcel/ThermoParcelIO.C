/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "ThermoParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ThermoParcel<ParcelType>::propertyList_ =
    Foam::ThermoParcel<ParcelType>::propertyList();


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    T_(0.0),
    Cp_(0.0),
    Tc_(0.0),
    Cpc_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            T_ = readScalar(is);
            Cp_ = readScalar(is);
        }
        else
        {
            label size = long(&Cp_) - long(&T_) + sizeof(Cp_);
            is.read(reinterpret_cast<char*>(&T_), size);
        }
    }

    // Check state of Istream
    is.check
    (
        "ThermoParcel::ThermoParcel(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::ThermoParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    IOField<scalar> T(c.fieldIOobject("T", IOobject::MUST_READ));
    c.checkFieldIOobject(c, T);

    IOField<scalar> Cp(c.fieldIOobject("Cp", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Cp);


    label i = 0;
    forAllIter(typename Cloud<ThermoParcel<ParcelType> >, c, iter)
    {
        ThermoParcel<ParcelType>& p = iter();

        p.T_ = T[i];
        p.Cp_ = Cp[i];
        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ThermoParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np);
    IOField<scalar> Cp(c.fieldIOobject("Cp", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<ThermoParcel<ParcelType> >, c, iter)
    {
        const ThermoParcel<ParcelType>& p = iter();

        T[i] = p.T_;
        Cp[i] = p.Cp_;
        i++;
    }

    T.write();
    Cp.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ThermoParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.T()
            << token::SPACE << p.Cp();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);

        label size = long(&p.Cp_) - long(&p.T_) + sizeof(p.Cp_);
        os.write(reinterpret_cast<const char*>(&p.T_), size);
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const ThermoParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
