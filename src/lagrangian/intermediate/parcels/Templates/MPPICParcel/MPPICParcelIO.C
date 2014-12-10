/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "MPPICParcel.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::MPPICParcel<ParcelType>::propertyList_ =
    Foam::MPPICParcel<ParcelType>::propertyList();


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    UCorrect_(vector::zero)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> UCorrect_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&UCorrect_),
              + sizeof(UCorrect_)
            );
        }
    }

    is.check
    (
        "MPPICParcel<ParcelType>::Collisions"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::MPPICParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    IOField<vector> UCorrect(c.fieldIOobject("UCorrect", IOobject::MUST_READ));
    c.checkFieldIOobject(c, UCorrect);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        MPPICParcel<ParcelType>& p = iter();

        p.UCorrect_ = UCorrect[i];

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::MPPICParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np =  c.size();

    IOField<vector>
        UCorrect(c.fieldIOobject("UCorrect", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const MPPICParcel<ParcelType>& p = iter();

        UCorrect[i] = p.UCorrect();

        i++;
    }

    UCorrect.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MPPICParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.UCorrect();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.UCorrect_),
          + sizeof(p.UCorrect())
        );
    }

    os.check
    (
        "Ostream& operator<<(Ostream&, const MPPICParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
