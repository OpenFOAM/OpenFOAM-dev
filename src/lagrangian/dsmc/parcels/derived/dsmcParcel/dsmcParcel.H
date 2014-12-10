/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::dsmcParcel

Description
    Declaration of dsmc parcel type

SourceFiles
    dsmcParcel.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcParcel_H
#define dsmcParcel_H

#include "particle.H"
#include "DsmcParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    typedef DsmcParcel<particle> dsmcParcel;
/*---------------------------------------------------------------------------*\
                         Class dsmcParcel Declaration
\*---------------------------------------------------------------------------*/
/*
class dsmcParcel
:
    public DsmcParcel<dsmcParcel>
{

public:

    //- Run-time type information
    TypeName("dsmcParcel");


    // Constructors

        //- Construct from components
        dsmcParcel
        (
            DsmcCloud<dsmcParcel>& owner,
            const vector& position,
            const vector& U,
            const scalar Ei,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeId
        );

        //- Construct from Istream
        dsmcParcel
        (
            const Cloud<dsmcParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<Particle<dsmcParcel> > clone() const
        {
            return autoPtr<Particle<dsmcParcel> >(new dsmcParcel(*this));
        }


    //- Destructor
    virtual ~dsmcParcel();
};


template<>
inline bool contiguous<dsmcParcel>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
*/
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
