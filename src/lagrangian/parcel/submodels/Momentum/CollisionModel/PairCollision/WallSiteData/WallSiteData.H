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

Class
    Foam::WallSiteData

Description
    Stores the patch ID and templated data to represent a collision
    with a wall to be passed to the wall model.

SourceFiles
    WallSiteDataI.H
    WallSiteData.C
    WallSiteDataIO.C

\*---------------------------------------------------------------------------*/

#ifndef WallSiteData_H
#define WallSiteData_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class WallSiteData;

template<class Type>
Istream& operator>>(Istream&, WallSiteData<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const WallSiteData<Type>&);


/*---------------------------------------------------------------------------*\
                        Class WallSiteData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WallSiteData
{
    // Private Data

        //- Index of originating patch
        label patchi_;

        //- Wall data
        Type wallData_;


public:

    // Constructors

        //- Construct null
        WallSiteData();

        //- Construct from components
        WallSiteData
        (
            label patchi,
            const Type& wallData
        );


    //- Destructor
    ~WallSiteData();


    // Member Functions

            //- Return access to the patch index
            inline label patchIndex() const;

            //- Return non-const access to the patch index
            inline label& patchIndex();

            //- Return access to wall data
            inline const Type& wallData() const;

            //- Return non-const access to wall data
            inline Type& wallData();


    // Member Operators

        bool operator==(const WallSiteData<Type>&) const;
        bool operator!=(const WallSiteData<Type>&) const;


    // IOstream Operators

        friend Istream& operator>> <Type>
        (Istream&, WallSiteData<Type>&);

        friend Ostream& operator<< <Type>
        (Ostream&, const WallSiteData<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "WallSiteDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "WallSiteData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
