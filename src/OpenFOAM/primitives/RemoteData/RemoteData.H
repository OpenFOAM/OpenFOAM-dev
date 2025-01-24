/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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
    Foam::RemoteData

Description
    Struct for keeping processor, element (cell, face, point) and a piece of
    data. Used for finding minimum values across multiple processes.

SourceFiles
    RemoteDataI.H

\*---------------------------------------------------------------------------*/

#ifndef RemoteData_H
#define RemoteData_H

#include "remote.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type>
class RemoteData;
template<class Type>
inline bool operator==(const RemoteData<Type>&, const RemoteData<Type>&);
template<class Type>
inline bool operator!=(const RemoteData<Type>&, const RemoteData<Type>&);
template<class Type>
inline Istream& operator>>(Istream&, RemoteData<Type>&);
template<class Type>
inline Ostream& operator<<(Ostream&, const RemoteData<Type>&);

/*---------------------------------------------------------------------------*\
                          Class RemoteData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class RemoteData
:
    public remote
{
public:

    //- Component type
    typedef RemoteData<Type> cmptType;


    // Public Data

        //- Data
        Type data;


    // Constructors

        //- Construct null
        inline RemoteData();

        //- Construct from components
        inline RemoteData(const label, const label, const Type&);

        //- Construct from element index and data. Pstream will be used to
        //  determine the processor index.
        inline RemoteData(const label, const Type&);

        //- Construct from stream
        inline RemoteData(Istream& is);


    // Public Classes

        //- Operator to take the first valid process
        struct firstProcOp
        {
            inline RemoteData<Type> operator()
            (
                const RemoteData<Type>& a,
                const RemoteData<Type>& b
            ) const;
        };

        //- Operator to take the first valid process
        struct firstProcEqOp
        {
            inline void operator()
            (
                RemoteData<Type>& a,
                const RemoteData<Type>& b
            ) const;
        };

        //- Operator to take smallest valid value
        struct smallestEqOp
        {
            inline void operator()
            (
                RemoteData<Type>& a,
                const RemoteData<Type>& b
            ) const;
        };

        //- Operator to take greatest valid value
        struct greatestEqOp
        {
            inline void operator()
            (
                RemoteData<Type>& a,
                const RemoteData<Type>& b
            ) const;
        };

        //- Operator to take smallest first valid value
        struct smallestFirstEqOp
        {
            inline void operator()
            (
                RemoteData<Type>& a,
                const RemoteData<Type>& b
            ) const;
        };

        //- Operator to take greatest first valid value
        struct greatestFirstEqOp
        {
            inline void operator()
            (
                RemoteData<Type>& a,
                const RemoteData<Type>& b
            ) const;
        };


    // Friend Operators

        //- Equality comparison
        friend bool operator== <Type>
        (
            const RemoteData<Type>& a,
            const RemoteData<Type>& b
        );

        //- Inequality comparison
        friend bool operator!= <Type>
        (
            const RemoteData<Type>& a,
            const RemoteData<Type>& b
        );


    // IOstream Operators

        //- Write to stream
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const RemoteData<Type>& p
        );

        //- Read from stream
        friend Istream& operator>> <Type>
        (
            Istream& is,
            RemoteData<Type>& p
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RemoteDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
