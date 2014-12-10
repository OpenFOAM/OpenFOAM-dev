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

\*---------------------------------------------------------------------------*/

#include "openFoamTableReader.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::openFoamTableReader<Type>::openFoamTableReader(const dictionary& dict)
:
    tableReader<Type>(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::openFoamTableReader<Type>::~openFoamTableReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::openFoamTableReader<Type>::operator()
(
    const fileName& fName,
    List<Tuple2<scalar, Type> >& data
)
{
    // Read data from file
    IFstream(fName)() >> data;
}


template<class Type>
void Foam::openFoamTableReader<Type>::operator()
(
    const fileName& fName,
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >& data
)
{
    // Read data from file
    IFstream(fName)() >> data;
}


// ************************************************************************* //
