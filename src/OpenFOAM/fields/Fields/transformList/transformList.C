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

#include "transformList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T> Foam::transform
(
    const tensor& rotTensor,
    const UList<T>& field
)
{
    List<T> newField(field.size());

    forAll(field, i)
    {
        newField[i] = transform(rotTensor, field[i]);
    }

    return newField;
}


template<class T>
void Foam::transformList(const tensor& rotTensor, UList<T>& field)
{
    forAll(field, i)
    {
        field[i] = transform(rotTensor, field[i]);
    }
}


template<class T>
void Foam::transformList(const tensorField& rotTensor, UList<T>& field)
{
    if (rotTensor.size() == 1)
    {
        forAll(field, i)
        {
            field[i] = transform(rotTensor[0], field[i]);
        }
    }
    else if (rotTensor.size() == field.size())
    {
        forAll(field, i)
        {
            field[i] = transform(rotTensor[i], field[i]);
        }
    }
    else
    {
        FatalErrorIn
        (
            "transformList(const tensorField&, UList<T>&)"
        )   << "Sizes of field and transformation not equal. field:"
            << field.size() << " transformation:" << rotTensor.size()
            << abort(FatalError);
    }
}


template<class T>
void Foam::transformList(const tensor& rotTensor, Map<T>& field)
{
    forAllIter(typename Map<T>, field, iter)
    {
        iter() = transform(rotTensor[0], iter());
    }
}


template<class T>
void Foam::transformList(const tensorField& rotTensor, Map<T>& field)
{
    if (rotTensor.size() == 1)
    {
        forAllIter(typename Map<T>, field, iter)
        {
            iter() = transform(rotTensor[0], iter());
        }
    }
    else
    {
        FatalErrorIn
        (
            "transformList(const tensorField&, Map<T>&)"
        )   << "Multiple transformation tensors not supported. field:"
            << field.size() << " transformation:" << rotTensor.size()
            << abort(FatalError);
    }
}


template<class T>
void Foam::transformList(const tensor& rotTensor, EdgeMap<T>& field)
{
    forAllIter(typename EdgeMap<T>, field, iter)
    {
        iter() = transform(rotTensor[0], iter());
    }
}


template<class T>
void Foam::transformList(const tensorField& rotTensor, EdgeMap<T>& field)
{
    if (rotTensor.size() == 1)
    {
        forAllIter(typename EdgeMap<T>, field, iter)
        {
            iter() = transform(rotTensor[0], iter());
        }
    }
    else
    {
        FatalErrorIn
        (
            "transformList(const tensorField&, EdgeMap<T>&)"
        )   << "Multiple transformation tensors not supported. field:"
            << field.size() << " transformation:" << rotTensor.size()
            << abort(FatalError);
    }
}


// ************************************************************************* //
