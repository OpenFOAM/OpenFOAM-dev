/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

template<class Type>
void Foam::mappedPatchBase::distribute(List<Type>& lst) const
{
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            lst = AMI().interpolateToSource(Field<Type>(lst.xfer()));
            break;
        }
        default:
        {
            map().distribute(lst);
        }
    }
}


template<class Type, class CombineOp>
void Foam::mappedPatchBase::distribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            lst = AMI().interpolateToSource
                (
                    Field<Type>(lst.xfer()),
                    cop
                );
            break;
        }
        default:
        {
            mapDistributeBase::distribute
            (
                Pstream::defaultCommsType,
                map().schedule(),
                map().constructSize(),
                map().subMap(),
                false,
                map().constructMap(),
                false,
                lst,
                cop,
                flipOp(),
                Type(Zero)
            );
        }
    }
}


template<class Type>
void Foam::mappedPatchBase::reverseDistribute(List<Type>& lst) const
{
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            lst = AMI().interpolateToTarget(Field<Type>(lst.xfer()));
            break;
        }
        default:
        {
            map().reverseDistribute(sampleSize(), lst);
            break;
        }
    }
}


template<class Type, class CombineOp>
void Foam::mappedPatchBase::reverseDistribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            lst = AMI().interpolateToTarget
                (
                    Field<Type>(lst.xfer()),
                    cop
                );
            break;
        }
        default:
        {
            label cSize = sampleSize();
            mapDistributeBase::distribute
            (
                Pstream::defaultCommsType,
                map().schedule(),
                cSize,
                map().constructMap(),
                false,
                map().subMap(),
                false,
                lst,
                cop,
                flipOp(),
                Type(Zero)
            );
            break;
        }
    }
}


// ************************************************************************* //
