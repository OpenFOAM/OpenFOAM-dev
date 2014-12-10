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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::cyclicAMIPolyPatch::interpolate
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    if (owner())
    {
        return AMI().interpolateToSource(fld, defaultValues);
    }
    else
    {
        return neighbPatch().AMI().interpolateToTarget(fld, defaultValues);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::cyclicAMIPolyPatch::interpolate
(
    const tmp<Field<Type> >& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), defaultValues);
}


template<class Type, class CombineOp>
void Foam::cyclicAMIPolyPatch::interpolate
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    if (owner())
    {
        AMI().interpolateToSource
        (
            fld,
            cop,
            result,
            defaultValues
        );
    }
    else
    {
        neighbPatch().AMI().interpolateToTarget
        (
            fld,
            cop,
            result,
            defaultValues
        );
    }
}


// ************************************************************************* //
