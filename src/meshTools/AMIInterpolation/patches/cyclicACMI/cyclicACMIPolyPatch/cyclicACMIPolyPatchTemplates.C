/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
Foam::tmp<Foam::Field<Type> > Foam::cyclicACMIPolyPatch::interpolate
(
    const Field<Type>& fldCouple,
    const Field<Type>& fldNonOverlap
) const
{
    // note: do not scale AMI field as face areas have already been taken
    // into account

    if (owner())
    {
        const scalarField& w = srcMask_;

        tmp<Field<Type> > interpField(AMI().interpolateToSource(fldCouple));

        return interpField + (1.0 - w)*fldNonOverlap;
    }
    else
    {
        const scalarField& w = neighbPatch().tgtMask();

        tmp<Field<Type> > interpField
        (
            neighbPatch().AMI().interpolateToTarget(fldCouple)
        );

        return interpField + (1.0 - w)*fldNonOverlap;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::cyclicACMIPolyPatch::interpolate
(
    const tmp<Field<Type> >& tFldCouple,
    const tmp<Field<Type> >& tFldNonOverlap
) const
{
    return interpolate(tFldCouple(), tFldNonOverlap());
}


template<class Type, class CombineOp>
void Foam::cyclicACMIPolyPatch::interpolate
(
    const UList<Type>& fldCouple,
    const UList<Type>& fldNonOverlap,
    const CombineOp& cop,
    List<Type>& result
) const
{
    // note: do not scale AMI field as face areas have already been taken
    // into account

    if (owner())
    {
        const scalarField& w = srcMask_;

        AMI().interpolateToSource(fldCouple, cop, result);
        result = result + (1.0 - w)*fldNonOverlap;
    }
    else
    {
        const scalarField& w = neighbPatch().tgtMask();

        neighbPatch().AMI().interpolateToTarget(fldCouple, cop, result);
        result = result + (1.0 - w)*fldNonOverlap;
    }
}


// ************************************************************************* //
