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

template<class SourcePatch, class TargetPatch>
inline Foam::label
Foam::AMIInterpolation<SourcePatch, TargetPatch>::singlePatchProc() const
{
    return singlePatchProc_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalar
Foam::AMIInterpolation<SourcePatch, TargetPatch>::lowWeightCorrection() const
{
    return lowWeightCorrection_;
}


template<class SourcePatch, class TargetPatch>
inline bool
Foam::AMIInterpolation<SourcePatch, TargetPatch>::
applyLowWeightCorrection() const
{
    return lowWeightCorrection_ > 0;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcMagSf() const
{
    return srcMagSf_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::labelListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcAddress() const
{
    return srcAddress_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcWeights() const
{
    return srcWeights_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcWeights()
{
    return srcWeights_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcWeightsSum() const
{
    return srcWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcWeightsSum()
{
    return srcWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::mapDistribute&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcMap() const
{
    return srcMapPtr_();
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtMagSf() const
{
    return tgtMagSf_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::labelListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtAddress() const
{
    return tgtAddress_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtWeights() const
{
    return tgtWeights_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtWeights()
{
    return tgtWeights_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtWeightsSum() const
{
    return tgtWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtWeightsSum()
{
    return tgtWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::mapDistribute&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtMap() const
{
    return tgtMapPtr_();
}


// ************************************************************************* //
