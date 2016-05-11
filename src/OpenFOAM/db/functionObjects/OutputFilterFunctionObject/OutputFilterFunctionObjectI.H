/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class OutputFilter>
inline const Foam::Time&
Foam::OutputFilterFunctionObject<OutputFilter>::time() const
{
    return time_;
}


template<class OutputFilter>
inline const Foam::dictionary&
Foam::OutputFilterFunctionObject<OutputFilter>::dict() const
{
    return dict_;
}


template<class OutputFilter>
inline const Foam::word&
Foam::OutputFilterFunctionObject<OutputFilter>::regionName() const
{
    return regionName_;
}


template<class OutputFilter>
inline const Foam::outputFilterOutputControl&
Foam::OutputFilterFunctionObject<OutputFilter>::outputControl() const
{
    return outputControl_;
}


template<class OutputFilter>
inline const OutputFilter&
Foam::OutputFilterFunctionObject<OutputFilter>::outputFilter() const
{
    return filter_;
}


// ************************************************************************* //
