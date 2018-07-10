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

#ifndef steadyParticleTracksTemplates_H
#define steadyParticleTracksTemplates_H

#include "OFstream.H"
#include "SortableList.H"
#include "IOobjectList.H"
#include "PtrList.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool fieldOk(const IOobjectList& cloudObjs, const word& name);

template<class Type>
tmp<Field<Type>> readParticleField
(
    const word& name,
    const IOobjectList cloudObjs
);

template<class Type>
void readFields
(
    PtrList<List<Type>>& values,
    const List<word>& fields,
    const IOobjectList& cloudObjs
);

template<class Type>
void writeVTK(OFstream& os, const Type& value);

template<class Type>
void writeVTKFields
(
    OFstream& os,
    const PtrList<List<Type>>& values,
    const List<List<label>>& addr,
    const List<word>& fieldNames
);

template<class Type>
void processFields
(
    OFstream& os,
    const List<List<label>>& addr,
    const List<word>& userFieldNames,
    const IOobjectList& cloudObjs
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "steadyParticleTracksTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
