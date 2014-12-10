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

InNamespace
    Foam

Description
    Spatial transformation functions for primitive fields.

SourceFiles
    transformField.C

\*---------------------------------------------------------------------------*/

#ifndef transformField_H
#define transformField_H

#include "transform.H"
#include "quaternion.H"
#include "septernion.H"
#include "vectorField.H"
#include "tensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void transform(Field<Type>&, const tensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tensorField&, const tmp<Field<Type> >&);

template<class Type>
tmp<Field<Type> > transform(const tmp<tensorField>&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tmp<tensorField>&, const tmp<Field<Type> >&);


template<class Type>
void transform(Field<Type>&, const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tensor&, const tmp<Field<Type> >&);


template<class Type1, class Type2>
tmp<Field<Type1> > transformFieldMask(const Field<Type2>&);

template<class Type1, class Type2>
tmp<Field<Type1> > transformFieldMask(const tmp<Field<Type2> >&);


template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tensorField&
);

template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tmp<tensorField>&
);


template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tensorField&
);

template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tmp<tensorField>&
);


//- Rotate given vectorField with the given quaternion
void transform(vectorField&, const quaternion&, const vectorField&);

//- Rotate given vectorField with the given quaternion
tmp<vectorField> transform(const quaternion&, const vectorField&);

//- Rotate given tmp<vectorField> with the given quaternion
tmp<vectorField> transform(const quaternion&, const tmp<vectorField>&);


//- Transform given vectorField with the given septernion
void transform(vectorField&, const septernion&, const vectorField&);

//- Transform given vectorField with the given septernion
tmp<vectorField> transform(const septernion&, const vectorField&);

//- Transform given tmp<vectorField> with the given septernion
tmp<vectorField> transform(const septernion&, const tmp<vectorField>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "transformFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
