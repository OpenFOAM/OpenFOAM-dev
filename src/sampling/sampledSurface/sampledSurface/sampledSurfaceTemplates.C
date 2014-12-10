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

#include "sampledSurface.H"

template<class Type>
bool Foam::sampledSurface::checkFieldSize(const Field<Type>& field) const
{
    if (faces().empty() || field.empty())
    {
        return false;
    }

    if (field.size() != faces().size())
    {
        FatalErrorIn
        (
            "sampledSurface::checkFieldSize(const Field<Type>&) const"
        )
            << "size mismatch: "
            << "field (" << field.size()
            << ") != surface (" << faces().size() << ")"
            << exit(FatalError);
    }

    return true;
}


template<class Type>
Type Foam::sampledSurface::integrate(const Field<Type>& field) const
{
    Type value = pTraits<Type>::zero;

    if (checkFieldSize(field))
    {
        value = sum(field*magSf());
    }

    reduce(value, sumOp<Type>());
    return value;
}


template<class Type>
Type Foam::sampledSurface::integrate(const tmp<Field<Type> >& field) const
{
    Type value = integrate(field());
    field.clear();
    return value;
}


template<class Type>
Type Foam::sampledSurface::average(const Field<Type>& field) const
{
    Type value = pTraits<Type>::zero;

    if (checkFieldSize(field))
    {
        value = sum(field*magSf());
    }

    reduce(value, sumOp<Type>());

    // avoid divide-by-zero
    if (area())
    {
        return value/area();
    }
    else
    {
        return pTraits<Type>::zero;
    }
}


template<class Type>
Type Foam::sampledSurface::average(const tmp<Field<Type> >& field) const
{
    Type value = average(field());
    field.clear();
    return value;
}


template<class ReturnType, class Type>
void Foam::sampledSurface::project
(
    Field<ReturnType>& res,
    const Field<Type>& field
) const
{
    if (checkFieldSize(field))
    {
        const vectorField& norm = Sf();

        forAll(norm, faceI)
        {
            res[faceI] = field[faceI] & (norm[faceI]/mag(norm[faceI]));
        }
    }
    else
    {
        res.clear();
    }
}


template<class ReturnType, class Type>
void Foam::sampledSurface::project
(
    Field<ReturnType>& res,
    const tmp<Field<Type> >& field
) const
{
    project(res, field());
    field.clear();
}


template<class ReturnType, class Type>
Foam::tmp<Foam::Field<ReturnType> >
Foam::sampledSurface::project
(
    const tmp<Field<Type> >& field
) const
{
    tmp<Field<ReturnType> > tRes(new Field<ReturnType>(faces().size()));
    project(tRes(), field);
    return tRes;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::sampledSurface::pointAverage
(
    const GeometricField<Type, pointPatchField, pointMesh>& pfld
) const
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(pfld.mesh()());

    tmp<GeometricField<Type, fvPatchField, volMesh> > tcellAvg
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "cellAvg",
                mesh.time().timeName(),
                pfld.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& cellAvg = tcellAvg();

    labelField nPointCells(mesh.nCells(), 0);
    {
        for (label pointI = 0; pointI < mesh.nPoints(); pointI++)
        {
            const labelList& pCells = mesh.pointCells(pointI);

            forAll(pCells, i)
            {
                label cellI = pCells[i];

                cellAvg[cellI] += pfld[pointI];
                nPointCells[cellI]++;
            }
        }
    }
    forAll(cellAvg, cellI)
    {
        cellAvg[cellI] /= nPointCells[cellI];
    }
    // Give value to calculatedFvPatchFields
    cellAvg.correctBoundaryConditions();

    return tcellAvg;
}


// ************************************************************************* //
