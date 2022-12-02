/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

//extern "C"
//{
    #include "MASTER.h"
    #include "GLOBAL.h"
//}

#include "tecplotWriter.H"

#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::tecplotWriter::writeField(const Field<Type>& fld) const
{
    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        scalarField cmptFld(fld.component(cmpt));

        // Convert to float
        Field<float> floats(cmptFld.size());
        forAll(cmptFld, i)
        {
            floats[i] = float(cmptFld[i]);
        }

        INTEGER4 size = INTEGER4(floats.size());
        INTEGER4 IsDouble = 0;  // float

        // Pout<< "Writing component:" << cmpt << " of size:" << size
        //    << " floats." << endl;

        if (!TECDAT112(&size, floats.begin(), &IsDouble))
        {
//            FatalErrorInFunction
//                << "Error in TECDAT112." << exit(FatalError);
        }
    }
}


template<class Type>
Foam::tmp<Field<Type>> Foam::tecplotWriter::getPatchField
(
    const bool nearCellValue,
    const VolField<Type>& vfld,
    const label patchi
) const
{
    if (nearCellValue)
    {
        return vfld.boundaryField()[patchi].patchInternalField();
    }
    else
    {
        return vfld.boundaryField()[patchi];
    }
}


template<class Type>
Foam::tmp<Field<Type>> Foam::tecplotWriter::getFaceField
(
    const SurfaceField<Type>& sfld,
    const labelList& faceLabels
) const
{
    const polyBoundaryMesh& patches = sfld.mesh().boundaryMesh();

    tmp<Field<Type>> tfld(new Field<Type>(faceLabels.size()));
    Field<Type>& fld = tfld.ref();

    forAll(faceLabels, i)
    {
        label facei = faceLabels[i];

        label patchi = patches.whichPatch(facei);

        if (patchi == -1)
        {
            fld[i] = sfld[facei];
        }
        else
        {
            label localFacei = facei - patches[patchi].start();
            fld[i] = sfld.boundaryField()[patchi][localFacei];
        }
    }

    return tfld;
}


template<class GeoField>
Foam::wordList Foam::tecplotWriter::getNames
(
    const PtrList<GeoField>& flds
)
{
    wordList names(flds.size());
    forAll(flds, i)
    {
        names[i] = flds[i].name();
    }
    return names;
}


template<class Type>
void Foam::tecplotWriter::getTecplotNames
(
    const wordList& names,
    const INTEGER4 loc,
    string& varNames,
    DynamicList<INTEGER4>& varLocation
)
{
    forAll(names, i)
    {
        if (!varNames.empty())
        {
            varNames += " ";
        }

        label nCmpts = pTraits<Type>::nComponents;

        if (nCmpts == 1)
        {
            varNames += names[i];
            varLocation.append(loc);
        }
        else
        {
            for
            (
                direction cmpt = 0;
                cmpt < nCmpts;
                cmpt++
            )
            {
                string fldName =
                    (cmpt != 0 ? " " : string::null)
                  + names[i]
                  + "_"
                  + pTraits<Type>::componentNames[cmpt];
                varNames += fldName;
                varLocation.append(loc);
            }
        }
    }
}


template<class GeoField>
void Foam::tecplotWriter::getTecplotNames
(
    const PtrList<GeoField>& flds,
    const INTEGER4 loc,
    string& varNames,
    DynamicList<INTEGER4>& varLocation
)
{
    getTecplotNames<typename GeoField::value_type>
    (
        getNames(flds),
        loc,
        varNames,
        varLocation
    );
}


// ************************************************************************* //
