/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "surfaceFieldValue.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::validField
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if
    (
        regionType_ != regionTypes::sampledSurface
     && obr_.foundObject<sf>(fieldName)
    )
    {
        return true;
    }
    else if (obr_.foundObject<vf>(fieldName))
    {
        return true;
    }

    return false;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::getFieldValues
(
    const word& fieldName,
    const bool mustGet,
    const bool applyOrientation
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if
    (
        regionType_ != regionTypes::sampledSurface
     && obr_.foundObject<sf>(fieldName)
    )
    {
        return filterField(obr_.lookupObject<sf>(fieldName), applyOrientation);
    }
    else if (obr_.foundObject<vf>(fieldName))
    {
        const vf& fld = obr_.lookupObject<vf>(fieldName);

        if (surfacePtr_.valid())
        {
            if (surfacePtr_().interpolate())
            {
                const interpolationCellPoint<Type> interp(fld);
                tmp<Field<Type>> tintFld(surfacePtr_().interpolate(interp));
                const Field<Type>& intFld = tintFld();

                // Average
                const faceList& faces = surfacePtr_().faces();
                tmp<Field<Type>> tavg
                (
                    new Field<Type>(faces.size(), Zero)
                );
                Field<Type>& avg = tavg.ref();

                forAll(faces, facei)
                {
                    const face& f = faces[facei];
                    forAll(f, fp)
                    {
                        avg[facei] += intFld[f[fp]];
                    }
                    avg[facei] /= f.size();
                }

                return tavg;
            }
            else
            {
                return surfacePtr_().sample(fld);
            }
        }
        else
        {
            return filterField(fld, applyOrientation);
        }
    }

    if (mustGet)
    {
        FatalErrorInFunction
            << "Field " << fieldName << " not found in database"
            << abort(FatalError);
    }

    return tmp<Field<Type>>(new Field<Type>(0));
}


template<class Type>
Type Foam::functionObjects::fieldValues::surfaceFieldValue::
processSameTypeValues
(
    const Field<Type>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    Type result = Zero;
    switch (operation_)
    {
        case operationType::sum:
        {
            result = sum(values);
            break;
        }
        case operationType::weightedSum:
        {
            if (weightField.size())
            {
                result = sum(weightField*values);
            }
            else
            {
                result = sum(values);
            }
            break;
        }
        case operationType::sumMag:
        {
            result = sum(cmptMag(values));
            break;
        }
        case operationType::sumDirection:
        {
            FatalErrorInFunction
                << "Operation " << operationTypeNames_[operation_]
                << " not available for values of type "
                << pTraits<Type>::typeName
                << exit(FatalError);

            result = Zero;
            break;
        }
        case operationType::sumDirectionBalance:
        {
            FatalErrorInFunction
                << "Operation " << operationTypeNames_[operation_]
                << " not available for values of type "
                << pTraits<Type>::typeName
                << exit(FatalError);

            result = Zero;
            break;
        }
        case operationType::average:
        {
            result = sum(values)/values.size();
            break;
        }
        case operationType::weightedAverage:
        {
            if (weightField.size())
            {
                result =
                    sum(weightField*values)
                   /stabilise(sum(weightField), vSmall);
            }
            else
            {
                result = sum(values)/values.size();
            }
            break;
        }
        case operationType::areaAverage:
        {
            const scalarField magSf(mag(Sf));

            result = sum(magSf*values)/sum(magSf);
            break;
        }
        case operationType::weightedAreaAverage:
        {
            const scalarField magSf(mag(Sf));

            if (weightField.size())
            {
                result =
                    sum(weightField*magSf*values)
                   /stabilise(sum(magSf*weightField), vSmall);
            }
            else
            {
                result = sum(magSf*values)/sum(magSf);
            }
            break;
        }
        case operationType::areaIntegrate:
        {
            const scalarField magSf(mag(Sf));

            result = sum(magSf*values);
            break;
        }
        case operationType::weightedAreaIntegrate:
        {
            const scalarField magSf(mag(Sf));

            if (weightField.size())
            {
                result = sum(weightField*magSf*values);
            }
            else
            {
                result = sum(magSf*values);
            }
            break;
        }
        case operationType::min:
        {
            result = min(values);
            break;
        }
        case operationType::max:
        {
            result = max(values);
            break;
        }
        case operationType::CoV:
        {
            const scalarField magSf(mag(Sf));

            Type meanValue = sum(values*magSf)/sum(magSf);

            const label nComp = pTraits<Type>::nComponents;

            for (direction d=0; d<nComp; ++d)
            {
                scalarField vals(values.component(d));
                scalar mean = component(meanValue, d);
                scalar& res = setComponent(result, d);

                res = sqrt(sum(magSf*sqr(vals - mean))/sum(magSf))/mean;
            }

            break;
        }
        case operationType::areaNormalAverage:
        {}
        case operationType::areaNormalIntegrate:
        {}
        case operationType::none:
        {}
    }

    return result;
}


template<class Type>
Type Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<Type>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    return processSameTypeValues(values, Sf, weightField);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::writeValues
(
    const word& fieldName,
    const scalarField& weightField,
    const bool orient
)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> values(getFieldValues<Type>(fieldName, true, orient));

        vectorField Sf;
        if (surfacePtr_.valid())
        {
            // Get oriented Sf
            Sf = surfacePtr_().Sf();
        }
        else
        {
            // Get oriented Sf
            Sf = filterField(mesh_.Sf(), true);
        }

        // Combine onto master
        combineFields(values);
        combineFields(Sf);

        // Write raw values on surface if specified
        if (surfaceWriterPtr_.valid())
        {
            faceList faces;
            pointField points;

            if (surfacePtr_.valid())
            {
                combineSurfaceGeometry(faces, points);
            }
            else
            {
                combineMeshGeometry(faces, points);
            }

            if (Pstream::master())
            {
                surfaceWriterPtr_->write
                (
                    outputDir(),
                    regionTypeNames_[regionType_] + ("_" + regionName_),
                    points,
                    faces,
                    fieldName,
                    values,
                    false
                );
            }
        }

        if (operation_ != operationType::none)
        {
            // Apply scale factor
            values *= scaleFactor_;

            if (Pstream::master())
            {
                Type result = processValues(values, Sf, weightField);

                // Add to result dictionary, over-writing any previous entry
                resultDict_.add(fieldName, result, true);

                file() << tab << result;

                Log << "    " << operationTypeNames_[operation_]
                    << "(" << regionName_ << ") of " << fieldName
                    <<  " = " << result << endl;
            }
        }
    }

    return ok;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::filterField
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const bool applyOrientation
) const
{
    tmp<Field<Type>> tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues.ref();

    forAll(values, i)
    {
        label facei = faceId_[i];
        label patchi = facePatchId_[i];
        if (patchi >= 0)
        {
            values[i] = field.boundaryField()[patchi][facei];
        }
        else
        {
            FatalErrorInFunction
                << type() << " " << name() << ": "
                << regionTypeNames_[regionType_] << "(" << regionName_ << "):"
                << nl
                << "    Unable to process internal faces for volume field "
                << field.name() << nl << abort(FatalError);
        }
    }

    if (applyOrientation)
    {
        forAll(values, i)
        {
            values[i] *= faceSign_[i];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::filterField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
    const bool applyOrientation
) const
{
    tmp<Field<Type>> tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues.ref();

    forAll(values, i)
    {
        label facei = faceId_[i];
        label patchi = facePatchId_[i];
        if (patchi >= 0)
        {
            values[i] = field.boundaryField()[patchi][facei];
        }
        else
        {
            values[i] = field[facei];
        }
    }

    if (applyOrientation)
    {
        forAll(values, i)
        {
            values[i] *= faceSign_[i];
        }
    }

    return tvalues;
}


// ************************************************************************* //
