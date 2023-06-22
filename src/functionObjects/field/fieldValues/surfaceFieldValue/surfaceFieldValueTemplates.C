/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    if
    (
        selectionType_ != selectionTypes::sampledSurface
     && obr_.foundObject<SurfaceField<Type>>(fieldName)
    )
    {
        return true;
    }
    else if (obr_.foundObject<VolField<Type>>(fieldName))
    {
        return true;
    }

    return false;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::getFieldValues
(
    const word& fieldName
) const
{
    if (selectionType_ == selectionTypes::sampledSurface)
    {
        if (obr_.foundObject<VolField<Type>>(fieldName))
        {
            const VolField<Type>& fld =
                obr_.lookupObject<VolField<Type>>(fieldName);

            if (surfacePtr_().interpolate())
            {
                // Interpolate the field to the surface points
                const interpolationCellPoint<Type> interp(fld);
                tmp<Field<Type>> tintFld(surfacePtr_().interpolate(interp));
                const Field<Type>& intFld = tintFld();

                // Average the interpolated field onto the surface faces
                const faceList& faces = surfacePtr_().faces();
                tmp<Field<Type>> tavg(new Field<Type>(faces.size(), Zero));
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
        else if (obr_.foundObject<SurfaceField<Type>>(fieldName))
        {
            FatalErrorInFunction
                << "Surface field " << fieldName
                << " cannot be sampled onto surface " << surfacePtr_().name()
                << ". Only vol fields can be sampled onto surfaces."
                << abort(FatalError);
        }
    }
    else
    {
        if (obr_.foundObject<VolField<Type>>(fieldName))
        {
            const VolField<Type>& fld =
                obr_.lookupObject<VolField<Type>>(fieldName);
            return filterField(fld);
        }
        else if (obr_.foundObject<SurfaceField<Type>>(fieldName))
        {
            const SurfaceField<Type>& fld =
                obr_.lookupObject<SurfaceField<Type>>(fieldName);
            return filterField(fld);
        }
    }

    FatalErrorInFunction
        << "Field " << fieldName << " not found in database"
        << abort(FatalError);

    return tmp<Field<Type>>(nullptr);
}


template<class Type, class ResultType>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<Type>& values,
    const scalarField& signs,
    const scalarField& weights,
    const vectorField& Sf,
    ResultType& result
) const
{
    return false;
}


template<class Type>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<Type>& values,
    const scalarField& signs,
    const scalarField& weights,
    const vectorField& Sf,
    Type& result
) const
{
    return processValuesTypeType(values, signs, weights, Sf, result);
}


template<class Type>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<Type>& values,
    const scalarField& signs,
    const scalarField& weights,
    const vectorField& Sf,
    scalar& result
) const
{
    switch (operation_)
    {
        case operationType::minMag:
        {
            result = gMin(mag(values));
            return true;
        }
        case operationType::maxMag:
        {
            result = gMax(mag(values));
            return true;
        }
        default:
        {
            // No fall through
            return false;
        }
    }
}


template<class Type>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::
processValuesTypeType
(
    const Field<Type>& values,
    const scalarField& signs,
    const scalarField& weights,
    const vectorField& Sf,
    Type& result
) const
{
    switch (operation_)
    {
        case operationType::sum:
        {
            result = gSum(weights*values);
            return true;
        }
        case operationType::sumMag:
        {
            result = gSum(weights*cmptMag(values));
            return true;
        }
        case operationType::orientedSum:
        {
            result = gSum(signs*weights*values);
            return true;
        }
        case operationType::average:
        {
            result =
                gSum(weights*values)
               /stabilise(gSum(weights), vSmall);
            return true;
        }
        case operationType::areaAverage:
        {
            const scalarField magSf(mag(Sf));
            result =
                gSum(weights*magSf*values)
               /stabilise(gSum(weights*magSf), vSmall);
            return true;
        }
        case operationType::areaIntegrate:
        {
            const scalarField magSf(mag(Sf));
            result = gSum(weights*magSf*values);
            return true;
        }
        case operationType::min:
        {
            result = gMin(values);
            return true;
        }
        case operationType::max:
        {
            result = gMax(values);
            return true;
        }
        case operationType::CoV:
        {
            const scalarField magSf(mag(Sf));

            Type meanValue = gSum(values*magSf)/gSum(magSf);

            const label nComp = pTraits<Type>::nComponents;

            for (direction d=0; d<nComp; ++d)
            {
                scalarField vals(values.component(d));
                scalar mean = component(meanValue, d);
                scalar& res = setComponent(result, d);

                res = sqrt(gSum(magSf*sqr(vals - mean))/gSum(magSf))/mean;
            }

            return true;
        }
        case operationType::none:
        {
            return true;
        }
        default:
        {
            return false;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::writeValues
(
    const word& fieldName,
    const scalarField& signs,
    const scalarField& weights,
    const vectorField& Sf
)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        // Get the values
        Field<Type> values(getFieldValues<Type>(fieldName));

        // Write raw values on surface if specified
        if (writeFields_)
        {
            faceList faces;
            pointField points;

            if (selectionType_ == selectionTypes::sampledSurface)
            {
                combineSurfaceGeometry(faces, points);
            }
            else
            {
                combineMeshGeometry(faces, points);
            }

            Field<Type> writeValues(weights*values);
            combineFields(writeValues);

            if (Pstream::master())
            {
                surfaceWriterPtr_->write
                (
                    outputDir(),
                    fieldName
                  + '_' + selectionTypeNames[selectionType_]
                  + '_' + selectionName_,
                    points,
                    faces,
                    false,
                    fieldName,
                    writeValues
                );
            }
        }

        // Do the operation
        if (operation_ != operationType::none)
        {
            // Apply scale factor
            values *= scaleFactor_;

            bool ok = false;

            #define writeValuesFieldType(fieldType, none)                      \
                ok =                                                           \
                    ok                                                         \
                 || writeValues<Type, fieldType>                               \
                    (                                                          \
                        fieldName,                                             \
                        values,                                                \
                        signs,                                                 \
                        weights,                                               \
                        Sf                                                     \
                    );
            FOR_ALL_FIELD_TYPES(writeValuesFieldType);
            #undef writeValuesFieldType

            if (!ok)
            {
                FatalErrorInFunction
                    << "Operation " << operationTypeNames_[operation_]
                    << " not available for values of type "
                    << pTraits<Type>::typeName
                    << exit(FatalError);
            }
        }
    }

    return ok;
}


template<class Type, class ResultType>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::writeValues
(
    const word& fieldName,
    const Field<Type>& values,
    const scalarField& signs,
    const scalarField& weights,
    const vectorField& Sf
)
{
    ResultType result;

    if (processValues(values, signs, weights, Sf, result))
    {
        // Add to result dictionary, over-writing any previous entry
        resultDict_.add(fieldName, result, true);

        if (Pstream::master())
        {
            file() << tab << result;

            Log << "    " << operationTypeNames_[operation_]
                << "(" << selectionName_ << ") of " << fieldName
                <<  " = " << result << endl;
        }

        return true;
    }

    return false;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::filterField
(
    const VolField<Type>& field
) const
{
    tmp<Field<Type>> tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues.ref();

    forAll(values, i)
    {
        const label facei = faceId_[i];
        const label patchi = facePatchId_[i];

        if (patchi >= 0)
        {
            values[i] = field.boundaryField()[patchi][facei];
        }
        else
        {
            FatalErrorInFunction
                << type() << " " << name() << ": "
                << selectionTypeNames[selectionType_]
                << "(" << selectionName_ << "):"
                << nl
                << "    Unable to process internal faces for volume field "
                << field.name() << nl << abort(FatalError);
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::filterField
(
    const SurfaceField<Type>& field
) const
{
    tmp<Field<Type>> tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues.ref();

    forAll(values, i)
    {
        const label facei = faceId_[i];
        const label patchi = facePatchId_[i];

        if (patchi >= 0)
        {
            values[i] = field.boundaryField()[patchi][facei];
        }
        else
        {
            values[i] = field[facei];
        }
    }

    return tvalues;
}


// ************************************************************************* //
