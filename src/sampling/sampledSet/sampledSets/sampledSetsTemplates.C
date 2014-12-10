/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "sampledSets.H"
#include "volFields.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::sampledSets::volFieldSampler<Type>::volFieldSampler
(
    const word& interpolationScheme,
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const PtrList<sampledSet>& samplers
)
:
    List<Field<Type> >(samplers.size()),
    name_(field.name())
{
    autoPtr<interpolation<Type> > interpolator
    (
        interpolation<Type>::New(interpolationScheme, field)
    );

    forAll(samplers, setI)
    {
        Field<Type>& values = this->operator[](setI);
        const sampledSet& samples = samplers[setI];

        values.setSize(samples.size());
        forAll(samples, sampleI)
        {
            const point& samplePt = samples[sampleI];
            label cellI = samples.cells()[sampleI];
            label faceI = samples.faces()[sampleI];

            if (cellI == -1 && faceI == -1)
            {
                // Special condition for illegal sampling points
                values[sampleI] = pTraits<Type>::max;
            }
            else
            {
                values[sampleI] = interpolator().interpolate
                (
                    samplePt,
                    cellI,
                    faceI
                );
            }
        }
    }
}


template<class Type>
Foam::sampledSets::volFieldSampler<Type>::volFieldSampler
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const PtrList<sampledSet>& samplers
)
:
    List<Field<Type> >(samplers.size()),
    name_(field.name())
{
    forAll(samplers, setI)
    {
        Field<Type>& values = this->operator[](setI);
        const sampledSet& samples = samplers[setI];

        values.setSize(samples.size());
        forAll(samples, sampleI)
        {
            label cellI = samples.cells()[sampleI];

            if (cellI ==-1)
            {
                values[sampleI] = pTraits<Type>::max;
            }
            else
            {
                values[sampleI] = field[cellI];
            }
        }
    }
}


template<class Type>
Foam::sampledSets::volFieldSampler<Type>::volFieldSampler
(
    const List<Field<Type> >& values,
    const word& name
)
:
    List<Field<Type> >(values),
    name_(name)
{}


template<class Type>
void Foam::sampledSets::writeSampleFile
(
    const coordSet& masterSampleSet,
    const PtrList<volFieldSampler<Type> >& masterFields,
    const label setI,
    const fileName& timeDir,
    const writer<Type>& formatter
)
{
    wordList valueSetNames(masterFields.size());
    List<const Field<Type>*> valueSets(masterFields.size());

    forAll(masterFields, fieldi)
    {
        valueSetNames[fieldi] = masterFields[fieldi].name();
        valueSets[fieldi] = &masterFields[fieldi][setI];
    }

    fileName fName
    (
        timeDir/formatter.getFileName(masterSampleSet, valueSetNames)
    );

    OFstream ofs(fName);
    if (ofs.opened())
    {
        formatter.write
        (
            masterSampleSet,
            valueSetNames,
            valueSets,
            ofs
        );
    }
    else
    {
        WarningIn
        (
            "void Foam::sampledSets::writeSampleFile"
            "("
                "const coordSet&, "
                "const PtrList<volFieldSampler<Type> >&, "
                "const label, "
                "const fileName&, "
                "const writer<Type>&"
            ")"
        )   << "File " << ofs.name() << " could not be opened. "
            << "No data will be written" << endl;
    }
}


template<class T>
void Foam::sampledSets::combineSampledValues
(
    const PtrList<volFieldSampler<T> >& sampledFields,
    const labelListList& indexSets,
    PtrList<volFieldSampler<T> >& masterFields
)
{
    forAll(sampledFields, fieldi)
    {
        List<Field<T> > masterValues(indexSets.size());

        forAll(indexSets, setI)
        {
            // Collect data from all processors
            List<Field<T> > gatheredData(Pstream::nProcs());
            gatheredData[Pstream::myProcNo()] = sampledFields[fieldi][setI];
            Pstream::gatherList(gatheredData);

            if (Pstream::master())
            {
                Field<T> allData
                (
                    ListListOps::combine<Field<T> >
                    (
                        gatheredData,
                        Foam::accessOp<Field<T> >()
                    )
                );

                masterValues[setI] = UIndirectList<T>
                (
                    allData,
                    indexSets[setI]
                )();
            }
        }

        masterFields.set
        (
            fieldi,
            new volFieldSampler<T>
            (
                masterValues,
                sampledFields[fieldi].name()
            )
        );
    }
}


template<class Type>
void Foam::sampledSets::sampleAndWrite
(
    fieldGroup<Type>& fields
)
{
    if (fields.size())
    {
        bool interpolate = interpolationScheme_ != "cell";

        // Create or use existing writer
        if (fields.formatter.empty())
        {
            fields = writeFormat_;
        }

        // Storage for interpolated values
        PtrList<volFieldSampler<Type> > sampledFields(fields.size());

        forAll(fields, fieldi)
        {
            if (Pstream::master() && verbose_)
            {
                Pout<< "sampledSets::sampleAndWrite: "
                    << fields[fieldi] << endl;
            }

            if (loadFromFiles_)
            {
                GeometricField<Type, fvPatchField, volMesh> vf
                (
                    IOobject
                    (
                        fields[fieldi],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                );

                if (interpolate)
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<Type>
                        (
                            interpolationScheme_,
                            vf,
                            *this
                        )
                    );
                }
                else
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<Type>(vf, *this)
                    );
                }
            }
            else
            {
                if (interpolate)
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<Type>
                        (
                            interpolationScheme_,
                            mesh_.lookupObject
                            <GeometricField<Type, fvPatchField, volMesh> >
                            (fields[fieldi]),
                            *this
                        )
                    );
                }
                else
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<Type>
                        (
                            mesh_.lookupObject
                            <GeometricField<Type, fvPatchField, volMesh> >
                            (fields[fieldi]),
                            *this
                        )
                    );
                }
            }
        }

        // Combine sampled fields from processors.
        // Note: only master results are valid

        PtrList<volFieldSampler<Type> > masterFields(sampledFields.size());
        combineSampledValues(sampledFields, indexSets_, masterFields);

        if (Pstream::master())
        {
            forAll(masterSampledSets_, setI)
            {
                writeSampleFile
                (
                    masterSampledSets_[setI],
                    masterFields,
                    setI,
                    outputPath_/mesh_.time().timeName(),
                    fields.formatter()
                );
            }
        }
    }
}


// ************************************************************************* //
