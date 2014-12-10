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

Description
    Set values on a selected set of cells/patchfaces through a dictionary.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "faceSet.H"
#include "volFields.H"

using namespace Foam;

template<class Type>
bool setCellFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedCells,
    Istream& fieldValueStream
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldTypeDesc != fieldType::typeName + "Value")
    {
        return false;
    }

    word fieldName(fieldValueStream);

    // Check the current time directory
    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check the "constant" directory
    if (!fieldHeader.headerOk())
    {
        fieldHeader = IOobject
        (
            fieldName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ
        );
    }

    // Check field exists
    if (fieldHeader.headerOk())
    {
        Info<< "    Setting internal values of "
            << fieldHeader.headerClassName()
            << " " << fieldName << endl;

        fieldType field(fieldHeader, mesh);

        const Type& value = pTraits<Type>(fieldValueStream);

        if (selectedCells.size() == field.size())
        {
            field.internalField() = value;
        }
        else
        {
            forAll(selectedCells, celli)
            {
                field[selectedCells[celli]] = value;
            }
        }

        forAll(field.boundaryField(), patchi)
        {
            field.boundaryField()[patchi] =
                field.boundaryField()[patchi].patchInternalField();
        }

        if (!field.write())
        {
            FatalErrorIn
            (
                "void setCellFieldType"
                "(const fvMesh& mesh, const labelList& selectedCells,"
                "Istream& fieldValueStream)"
            ) << "Failed writing field " << fieldName << endl;
        }
    }
    else
    {
        WarningIn
        (
            "void setCellFieldType"
            "(const fvMesh& mesh, const labelList& selectedCells,"
            "Istream& fieldValueStream)"
        ) << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    return true;
}


class setCellField
{

public:

    setCellField()
    {}

    autoPtr<setCellField> clone() const
    {
        return autoPtr<setCellField>(new setCellField());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedCells_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedCells)
        :
            mesh_(mesh),
            selectedCells_(selectedCells)
        {}

        autoPtr<setCellField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
               !(
                    setCellFieldType<scalar>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<vector>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<symmTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<tensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                )
            )
            {
                WarningIn("setCellField::iNew::operator()(Istream& is)")
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setCellField>(new setCellField());
        }
    };
};


template<class Type>
bool setFaceFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedFaces,
    Istream& fieldValueStream
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldTypeDesc != fieldType::typeName + "Value")
    {
        return false;
    }

    word fieldName(fieldValueStream);

    // Check the current time directory
    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check the "constant" directory
    if (!fieldHeader.headerOk())
    {
        fieldHeader = IOobject
        (
            fieldName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ
        );
    }

    // Check field exists
    if (fieldHeader.headerOk())
    {
        Info<< "    Setting patchField values of "
            << fieldHeader.headerClassName()
            << " " << fieldName << endl;

        fieldType field(fieldHeader, mesh);

        const Type& value = pTraits<Type>(fieldValueStream);

        // Create flat list of selected faces and their value.
        Field<Type> allBoundaryValues(mesh.nFaces()-mesh.nInternalFaces());
        forAll(field.boundaryField(), patchi)
        {
            SubField<Type>
            (
                allBoundaryValues,
                field.boundaryField()[patchi].size(),
                field.boundaryField()[patchi].patch().start()
              - mesh.nInternalFaces()
            ).assign(field.boundaryField()[patchi]);
        }

        // Override
        bool hasWarned = false;
        labelList nChanged
        (
            returnReduce(field.boundaryField().size(), maxOp<label>()),
            0
        );
        forAll(selectedFaces, i)
        {
            label facei = selectedFaces[i];
            if (mesh.isInternalFace(facei))
            {
                if (!hasWarned)
                {
                    hasWarned = true;
                    WarningIn("setFaceFieldType(..)")
                        << "Ignoring internal face " << facei
                        << ". Suppressing further warnings." << endl;
                }
            }
            else
            {
                label bFaceI = facei-mesh.nInternalFaces();
                allBoundaryValues[bFaceI] = value;
                nChanged[mesh.boundaryMesh().patchID()[bFaceI]]++;
            }
        }

        Pstream::listCombineGather(nChanged, plusEqOp<label>());
        Pstream::listCombineScatter(nChanged);

        // Reassign.
        forAll(field.boundaryField(), patchi)
        {
            if (nChanged[patchi] > 0)
            {
                Info<< "    On patch "
                    << field.boundaryField()[patchi].patch().name()
                    << " set " << nChanged[patchi] << " values" << endl;
                field.boundaryField()[patchi] == SubField<Type>
                (
                    allBoundaryValues,
                    field.boundaryField()[patchi].size(),
                    field.boundaryField()[patchi].patch().start()
                  - mesh.nInternalFaces()
                );
            }
        }

        if (!field.write())
        {
            FatalErrorIn
            (
                "void setFaceFieldType"
                "(const fvMesh& mesh, const labelList& selectedFaces,"
                "Istream& fieldValueStream)"
            )   << "Failed writing field " << field.name() << exit(FatalError);
        }
    }
    else
    {
        WarningIn
        (
            "void setFaceFieldType"
            "(const fvMesh& mesh, const labelList& selectedFaces,"
            "Istream& fieldValueStream)"
        ) << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    return true;
}


class setFaceField
{

public:

    setFaceField()
    {}

    autoPtr<setFaceField> clone() const
    {
        return autoPtr<setFaceField>(new setFaceField());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedFaces_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedFaces)
        :
            mesh_(mesh),
            selectedFaces_(selectedFaces)
        {}

        autoPtr<setFaceField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
               !(
                    setFaceFieldType<scalar>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<vector>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<symmTensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<tensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                )
            )
            {
                WarningIn("setFaceField::iNew::operator()(Istream& is)")
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setFaceField>(new setFaceField());
        }
    };
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    Info<< "Reading setFieldsDict\n" << endl;

    IOdictionary setFieldsDict
    (
        IOobject
        (
            "setFieldsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    if (setFieldsDict.found("defaultFieldValues"))
    {
        Info<< "Setting field default values" << endl;
        PtrList<setCellField> defaultFieldValues
        (
            setFieldsDict.lookup("defaultFieldValues"),
            setCellField::iNew(mesh, labelList(mesh.nCells()))
        );
        Info<< endl;
    }


    Info<< "Setting field region values" << endl;

    PtrList<entry> regions(setFieldsDict.lookup("regions"));

    forAll(regions, regionI)
    {
        const entry& region = regions[regionI];

        autoPtr<topoSetSource> source =
            topoSetSource::New(region.keyword(), mesh, region.dict());

        if (source().setType() == topoSetSource::CELLSETSOURCE)
        {
            cellSet selectedCellSet
            (
                mesh,
                "cellSet",
                mesh.nCells()/10+1  // Reasonable size estimate.
            );

            source->applyToSet
            (
                topoSetSource::NEW,
                selectedCellSet
            );

            PtrList<setCellField> fieldValues
            (
                region.dict().lookup("fieldValues"),
                setCellField::iNew(mesh, selectedCellSet.toc())
            );
        }
        else if (source().setType() == topoSetSource::FACESETSOURCE)
        {
            faceSet selectedFaceSet
            (
                mesh,
                "faceSet",
                (mesh.nFaces()-mesh.nInternalFaces())/10+1
            );

            source->applyToSet
            (
                topoSetSource::NEW,
                selectedFaceSet
            );

            PtrList<setFaceField> fieldValues
            (
                region.dict().lookup("fieldValues"),
                setFaceField::iNew(mesh, selectedFaceSet.toc())
            );
        }
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
