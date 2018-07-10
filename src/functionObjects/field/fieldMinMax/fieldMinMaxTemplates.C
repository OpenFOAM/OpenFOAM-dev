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

#include "fieldMinMax.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldMinMax::output
(
    const word& fieldName,
    const word& outputName,
    const label minCell,
    const label maxCell,
    const vector& minC,
    const vector& maxC,
    const label minProci,
    const label maxProci,
    const Type& minValue,
    const Type& maxValue
)
{
    OFstream& file = this->file();

    if (location_)
    {
        writeTime(file());

        writeTabbed(file, fieldName);

        file<< tab << minValue
            << tab << minC;

        if (Pstream::parRun())
        {
            file<< tab << minProci;
        }

        file<< tab << maxValue
            << tab << maxC;

        if (Pstream::parRun())
        {
            file<< tab << maxProci;
        }

        file<< endl;

        Log << "    min(" << outputName << ") = " << minValue
            << " in cell " << minCell
            << " at location " << minC;

        if (Pstream::parRun())
        {
            Log << " on processor " << minProci;
        }

        Log << nl << "    max(" << outputName << ") = " << maxValue
            << " in cell " << maxCell
            << " at location " << maxC;

        if (Pstream::parRun())
        {
            Log << " on processor " << maxProci;
        }
    }
    else
    {
        file<< tab << minValue << tab << maxValue;

        Log << "    min/max(" << outputName << ") = "
            << minValue << ' ' << maxValue;
    }

    Log << endl;
}


template<class Type>
void Foam::functionObjects::fieldMinMax::calcMinMaxFields
(
    const word& fieldName,
    const modeType& mode
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const label proci = Pstream::myProcNo();

        const fieldType& field = obr_.lookupObject<fieldType>(fieldName);

        const volVectorField::Boundary& CfBoundary =
            mesh_.C().boundaryField();

        switch (mode)
        {
            case mdMag:
            {
                const volScalarField magField(mag(field));
                const volScalarField::Boundary& magFieldBoundary =
                    magField.boundaryField();

                scalarList minVs(Pstream::nProcs());
                labelList minCells(Pstream::nProcs());
                List<vector> minCs(Pstream::nProcs());
                label minProci = findMin(magField);
                minVs[proci] = magField[minProci];
                minCells[proci] = minProci;
                minCs[proci] = mesh_.C()[minProci];

                scalarList maxVs(Pstream::nProcs());
                labelList maxCells(Pstream::nProcs());
                List<vector> maxCs(Pstream::nProcs());
                label maxProci = findMax(magField);
                maxVs[proci] = magField[maxProci];
                maxCells[proci] = maxProci;
                maxCs[proci] = mesh_.C()[maxProci];

                forAll(magFieldBoundary, patchi)
                {
                    const scalarField& mfp = magFieldBoundary[patchi];
                    if (mfp.size())
                    {
                        const vectorField& Cfp = CfBoundary[patchi];

                        const labelList& faceCells =
                            magFieldBoundary[patchi].patch().faceCells();

                        label minPI = findMin(mfp);
                        if (mfp[minPI] < minVs[proci])
                        {
                            minVs[proci] = mfp[minPI];
                            minCells[proci] = faceCells[minPI];
                            minCs[proci] = Cfp[minPI];
                        }

                        label maxPI = findMax(mfp);
                        if (mfp[maxPI] > maxVs[proci])
                        {
                            maxVs[proci] = mfp[maxPI];
                            maxCells[proci] = faceCells[maxPI];
                            maxCs[proci] = Cfp[maxPI];
                        }
                    }
                }

                Pstream::gatherList(minVs);
                Pstream::gatherList(minCells);
                Pstream::gatherList(minCs);

                Pstream::gatherList(maxVs);
                Pstream::gatherList(maxCells);
                Pstream::gatherList(maxCs);

                if (Pstream::master())
                {
                    const label minI = findMin(minVs);
                    const scalar minValue = minVs[minI];
                    const label minCell = minCells[minI];
                    const vector& minC = minCs[minI];

                    const label maxI = findMax(maxVs);
                    const scalar maxValue = maxVs[maxI];
                    const label maxCell = maxCells[maxI];
                    const vector& maxC = maxCs[maxI];

                    output
                    (
                        fieldName,
                        word("mag(" + fieldName + ")"),
                        minCell,
                        maxCell,
                        minC,
                        maxC,
                        minI,
                        maxI,
                        minValue,
                        maxValue
                    );
                }
                break;
            }
            case mdCmpt:
            {
                const typename fieldType::Boundary&
                    fieldBoundary = field.boundaryField();

                List<Type> minVs(Pstream::nProcs());
                labelList minCells(Pstream::nProcs());
                List<vector> minCs(Pstream::nProcs());
                label minProci = findMin(field);
                minVs[proci] = field[minProci];
                minCells[proci] = minProci;
                minCs[proci] = mesh_.C()[minProci];

                List<Type> maxVs(Pstream::nProcs());
                labelList maxCells(Pstream::nProcs());
                List<vector> maxCs(Pstream::nProcs());
                label maxProci = findMax(field);
                maxVs[proci] = field[maxProci];
                maxCells[proci] = maxProci;
                maxCs[proci] = mesh_.C()[maxProci];

                forAll(fieldBoundary, patchi)
                {
                    const Field<Type>& fp = fieldBoundary[patchi];

                    if (fp.size())
                    {
                        const vectorField& Cfp = CfBoundary[patchi];

                        const labelList& faceCells =
                            fieldBoundary[patchi].patch().faceCells();

                        label minPI = findMin(fp);
                        if (fp[minPI] < minVs[proci])
                        {
                            minVs[proci] = fp[minPI];
                            minCells[proci] = faceCells[minPI];
                            minCs[proci] = Cfp[minPI];
                        }

                        label maxPI = findMax(fp);
                        if (fp[maxPI] > maxVs[proci])
                        {
                            maxVs[proci] = fp[maxPI];
                            maxCells[proci] = faceCells[maxPI];
                            maxCs[proci] = Cfp[maxPI];
                        }
                    }
                }

                Pstream::gatherList(minVs);
                Pstream::gatherList(minCells);
                Pstream::gatherList(minCs);

                Pstream::gatherList(maxVs);
                Pstream::gatherList(maxCells);
                Pstream::gatherList(maxCs);

                if (Pstream::master())
                {
                    label minI = findMin(minVs);
                    Type minValue = minVs[minI];
                    const label minCell = minCells[minI];
                    const vector& minC = minCs[minI];

                    label maxI = findMax(maxVs);
                    Type maxValue = maxVs[maxI];
                    const label maxCell = maxCells[maxI];
                    const vector& maxC = maxCs[maxI];

                    output
                    (
                        fieldName,
                        fieldName,
                        minCell,
                        maxCell,
                        minC,
                        maxC,
                        minI,
                        maxI,
                        minValue,
                        maxValue
                    );
                }
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown min/max mode: " << modeTypeNames_[mode_]
                    << exit(FatalError);
            }
        }
    }
}


// ************************************************************************* //
