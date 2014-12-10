/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "OFstream.H"
#include "IOmanip.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::nastranSurfaceWriter::writeFaceValue
(
    const word& nasFieldName,
    const Type& value,
    const label EID,
    OFstream& os
) const
{
    // Fixed short/long formats:
    // 1 Nastran distributed load type, e.g. PLOAD4
    // 2 SID  : load set ID
    // 3 EID  : element ID
    // 4 onwards: load values

    label SID = 1;

    Type scaledValue = scale_*value;

    switch (writeFormat_)
    {
        case wfShort:
        {
            os.setf(ios_base::left);
            os  << setw(8) << nasFieldName;
            os.unsetf(ios_base::left);
            os.setf(ios_base::right);
            os  << setw(8) << SID
                << setw(8) << EID;

            for (direction dirI = 0; dirI < pTraits<Type>::nComponents; dirI++)
            {
                os  << setw(8) << component(scaledValue, dirI);
            }

            os.unsetf(ios_base::right);

            break;
        }
        case wfLong:
        {
            os.setf(ios_base::left);
            os  << setw(8) << word(nasFieldName + "*");
            os.unsetf(ios_base::left);
            os.setf(ios_base::right);
            os  << setw(16) << SID
                << setw(16) << EID;

            for (direction dirI = 0; dirI < pTraits<Type>::nComponents; dirI++)
            {
                os  << setw(16) << component(scaledValue, dirI);
            }

            os.unsetf(ios_base::right);

            os  << nl;

            os.setf(ios_base::left);
            os  << '*';
            os.unsetf(ios_base::left);

            break;
        }
        case wfFree:
        {
            os  << nasFieldName << ','
                << SID << ','
                << EID;

            for (direction dirI = 0; dirI < pTraits<Type>::nComponents; dirI++)
            {
                os  << ',' << component(scaledValue, dirI);
            }

            break;
        }
        default:
        {
        }
    }

    os << nl;
}


template<class Type>
void Foam::nastranSurfaceWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    if (!fieldMap_.found(fieldName))
    {
        WarningIn
        (
            "void Foam::nastranSurfaceWriter::writeTemplate"
            "("
                "const fileName&, "
                "const fileName&, "
                "const pointField&, "
                "const faceList&, "
                "const word&, "
                "const Field<Type>&, "
                "const bool, "
                "const bool"
            ") const"
        )
            << "No mapping found between field " << fieldName
            << " and corresponding Nastran field.  Available types are:"
            << fieldMap_
            << exit(FatalError);

        return;
    }

    const word& nasFieldName(fieldMap_[fieldName]);

    if (!isDir(outputDir/fieldName))
    {
        mkDir(outputDir/fieldName);
    }

    // const scalar timeValue = Foam::name(this->mesh().time().timeValue());
    const scalar timeValue = 0.0;

    OFstream os(outputDir/fieldName/surfaceName + ".dat");
    formatOS(os);

    if (verbose)
    {
        Info<< "Writing nastran file to " << os.name() << endl;
    }

    os  << "TITLE=OpenFOAM " << surfaceName.c_str() << " " << fieldName
        << " data" << nl
        << "$" << nl
        << "TIME " << timeValue << nl
        << "$" << nl
        << "BEGIN BULK" << nl;

    List<DynamicList<face> > decomposedFaces(faces.size());

    writeGeometry(points, faces, decomposedFaces, os);


    os  << "$" << nl
        << "$ Field data" << nl
        << "$" << nl;

    if (isNodeValues)
    {
        label n = 0;

        forAll(decomposedFaces, i)
        {
            const DynamicList<face>& dFaces = decomposedFaces[i];
            forAll(dFaces, faceI)
            {
                Type v = pTraits<Type>::zero;
                const face& f = dFaces[faceI];

                forAll(f, fptI)
                {
                    v += values[f[fptI]];
                }
                v /= f.size();

                writeFaceValue(nasFieldName, v, ++n, os);
            }
        }
    }
    else
    {
        label n = 0;

        forAll(decomposedFaces, i)
        {
            const DynamicList<face>& dFaces = decomposedFaces[i];

            forAll(dFaces, faceI)
            {
                writeFaceValue(nasFieldName, values[faceI], ++n, os);
            }
        }
    }

    os  << "ENDDATA" << endl;
}


// ************************************************************************* //
