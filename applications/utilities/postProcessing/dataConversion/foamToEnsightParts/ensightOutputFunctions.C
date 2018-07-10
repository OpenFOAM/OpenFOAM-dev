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

#include "ensightOutputFunctions.H"

#include "passiveParticle.H"
#include "IOField.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "OFstream.H"
#include "IOmanip.H"


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::ensightCaseEntry
(
    OFstream& caseFile,
    const string& ensightType,
    const word& fieldName,
    const fileName& dataMask,
    const fileName& local,
    const label cloudNo,
    const label timeSet
)
{
    caseFile.setf(ios_base::left);

    fileName dirName(dataMask);
    if (local.size())
    {
        dirName = dirName/local;
    }

    if (cloudNo >= 0)
    {
        label ts = 1;
        if (timeSet > ts)
        {
            ts = timeSet;
        }

        // prefix variables with 'c' (cloud)
        caseFile
            << ensightType.c_str()
            << " per measured node: " << ts << " "
            << setw(15)
            << ("c" + Foam::name(cloudNo) + fieldName).c_str()
            << " "
            << (dirName/fieldName).c_str()
            << nl;
    }
    else
    {
        caseFile
            << ensightType.c_str()
            << " per element: "
            << setw(15) << fieldName
            << " "
            << (dirName/fieldName).c_str()
            << nl;
    }
}


void Foam::ensightParticlePositions
(
    const polyMesh& mesh,
    const fileName& dataDir,
    const fileName& subDir,
    const word& cloudName,
    IOstream::streamFormat format
)
{
    Cloud<passiveParticle> parcels(mesh, cloudName, false);

    fileName cloudDir = subDir/cloud::prefix/cloudName;
    fileName postFileName = cloudDir/"positions";

    // the ITER/lagrangian subdirectory must exist
    mkDir(dataDir/cloudDir);
    ensightFile os(dataDir/postFileName, format);

    // tag binary format (just like geometry files)
    os.writeBinaryHeader();
    os.write(postFileName);
    os.newline();
    os.write("particle coordinates");
    os.newline();
    os.write(parcels.size(), 8);   // unusual width
    os.newline();

    // binary write is Ensight6 - first ids, then positions
    if (format == IOstream::BINARY)
    {
        forAll(parcels, i)
        {
            os.write(i+1);
        }

        forAllConstIter(Cloud<passiveParticle>, parcels, elmnt)
        {
            const vector& p = elmnt().position();

            os.write(p.x());
            os.write(p.y());
            os.write(p.z());
        }
    }
    else
    {
        label nParcels = 0;

        forAllConstIter(Cloud<passiveParticle>, parcels, elmnt)
        {
            const vector& p = elmnt().position();

            os.write(++nParcels, 8);    // unusual width
            os.write(p.x());
            os.write(p.y());
            os.write(p.z());
            os.newline();
        }
    }
}



template<class Type>
void Foam::ensightLagrangianField
(
    const IOobject& fieldObject,
    const fileName& dataDir,
    const fileName& subDir,
    const word& cloudName,
    IOstream::streamFormat format
)
{
    Info<< " " << fieldObject.name() << flush;

    fileName cloudDir = subDir/cloud::prefix/cloudName;
    fileName postFileName = cloudDir/fieldObject.name();

    string title =
        postFileName + " with " + pTraits<Type>::typeName + " values";

    ensightFile os(dataDir/postFileName, format);
    os.write(title);
    os.newline();

    IOField<Type> field(fieldObject);

    // 6 values per line
    label count = 0;

    forAll(field, i)
    {
        Type val = field[i];

        if (mag(val) < 1.0e-90)
        {
            val = Zero;
        }

        for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            os.write( component(val, cmpt) );
        }

        count += pTraits<Type>::nComponents;

        if (count % 6 == 0)
        {
            os.newline();
        }
    }

    // add final newline if required
    if (count % 6)
    {
        os.newline();
    }
}


template<class Type>
void Foam::ensightVolField
(
    const ensightParts& partsList,
    const IOobject& fieldObject,
    const fvMesh& mesh,
    const fileName& dataDir,
    const fileName& subDir,
    IOstream::streamFormat format
)
{
    Info<< " " << fieldObject.name() << flush;

    fileName postFileName = subDir/fieldObject.name();

    ensightFile os(dataDir/postFileName, format);
    os.write(postFileName);
    os.newline();

    // ie, volField<Type>
    partsList.writeField
    (
        os,
        GeometricField<Type, fvPatchField, volMesh>
        (
            fieldObject,
            mesh
        )
    );
}


// ************************************************************************* //
