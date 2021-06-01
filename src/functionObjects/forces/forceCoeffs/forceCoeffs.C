/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "forceCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forceCoeffs, 0);
    addToRunTimeSelectionTable(functionObject, forceCoeffs, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forceCoeffs::writeFileHeader(const label i)
{
    switch (fileID(i))
    {
        case fileID::mainFile:
        {
            // force coeff data

            writeHeader(file(i), "Force coefficients");
            writeHeaderValue(file(i), "liftDir", liftDir_);
            writeHeaderValue(file(i), "dragDir", dragDir_);
            writeHeaderValue(file(i), "pitchAxis", pitchAxis_);
            writeHeaderValue(file(i), "magUInf", magUInf_);
            writeHeaderValue(file(i), "lRef", lRef_);
            writeHeaderValue(file(i), "Aref", Aref_);
            writeHeaderValue(file(i), "CofR", coordSys_.origin());
            writeCommented(file(i), "Time");
            writeTabbed(file(i), "Cm");
            writeTabbed(file(i), "Cd");
            writeTabbed(file(i), "Cl");
            writeTabbed(file(i), "Cl(f)");
            writeTabbed(file(i), "Cl(r)");

            break;
        }
        case fileID::binsFile:
        {
            // bin coeff data

            writeHeader(file(i), "Force coefficient bins");
            writeHeaderValue(file(i), "bins", nBin_);
            writeHeaderValue(file(i), "start", binMin_);
            writeHeaderValue(file(i), "delta", binDx_);
            writeHeaderValue(file(i), "direction", binDir_);

            vectorField binPoints(nBin_);
            writeCommented(file(i), "x co-ords  :");
            forAll(binPoints, pointi)
            {
                binPoints[pointi] = (binMin_ + (pointi + 1)*binDx_)*binDir_;
                file(i) << tab << binPoints[pointi].x();
            }
            file(i) << nl;

            writeCommented(file(i), "y co-ords  :");
            forAll(binPoints, pointi)
            {
                file(i) << tab << binPoints[pointi].y();
            }
            file(i) << nl;

            writeCommented(file(i), "z co-ords  :");
            forAll(binPoints, pointi)
            {
                file(i) << tab << binPoints[pointi].z();
            }
            file(i) << nl;

            writeCommented(file(i), "Time");

            for (label j = 0; j < nBin_; j++)
            {
                const word jn('(' + Foam::name(j) + ')');
                writeTabbed(file(i), "Cm" + jn);
                writeTabbed(file(i), "Cd" + jn);
                writeTabbed(file(i), "Cl" + jn);
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled file index: " << i
                << abort(FatalError);
        }
    }

    file(i)<< endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffs::forceCoeffs
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces(name, runTime, dict),
    liftDir_(Zero),
    dragDir_(Zero),
    pitchAxis_(Zero),
    magUInf_(0.0),
    lRef_(0.0),
    Aref_(0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffs::~forceCoeffs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forceCoeffs::read(const dictionary& dict)
{
    forces::read(dict);

    // Directions for lift and drag forces, and pitch moment
    // Normalise to ensure that the directions are unit vectors

    dict.lookup("liftDir") >> liftDir_;
    liftDir_ /= mag(liftDir_);

    dict.lookup("dragDir") >> dragDir_;
    dragDir_ /= mag(dragDir_);

    dict.lookup("pitchAxis") >> pitchAxis_;
    pitchAxis_ /= mag(pitchAxis_);

    // Free stream velocity magnitude
    dict.lookup("magUInf") >> magUInf_;

    // Reference (free stream) density
    dict.lookup("rhoInf") >> rhoRef_;

    // Reference length and area scales
    dict.lookup("lRef") >> lRef_;
    dict.lookup("Aref") >> Aref_;

    return true;
}


bool Foam::functionObjects::forceCoeffs::execute()
{
    return true;
}


bool Foam::functionObjects::forceCoeffs::write()
{
    forces::calcForcesMoment();

    if (Pstream::master())
    {
        logFiles::write();

        scalar pDyn = 0.5*rhoRef_*magUInf_*magUInf_;

        Field<vector> totForce(force_[0] + force_[1] + force_[2]);
        Field<vector> totMoment(moment_[0] + moment_[1] + moment_[2]);

        List<Field<scalar>> coeffs(3);
        coeffs[0].setSize(nBin_);
        coeffs[1].setSize(nBin_);
        coeffs[2].setSize(nBin_);

        // lift, drag and moment
        coeffs[0] = (totForce & liftDir_)/(Aref_*pDyn);
        coeffs[1] = (totForce & dragDir_)/(Aref_*pDyn);
        coeffs[2] = (totMoment & pitchAxis_)/(Aref_*lRef_*pDyn);

        scalar Cl = sum(coeffs[0]);
        scalar Cd = sum(coeffs[1]);
        scalar Cm = sum(coeffs[2]);

        scalar Clf = Cl/2.0 + Cm;
        scalar Clr = Cl/2.0 - Cm;

        writeTime(file(fileID::mainFile));
        file(fileID::mainFile)
            << tab << Cm << tab  << Cd
            << tab << Cl << tab << Clf << tab << Clr << endl;

        Log << type() << " " << name() << " write:" << nl
            << "    Cm    = " << Cm << nl
            << "    Cd    = " << Cd << nl
            << "    Cl    = " << Cl << nl
            << "    Cl(f) = " << Clf << nl
            << "    Cl(r) = " << Clr << endl;

        if (nBin_ > 1)
        {
            if (binCumulative_)
            {
                for (label i = 1; i < coeffs[0].size(); i++)
                {
                    coeffs[0][i] += coeffs[0][i-1];
                    coeffs[1][i] += coeffs[1][i-1];
                    coeffs[2][i] += coeffs[2][i-1];
                }
            }

            writeTime(file(fileID::binsFile));

            forAll(coeffs[0], i)
            {
                file(fileID::binsFile)
                    << tab << coeffs[2][i]
                    << tab << coeffs[1][i]
                    << tab << coeffs[0][i];
            }

            file(fileID::binsFile) << endl;
        }

        Log << endl;
    }

    return true;
}


// ************************************************************************* //
