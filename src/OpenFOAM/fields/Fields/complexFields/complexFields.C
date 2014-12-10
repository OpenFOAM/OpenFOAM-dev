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

Description
    Specialisation of Field\<T\> for complex and complexVector.

\*---------------------------------------------------------------------------*/

#include "complexFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineCompoundTypeName(List<complex>, complexList);
addCompoundToRunTimeSelectionTable(List<complex>, complexList);

complexField ComplexField(const UList<scalar>& re, const UList<scalar>& im)
{
    complexField cf(re.size());

    forAll(cf, i)
    {
        cf[i].Re() = re[i];
        cf[i].Im() = im[i];
    }

    return cf;
}


complexField ReComplexField(const UList<scalar>& sf)
{
    complexField cf(sf.size());

    forAll(cf, i)
    {
        cf[i].Re() = sf[i];
        cf[i].Im() = 0.0;
    }

    return cf;
}


complexField ImComplexField(const UList<scalar>& sf)
{
    complexField cf(sf.size());

    forAll(cf, i)
    {
        cf[i].Re() = 0.0;
        cf[i].Im() = sf[i];
    }

    return cf;
}


scalarField ReImSum(const UList<complex>& cf)
{
    scalarField sf(cf.size());

    forAll(sf, i)
    {
        sf[i] = cf[i].Re() + cf[i].Im();
    }

    return sf;
}


scalarField Re(const UList<complex>& cf)
{
    scalarField sf(cf.size());

    forAll(sf, i)
    {
        sf[i] = cf[i].Re();
    }

    return sf;
}


scalarField Im(const UList<complex>& cf)
{
    scalarField sf(cf.size());

    forAll(sf, i)
    {
        sf[i] = cf[i].Im();
    }

    return sf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineCompoundTypeName(List<complexVector>, complexVectorList);
addCompoundToRunTimeSelectionTable(List<complexVector>, complexVectorList);

complexVectorField ComplexField
(
    const UList<vector>& re,
    const UList<vector>& im
)
{
    complexVectorField cvf(re.size());

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        forAll(cvf, i)
        {
            cvf[i].component(cmpt).Re() = re[i].component(cmpt);
            cvf[i].component(cmpt).Im() = im[i].component(cmpt);
        }
    }

    return cvf;
}


complexVectorField ReComplexField(const UList<vector>& vf)
{
    complexVectorField cvf(vf.size());

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        forAll(cvf, i)
        {
            cvf[i].component(cmpt).Re() = vf[i].component(cmpt);
            cvf[i].component(cmpt).Im() = 0.0;
        }
    }

    return cvf;
}


complexVectorField ImComplexField(const UList<vector>& vf)
{
    complexVectorField cvf(vf.size());

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        forAll(cvf, i)
        {
            cvf[i].component(cmpt).Re() = 0.0;
            cvf[i].component(cmpt).Im() = vf[i].component(cmpt);
        }
    }

    return cvf;
}


vectorField ReImSum(const UList<complexVector>& cvf)
{
    vectorField vf(cvf.size());

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        forAll(cvf, i)
        {
            vf[i].component(cmpt) =
                cvf[i].component(cmpt).Re() + cvf[i].component(cmpt).Im();
        }
    }

    return vf;
}


vectorField Re(const UList<complexVector>& cvf)
{
    vectorField vf(cvf.size());

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        forAll(cvf, i)
        {
            vf[i].component(cmpt) = cvf[i].component(cmpt).Re();
        }
    }

    return vf;
}


vectorField Im(const UList<complexVector>& cvf)
{
    vectorField vf(cvf.size());

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        forAll(cvf, i)
        {
            vf[i].component(cmpt) = cvf[i].component(cmpt).Im();
        }
    }

    return vf;
}


complexVectorField operator^
(
    const UList<vector>& vf,
    const UList<complexVector>& cvf
)
{
    return ComplexField(vf^Re(cvf), vf^Im(cvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
