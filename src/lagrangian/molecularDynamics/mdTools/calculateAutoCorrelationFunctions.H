/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

if (mesh.time().timeIndex() % vacf.sampleSteps() == 0)
{
    Field<vector> uVals(molecules.size());

    label uV = 0;

    forAllConstIter(IDLList<molecule>, molecules, mol)
    {
        uVals[uV++] = mol().U();
    }

    vacf.calculateCorrelationFunction(uVals);
}

if (mesh.time().timeIndex() % pacf.sampleSteps() == 0)
{
    vector p = Zero;

    forAllConstIter(IDLList<molecule>, molecules, mol)
    {
        p.x() +=
            mol().mass() * mol().U().y() * mol().U().z()
          + 0.5*mol().rf().yz();

        p.y() +=
            mol().mass() * mol().U().z() * mol().U().x()
          + 0.5*mol().rf().zx();

        p.z() +=
            mol().mass() * mol().U().x() * mol().U().y()
          + 0.5*mol().rf().xy();
    }

    pacf.calculateCorrelationFunction(p);
}

if (mesh.time().timeIndex() % hfacf.sampleSteps() == 0)
{
    vector s = Zero;

    forAllConstIter(IDLList<molecule>, molecules, mol)
    {
        s +=
        (
            0.5*mol().mass()*magSqr(mol().U())
          + mol().potentialEnergy()
        )*mol().U()
      + 0.5*(mol().rf() & mol().U());
    }

    hfacf.calculateCorrelationFunction(s);
}
