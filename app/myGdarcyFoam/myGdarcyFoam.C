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

Application
    mydarcyFoam  (laplacianFoam)

Description
    Solves a simple Darcy equation, e.g. for flow in porous media.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

#include "fvIOoptionList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    #include "readGravitationalAcceleration.H"
    #include "createFvOptions.H" //initialize  a global list of fvOptions to be used with this solver.


    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating pressure distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {

            fvScalarMatrix pEqn
            (
                        fvm::laplacian(k/mu, p)
                                - fvc::div(RHO*k/mu*g)  //DENSITY = 1
                                ==
                                -nS*fvOptions(p) //S
            );

            pEqn.setReference(pRefCell, pRefValue);

            fvOptions.constrain(pEqn);

            pEqn.solve();

            fvOptions.correct(p);

            U = fvc::reconstruct(-fvc::interpolate(k)/mu*
                                  (fvc::snGrad(p)*mesh.magSf()-
                                   ((RHO*g)
                                               & mesh.Sf()))
                                 );

        }


        //U = fvc::reconstruct(-fvc::interpolate(k)/mu*fvc::snGrad(p)*mesh.magSf());



        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
