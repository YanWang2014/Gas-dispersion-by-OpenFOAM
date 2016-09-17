    Info<< "Reading field k\n" << endl;

    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p\n" << endl;

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
		dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0), vector::zero) //additional argument
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


//    dimensionedScalar k
//    (
//        transportProperties.lookup("k")   //permeability
//    );

    dimensionedScalar mu
    (
        transportProperties.lookup("mu")   //fluid dynamic vicosity
    );
    
        dimensionedScalar S
    (
        transportProperties.lookup("S")  
    );

        dimensionedScalar RHO
    (
        transportProperties.lookup("RHO")
    );
    
        dimensionedScalar nS
    (
        transportProperties.lookup("nS")
    );

//        volScalarField n
//        (
//            IOobject
//            (
//                "n",
//                runTime.timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mesh,
//            dimensionedScalar("n", dimensionSet(0,0,0,0,0,0,0), 1.0) //additional argument
//        );

        label pRefCell = 0;
        scalar pRefValue = 0;
        setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);


