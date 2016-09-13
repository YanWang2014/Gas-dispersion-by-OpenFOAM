{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->nu()/Pr + turbulence->nut()/Prt
    );

    fvScalarMatrix TSEqn
    (
        //fvm::ddt(TS)
      fvm::div(phi, TS)
      - fvm::laplacian(kappaEff, TS)
      ==
      fvOptions(TS) // 
    );

    TSEqn.relax();

    fvOptions.constrain(TSEqn);

    TSEqn.solve();
    
    fvOptions.correct(TS);
    
    //TEqn.solve().initialResidual();

    //rhok = 1.0 - beta*(T - TRef);
}
