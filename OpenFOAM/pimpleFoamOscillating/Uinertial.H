  // scalar tt = runTime.value();
  dimensionedVector Vdim( "Vdim" , dimensionSet(0,1,-1,0,0,0,0), vector( 0 , 0 , Amp*Omega*Foam::sin(Omega*tt) ) );
  Info<< " just to make sure " << tt << endl;
  Uss = U + Vdim;