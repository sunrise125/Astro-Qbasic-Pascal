Program SWIFT2;     { Programma PRINCIPALE di ®Code SWIFT¯ di Duncan e Levison
                      mutuato dal FORTRAN degli stessi autori - Mag.97
    Routine per passaggio
    Vettori Posiz./Veloc. =>  Elem.Orbitali
}

  Uses  Crt,UNIT_2;

   Const
         rad=pi/180; kgauss=0.01720209895;
   Var
      x,y,z,vx,vy,vz,gmsum : Real;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
  End;

  Procedure Titolo;
   Begin;
     textcolor(10);
     gotoxy(8,2); writeln('CALCOLO del CODICE ®SWIFT¯ di Duncan/Levison (2^parte: ORBEL_XV2EL)');
     gotoxy(8,3); writeln('-------------------------------------------------------------------');
     textcolor(14);
  End;
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
       ColoriSchermo;
       Titolo;
{*****************************************************************************}
    gmsum:=SQR(kgauss);
    x:=-0.770320619;  y:=0.6167274476;  z:=-0.0000020354;
    vx:=-0.0110308591;  vy:=-0.0134877323;  vz:=-0.0000003624;

  gotoxy(3,4); writeln('1ø Esempio: Terra al 10/02/2001 Equin. DATA');
     textcolor(10);
  gotoxy(8,5); writeln('          x=',x:11:8,'   y=',y:11:8,'   z=',z:11:8);
  gotoxy(8,6); writeln('         vx=',vx:11:8,'  vy=',vy:11:8,'  vz=',vz:11:8);
     textcolor(14);

       ORBEL_XV2EL(x,y,z,vx,vy,vz,gmsum);
  writeln('     a= ',a:9:7,'  e= ',e:9:7,'   i=',inc/rad:9:5,'ø  ê=',capom/rad:9:5,'ø  w=',omega/rad:9:5,'ø');
  writeln('     M= ',capm/rad:9:5,'ø');
{*****************************************************************************}
    gmsum:=SQR(kgauss);
    x:=-0.7701813701; y:=0.6169433349;  z:=-0.0000011083;
    vx:=-0.0110363093;  vy:=-0.0134915447;  vz:=-0.0000000362;

  gotoxy(3,10); writeln('2ø Esempio: BARICENTRO (Terra+Luna) al 10/02/2001 Eq.J2000');
     textcolor(10);
  gotoxy(8,11); writeln('          x=',x:11:8,'   y=',y:11:8,'   z=',z:11:8);
  gotoxy(8,12); writeln('         vx=',vx:11:8,'  vy=',vy:11:8,'  vz=',vz:11:8);
     textcolor(14);

       ORBEL_XV2EL(x,y,z,vx,vy,vz,gmsum);
  writeln('     a= ',a:9:7,'  e= ',e:9:7,'   i=',inc/rad:9:5,'ø  ê=',capom/rad:9:5,'ø  w=',omega/rad:9:5,'ø');
  writeln('     M= ',capm/rad:9:5,'ø');
{*****************************************************************************}
    gmsum:=1;
  x:=0.75716461; y:=-0.67363699; z:=0;
  vx:=0.515622285; vy:=0.626855851; vz:=-0.284937043;
  gotoxy(3,16); writeln('3ø Esempio: Orbita ELIOcentrica METEORA)');
     textcolor(10);
  gotoxy(8,17); writeln('          x=',x:11:8,'   y=',y:11:8,'   z=',z:11:8);
  gotoxy(8,18); writeln('         vx=',vx:11:8,'  vy=',vy:11:8,'  vz=',vz:11:8);
     textcolor(14);

       ORBEL_XV2EL(x,y,z,vx,vy,vz,gmsum);
  writeln('     a= ',a:9:7,'  e= ',e:9:7,'   i=',inc/rad:9:5,'ø  ê=',capom/rad:9:5,'ø  w=',omega/rad:9:5,'ø');
  writeln('     M= ',capm/rad:9:5,'ø');
{*****************************************************************************}   
    repeat
     ch:=readkey;
       if ch='q' then
         begin
          exit;
       end;
  until keypressed;   

     
     End.
   {FINE Programma}
