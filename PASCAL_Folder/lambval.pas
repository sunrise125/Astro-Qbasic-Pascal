  Program LAMBVAL;
        { PROBLEMA di LAMBERT risolto con le VARIABILI UNIVERSALI (Metodo BATE)
          Determinazione degli elementi ORBITALI di una CONICA, NOTI:
               1   Due RAGGI Vettori r1,r2 in KM
               2 - L'angolo di trasferimento Theta in gradi
               3 - Il tempo di Percorrenza tau=t2-t1 in secondi
               (Esempio di Vallado: MU=398600.44018 km^3/s^2)
  ‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹
   Centro GRAVITAZIONALE
           TERRA
  ‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹
                   }
  Uses  Crt;

    Const
        
         MU=398600.44018;
         {kgauss=0.001239446678;         Parametri TERRESTRI}
         

         rad=pi/180;
   Var
         tr       : String;
         segno    : Integer;

         r1,r2,tetagradi,teta,tau,tsec,
         cot,Ag, x, y, zit,z,zz, S,C,Spri,Cpri,
                          f,g,f1,g1,
         ff,gg,ff1,gg1,
         V1,V2, sigma1, fi1, cosfi,
         E1,t1menoTp,FH1,
         a,e,p,q,Triv,amin       : Real;


   Function ArcoCOS(valore:Real): Real;
     Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(abs(1-SQR(valore))));
     End;

      Function ArcoTANG(num,den:Real) : Real;
   Var atan : Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
  End;

  Procedure LogoDatiLAMBERT;
   Begin
     textcolor(10);
     gotoxy(7,2); writeln('PROBLEMA di LAMBERT risolto con le VARIABILI UNIVERSALI (Metodo BATE)');
     gotoxy(7,3); writeln('---------------------------------------------------------------------');
     textcolor(14);
     gotoxy(25,5); write(' Raggio Vettore r1     [Km] = '); read(r1);
     gotoxy(25,6); write(' Raggio Vettore r2     [Km] = '); read(r2);
     gotoxy(25,7); write(' Angolo Trasf.teta [ gradi] = '); read(tetagradi);
     gotoxy(25,8); write(' Tempo di volo tau[secondi] = '); read(tsec);
  End;


  Procedure LAMBERT;
  Label Salta;
   Begin
       teta:=tetagradi*rad; tau:=tsec;
          if tetagradi<=180 then segno:=1 ELSE segno:=-1;
     Cot:=COS(teta);
     amin:=(r1+r2+SQRT(r1*r1+r2*r2-2*r1*r2*cot))/4;
     Ag:=segno*SQRT(r1*r2*(1+Cot));
textcolor (10);

     z:=15000;  { Valore iniziale di z=x^2/a = Parametro adimensionale }
repeat
    {Funzioni Coniche -> S,C,S',C'}
  S:=1/6*((((((((z/342-1)*z/272+1)*z/210-1)*z/156+1)*z/110-1)*z/72+1)*z/42-1)*z/20+1);
  C:=1/2*((((((((z/306-1)*z/240+1)*z/182-1)*z/132+1)*z/90-1)*z/56+1)*z/30-1)*z/12+1);
  Spri:=(C-3*S)/(2*z);
  Cpri:=(1-z*S-2*C)/(2*z);
     { Anomalia eccentrica x e parametro y }
  y:=r1+r2-Ag*(1-z*S)/SQRT(C);
     if y<0 then
        begin
          z:=z/10;
        goto Salta
     end;

  x:=SQRT(y/C);
    F:=x*x*x*S+Ag*SQRT(y)-tau*SQRT(MU);
    F1:=1/SQRT(MU)*(x*x*x*(Spri-3*S*Cpri/(2*C))+Ag/8 *(3*S*SQRT(y)/C + Ag/x));

     zit:= z - F / F1;   { Iterazione di NEWTON }

   zz:=ABS(zit-z);
   z:=zit;
     gotoxy(25,10); write('    z = ',zit:10:8,'   x =',x:10:8);
salta:
until zz < 0.00000001;
textcolor (14);
       { Semiasse maggiore a, parametro focale p, eccentricitÖ e }
       a:=x*x/z;  p:=r1*r2/y*(1-cot); e:=SQRT(1-p/a); q:=a*(1-e);

         { VelocitÖ agli estremi dell'arco }
       V1:=SQRT(2/r1-1/a);  V2:=SQRT(2/r2-1/a);

         { Parametri di Lagrange: f,g,f',g'}
      ff:=1-y/r1;  gg:=Ag*SQRT(y);
      gg1:=1-y/r2; ff1:=(ff*gg1-1)/gg;

 
    End;

{-----------------------------------------------------------------------------------------------
   INIZIO Programma}

   Begin
       ColoriSchermo;
       LogoDatiLAMBERT;
       LAMBERT;
       End.
   {FINE Programma}
