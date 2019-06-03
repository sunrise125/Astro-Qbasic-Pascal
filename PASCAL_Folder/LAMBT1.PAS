  Program LAMBT1;
        { PROBLEMA di LAMBERT risolto con le VARIABILI UNIVERSALI (Metodo BATE)
          Determinazione degli elementi ORBITALI di una CONICA, NOTI:
               1   Due RAGGI Vettori r1,r2
               2 - L'angolo di trasferimento È
               3 - Il tempo di Percorrenza Á=t2-t1
  ‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹
   Centro GRAVITAZIONALE       LAMBT1.PAS -> r1,r2 in Km invece che in DU
           TERRA
  ‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹
                   }
  Uses  Crt;

    Const
         {kgauss=0.01720209895;
         uv0=29.784691694;      }

         kgauss=0.001239446678;         {Parametri TERRESTRI}
         uv0=7.905364437;

         rad=pi/180;
   Var
         tr       : String;
         segno    : Integer;

         r1,r2,tetagradi,teta,tau,tsec,
         cot,Ag, x, y, zit,z,zz, S,C,Spri,Cpri,
                          f,g,f1,g1,
         ff,gg,ff1,gg1,
         V1,V2, sigma1, fi1, cosfi, q1,q2,RT,
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
     gotoxy(25,4); write(' Raggio LOCALE Terra   [Km] = '); read(RT);
     gotoxy(25,5); write(' Quota Vettore r1      [Km] = '); read(q1);
     gotoxy(25,6); write(' Quota Vettore r2      [Km] = '); read(q2);
     gotoxy(25,7); write(' Angolo Trasfer. È [ gradi] = '); read(tetagradi);
     gotoxy(25,8); write(' Tempo di volo   Á[secondi] = '); read(tsec);
  End;


  Procedure LAMBERT;
  Label Salta;
   Begin
     r1:=1+q1/RT;  r2:=1+q2/RT; {in UnitÖ Canoniche -> RT˜DU}
       teta:=tetagradi*rad; tau:=tsec*kgauss;
          if tetagradi<=180 then segno:=1 ELSE segno:=-1;
     Cot:=COS(teta);
     amin:=(r1+r2+SQRT(r1*r1+r2*r2-2*r1*r2*cot))/4;
     Ag:=segno*SQRT(r1*r2*(1+Cot));
textcolor (10);

     z:=1;  { Valore iniziale di z=x˝/a = Parametro adimensionale }
repeat
    {Funzioni Coniche -> S,C,S',C'}
  S:=1/6*(((((z/156+1)*z/110-1)*z/72+1)*z/42-1)*z/20+1);
  C:=1/2*(((((z/132+1)*z/90-1)*z/56+1)*z/30-1)*z/12+1);
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
    F:=x*x*x*S+Ag*SQRT(y)-tau;
    F1:=x*x*x*(Spri-3*S*Cpri/(2*C))+Ag/8 *(3*S*SQRT(y)/C+Ag/x);

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

                  { Â1 = r1˘V1 = prodotto scalare
      sigma1:=(r1*r2*cot-ff*r1*r1)/gg; formula uguale a quella di sotto }
      sigma1:=(gg+r1*x*(z*S-1))/(x*x*C);

         fi1:= ArcoCOS(1/e*(p/r1-1));
         cosfi:= 1/e*(p/r1-1);
            if sigma1<0 then fi1:=2*pi-fi1;

             if e <= 0.9995 then
              begin
               E1:=ArcoCOS(1/e*(1-r1/a));
               t1menoTp:=a*SQRT(a)*(E1-sigma1/SQRT(a));
               tr:='Orbita ELLITTICA';
               Triv:=a*SQRT(a)*2*pi/kgauss;
              end
                ELSE if (e < 1.0005) and (e > 0.9995)
                 then
                  begin
                   t1menoTp:=0.5*(p*sigma1+sigma1*sigma1*sigma1/3);
                    tr:='Orbita PARABOLICA';
                  end
                 ELSE
                   begin
                   FH1:=ln(1/e*(sigma1/SQRT(-a)+(1-r1/a)));
                    t1menoTp:=-a*SQRT(-a)*(sigma1/SQRT(-a)-FH1);
                    tr:='Orbita IPERBOLICA';
                   end;
{
   if ecc>1 then begin
   yb:= sigma * h/ecc;
   bb:= (-1 / alf) * sqrt(ecc * ecc - 1);
   sx:= yb/bb;
   F1:=LN(sx + sqrt(sx * sx + 1));
   anoM:=ecc*sx - F1;
   motog:= Kgauss/rad*(-alf)*sqrt(-alf);
 writeln('          e=',ecc:10:8,'    a=',amag:10:8,' UA    n=',motog:10:8,' ¯/giorno');
 writeln('          i=',incl/rad:7:4,'¯      Í=',lnodo/rad:7:4,'¯     w=',om/rad:7:4,'¯');
 writeln('          M=',anoM/rad:7:4,'¯');
}
     gotoxy(32,12);write('RISULTATI FINALI');
     gotoxy(32,13);write('================');
     writeln;
     writeln('                     Semiasse magg.  a = ',a:11:8,' RT =',a*RT:6:0,' km');
     writeln('                     EccentricitÖ    e = ',e:11:8);
if e <= 0.9995 then writeln('                     Periodo Rivol.  P = ',Triv:11:6,' giorni');
{     writeln('                     Param.focale    p = ',p:11:8,' RT =',p*RT:6:0,' km');}
     writeln('                     Dist._periel.   q = ',q:11:8,' RT =',q*RT:6:0,' km');
     writeln('                     Quota_periel.  Qp = ',(q-1):11:8,' RT =',(q-1)*RT:6:0,' km');
     writeln('                     VelocitÖ_1     V1 = ',V1:11:8,' DU/TU = ',V1*uv0:6:3,' km/s');
     writeln('                     VelocitÖ_2     V2 = ',V2:11:8,' DU/TU = ',V2*uv0:6:3,' km/s');
     writeln('                     Anomalia_P1    Ì1 = ',fi1/rad:9:4,'¯');
     writeln('                     Pass.Perie. t1-Tp = ',t1menoTp/kgauss:6:1,' sec');
     textcolor(1);
     writeln('                                '+tr); textcolor(14);

    End;

{-----------------------------------------------------------------------------------------------
   INIZIO Programma}

   Begin
       ColoriSchermo;
       LogoDatiLAMBERT;
       LAMBERT;
       End.
   {FINE Programma}
