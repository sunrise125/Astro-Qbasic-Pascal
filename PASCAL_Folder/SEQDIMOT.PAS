Program SEQDIMOT; {Calcolo della MINIMA Distanza tra 2 Curve nello SPAZIO
                   Serve per il calcolo dell'avvicinamento minimo POSSIBILE
                   alla Terra degli Asteroidi o altri corpi celesti

    Lettura SEQUENZIALE e risultati su file
    esterno (.RIS)
               }

  Uses Crt,PASLOGOT;
    Label Etic1,Etic2,Etic3,Etic4,Skip;  {... nell'ordine 170,180,530,200 }
    Const
         rad=pi/180;
         et = 0.01670862; lp = 102.937348*rad; at = 1.000001018;
   Var
      fl,ft                    : Text;
      flag,np,i      : Integer;
      af,xf,yf,zf,xtf,ytf   : Array[1..2] of Real;
      NumOrdin        : String[5];
      NumCat          : String[5];
      nomeastero      : String[27];
      om,lnodo,incl,ecc,asem,
      R11,R12,R21,R22,R31,R32,
      te,el,p,pt,te0,el0,te1,el1,ddte,ddel,DIMOT

         : Real;

      Function Q(te:Real): Real;
   Begin  Q:= 1+ecc*COS(te); End;

      Function T(el:Real): Real;
   Begin  T:= 1+et*COS(el-lp);  End;

      Function x(te:Real): Real;
   Begin  x:= p/Q(te)*(R11*COS(te)+R12*SIN(te));  End;

      Function y(te:Real): Real;
   Begin  y:= p/Q(te)*(R21*COS(te)+R22*SIN(te));  End;

      Function z(te:Real): Real;
   Begin  z:= p/Q(te)*(R31*COS(te)+R32*SIN(te));  End;

      Function xt(el:Real): Real;
   Begin  xt:= pt*COS(el)/T(el);  End;
      Function yt(el:Real): Real;
   Begin  yt:= pt*SIN(el)/T(el);  End;

      Function A(te:Real): Real;
   Begin  A:=p/SQR(Q(te))*(-R11*SIN(te)+R12*(COS(te)+ecc));  End;
      Function B(te:Real): Real;
   Begin  B:=p/SQR(Q(te))*(-R21*SIN(te)+R22*(COS(te)+ecc));  End;
      Function C(te:Real): Real;
   Begin  C:=p/SQR(Q(te))*(-R31*SIN(te)+R32*(COS(te)+ecc));  End;
      Function D(el:Real): Real;
   Begin  D:=pt/SQR(T(el))*(et*(-SIN(el)*COS(el-lp)+COS(el)*SIN(el-lp))-SIN(el));  End;
      Function E(el:Real): Real;
   Begin  E:=pt/SQR(T(el))*(et*(COS(el)*COS(el-lp)+SIN(el)*SIN(el-lp))+COS(el));  End;

      Function Atil(te:Real): Real;
   Begin  Atil:=SQR(Q(te))*A(te)/p;  End;
      Function Btil(te:Real): Real;
   Begin  Btil:=SQR(Q(te))*B(te)/p;  End;
      Function Ctil(te:Real): Real;
   Begin  Ctil:=SQR(Q(te))*C(te)/p;  End;

      Function Aast(te:Real): Real;
   Begin  Aast:=-R11*COS(te)-R12*SIN(te);  End;
      Function Bast(te:Real): Real;
   Begin  Bast:=-R21*COS(te)-R22*SIN(te);  End;
      Function Cast(te:Real): Real;
   Begin  Cast:=-R31*COS(te)-R32*SIN(te);  End;

      Function Dtil(el:Real): Real;
   Begin  Dtil:=et*(-SIN(el)*COS(el-lp)+COS(el)*SIN(el-lp))-SIN(el);  End;
      Function Etil(el:Real): Real;
   Begin  Etil:=et*(COS(el)*COS(el-lp)+SIN(el)*SIN(el-lp))+COS(el);  End;

      Function F1(te,el:Real): Real;
   Begin  F1:=Atil(te)*(x(te)-xt(el))+Btil(te)*(y(te)-yt(el))+Ctil(te)*z(te);  End;
      Function F2(te,el:Real): Real;
   Begin  F2:=Dtil(el)*(x(te)-xt(el))+Etil(el)*(y(te)-yt(el));  End;

      Function F11(te,el:Real): Real;
   Begin  F11:=Aast(te)*(x(te)-xt(el))+Atil(te)*A(te)+Bast(te)*(y(te)
               -yt(el))+Btil(te)*B(te)+Cast(te)*z(te)+Ctil(te)*C(te);  End;
      Function F22(te,el:Real): Real;
   Begin  F22:=-COS(el)*(x(te)-xt(el))-Dtil(el)*D(el)-SIN(el)*(y(te)
               -yt(el))-Etil(el)*E(el);  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
  End;
{------------------------------------------------------
        INIZIO Programmma }
   Begin
     ColoriSchermo;
     Icona;
     write('                   N_ro dei Records (Asteroidi) del FILE = ');
     readln (np);
     ColoriSchermo;

     ASSIGN(fl,'PASDIMOT.DAT');
          Reset(fl);

     ASSIGN(ft,'PASDIMOT.RIS');
          Rewrite(ft);

     ColoriSchermo;

      for i:=1 to np do begin
      readln(fl,NumOrdin,NumCat,nomeastero,om,lnodo,incl,ecc,asem);

  gotoxy(31,15); writeln('    Item:',i:3);
   textcolor(12);
  gotoxy(31,16); writeln('Elementi Orbitali');
   textcolor(14);
  gotoxy(15,17);
    writeln('   w= ',om:9:5,'ψ  κ= ',lnodo:9:5,'ψ  i= ',incl:8:5,'ψ');
  gotoxy(23,18);
    writeln('   e= ',ecc:9:7,'   a= ',asem:9:7,' UA');

lnodo:= lnodo*rad; om:= om*rad; incl:= incl*rad;

{  Angoli di EULERO (i, κ, om) --> Per Calcolo MATRICE di Trasformazione }
R11:= COS(lnodo) * COS(om) - SIN(lnodo) * COS(incl) * SIN(om);
R12:= -COS(lnodo) * SIN(om) - SIN(lnodo) * COS(incl) * COS(om);
R21:= SIN(lnodo) * COS(om) + COS(lnodo) * COS(incl) * SIN(om);
R22:= -SIN(lnodo) * SIN(om) + COS(lnodo) * COS(incl) * COS(om);
R31:= SIN(incl) * SIN(om);
R32:= SIN(incl) * COS(om);
{****************************************************************************}
  p:= asem*(1-ecc*ecc);  pt:= at*(1-et*et);
  textcolor (10);
  gotoxy(29,2); writeln(NumCat+'  '+nomeastero);
  textcolor (14);

  flag:=0;
 Etic1:
        if flag=0 then begin te0:=-om; el0:=lnodo; end;
        if flag=1 then begin te0:=-om+pi; el0:=lnodo+pi; end;
       flag:=flag+1;
 Etic2:
     ddte:= F1(te0,el0)/F11(te0,el0);
      te1:= te0 - ddte;

      ddel:= F2(te1,el0)/F22(te1,el0);
      el1:= el0-ddel;

 gotoxy(6,4);  writeln('    λι= ',ddte,'       λl= ',ddel);
  writeln('     Anom_Corpo= ',te1:9:6,' rad   Long_Terra= ',el1:9:6,' rad');

    IF (ABS(ddte) < 1e-5) AND (ABS(ddel) < 1e-5) THEN goto Etic3;
       te0:= te1; el0:= el1;
goto Etic2;

Etic3:
 writeln;
 writeln('  Soluzione:  Anom_Corpo = ',te1:10:6,' rad    Long_Terra = ',el1:10:6,' rad');

 { Distanza Minima o Massima }
  DIMOT:= SQRT(SQR(x(te1)-xt(el1)) + SQR(y(te1)-yt(el1))+ SQR(z(te1)));
  af[flag]:= DIMOT;
  writeln;
  writeln('           Dist. Minima o Massima [',flag:1,']= ',af[flag]:8:6);
   xf[flag]:= x(te1);    xtf[flag]:= xt(el1);
   yf[flag]:= y(te1);    ytf[flag]:= yt(el1);
   zf[flag]:= z(te1);

  IF flag=2 THEN goto Etic4;
goto Etic1;
   { Scelta della Minima Distanza }
Etic4:
   IF af[1] < af[2] THEN DIMOT:= af[1] ELSE DIMOT:= af[2];
   textcolor (10);
   writeln;writeln('    Distanza MINIMA tra le 2 Orbite -> ',DIMOT:8:6,' UA  (DIMOT)');
   writeln(ft,i:4,') ',copy(NumCat,1,4)+'-'+copy(nomeastero,1,16)+' Dimot -> ', DIMOT:8:6,' UA');
   textcolor (14);
  gotoxy(1,19);

  end; {fine ciclo FOR..END}

    Close(fl);
    Close(ft);
 End.
