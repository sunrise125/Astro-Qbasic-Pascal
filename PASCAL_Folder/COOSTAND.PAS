Program COOSTAND;    { Coordinate Standard per RIDUZIONE ASTROMETRICA
                       di un oggetto celeste }
{$N+}
  Uses Crt;
   Const
         rad=pi/180; RT=1;
   Var
      Xg,Yg,c, aux, aux2, alfag,deltg,a00,d00,
      deltalfa, dpri, deltadelt, Xgb, Ygb,
      a0,d0,alfa,delt: Double;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
  End;

      Function ArcoTANG(num,den:Real): Real;
    Var atan: Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

  Procedure Titolo;
   Begin;
     textcolor(10);
     gotoxy(13,2); writeln('RIDUZIONE ASTROMETRICA - Metodo di Montenbruck-Pfleger');
     gotoxy(13,3); writeln('        Coordinate Standard X,Y             ');
     textcolor(14);
  End;
{-------------------------------------------------------------------------------

   INIZIO Programma}
   Begin
      ColoriSchermo;
      Titolo;

 {------------  Dati: angoli in Gradi -----------}
   a00:=0; d00:=90;
   alfag:=18; deltg:=59;
 {-----------------------------------------------}

      a0:=a00*rad; d0:=d00*rad;
   alfa:=alfag*rad; delt:=deltg*rad;
      c:=COS(d0)*COS(delt)*COS(alfa-a0)+SIN(d0)*SIN(delt);
      Xg:=RT*(COS(delt)*SIN(alfa-a0))/c;
      Yg:=RT*(SIN(d0)*COS(delt)*COS(alfa-a0)-COS(d0)*SIN(delt))/c;

  writeln(' à= ',alfag:9:5,'  à0= ',a00:9:5,'      ë= ',deltg:9:5,'  ë0= ',d00:9:5);
  writeln;
  writeln(' X=',Xg:15:9,'  Y=',Yg:15:9);

     if (a00=0.0) and (d00=90.0) then begin
 {...... vale solo per Proiez.GNOMONICA Polare: à0=0ø, ë0=90ø }
   aux:=RT*sin(pi/2-delt)/cos(pi/2-delt);
    Xg:=aux*SIN(alfa);
    Yg:=aux*COS(alfa);
  writeln;
  writeln(' X=',Xg:15:9,'  Y=',Yg:15:9,'  ( Verifica: Proiez. GNOMONICA Polare)');
    end;

  writeln;
  writeln('                  Metodo Beltrame');
       deltalfa:=alfa-a0;
       aux2:=sin(delt)/COS(deltalfa);
       dpri:=ArcoTaNg(aux2,cos(delt));
       deltadelt:=dpri-d0;
       Ygb:=-RT*SIN(deltadelt)/COS(deltadelt);
       Xgb:=sqrt(RT*RT+Ygb*Ygb)*COS(delt)*sin(deltalfa)/cos(deltalfa);
  writeln(' X=',Xgb:15:9,'  Y=',Ygb:15:9);

     End.
   {FINE Programma}
