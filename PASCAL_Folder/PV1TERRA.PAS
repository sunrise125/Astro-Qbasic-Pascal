 Program PV1TERRA; { POSIZIONE e VELOCITA' ECLITTICHE Terra - Equin. della data
                    - Low Accuracy - Cap.30 Meeus, senza Eq. del Centro }
{$N+}
 Uses Crt;

 Const  rad=pi/180; Kgauss=0.01720209895;

 Var
    ddata      : String;
    T,Xdata,Ydata,Zdata,
    VXdata,VYdata,VZdata  :Double;
    ft     : Text;

  Function MODULO(x,y:Double) : Double;
   Var modu : Double;
   Begin
     modu:=(x/y-INT(x/y))*y;
     MODULO:=modu;
  End;

      Function ArcoTANG(num,den:Double): Double;
    Var atan: Double;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;


  Procedure GiornoGiuliano;
 Var m,ha,mn,sc: WORD;
        y: INTEGER;
  jd,g,mm,yy,a,bw,kk,c,d: Double;
  dd1,dd2,dd3,dd4,dd5,dd6,tm1,tm2: STRING;

    BEGIN
gotoxy(20,4); write('ÚÄÄÄÄÄÄÄ Vettore POS_VEL Terra ÄÄÄÄÄÄ¿');
gotoxy(20,5); write('³                                    ³');
gotoxy(20,6); write('³                                    ³');
gotoxy(20,7); write('³                                    ³');
gotoxy(20,8); write('³                                    ³');
gotoxy(20,9); write('³                                    ³');
gotoxy(20,10);write('³                                    ³');
gotoxy(20,11);write('³                                    ³');
gotoxy(20,12);write('ÀÄ Coordinate ECLITTICHE-Equin. DATAÄÙ');

repeat gotoxy(29,5);write('  Giorno = ');readln(g);
      until (g<32) and (g>=0);
repeat gotoxy(29,6);write('    Mese = ');readln(m);
      until m in [1..12];
       gotoxy(29,7);write('    Anno = ');readln(y);
 repeat gotoxy(29,8);write('  Ore TU = ');readln(ha);
      until ha in [0..23];
 repeat gotoxy(29,9);write('  Min    = ');readln(mn);
      until mn in [0..59];
repeat gotoxy(29,10);write('  Sec    = ');readln(sc);
      until sc in [0..59];

       yy:=y; mm:=m;
          IF m <= 2 then
           begin
           y:=y-1; m:=m+12;
           end;
              bw:=0; kk:=yy+mm/100+g/10000;
          IF kk > 1582.101499 then    {Data >= 15.10.1582--->Rif.Gregoriana}
           begin
          a:=int(y/100); bw:=2 - a + int(a/4);
           end;
              c:=int(365.25*y);
              d:=int(30.6001*(m+1));
      jd:=bw+c+d+g+1720994.5+ha/24+mn/1440+sc/86400;
  str(int(jd):7:0,tm1);
  str(int(frac(jd)*1e5+0.5)/1e5:7:5,tm2);
  delete(tm2,1,1);
       { Secoli Giuliani }
  T:=(jd-2451545)/36525;
     str( g:2:0,dd1);
     str(mm:2:0,dd2);
     str(yy:4:0,dd3);
     str(ha:0,dd4);
     str(mn:0,dd5);
     str(sc:0,dd6);
     ddata:=dd1+'/'+dd2+'/'+dd3+' - '+dd4+':'+dd5+':'+dd6+' TU';
     end;

  Procedure CancellaSchermo;
    Begin
    textcolor(14);textbackground(1); clrscr;
    End;

   Procedure FormatData;
   Begin
       gotoxy (27,1); textcolor(10);
       writeln(ddata);textcolor(14);
       writeln(ft,'                       '+ddata);textcolor(14);
   End;

   Procedure Kepler(m,e:Double; VAR v:Double);
   Var anomecc,variaz,tgv: Double;
   Begin
   anomecc:= m;  {Valore iniziale di E assunto pari all'anom.media M}
 repeat
    variaz:= (m + e * SIN(anomecc) - anomecc) / (1 - e * COS(anomecc));
    anomecc:= anomecc + variaz;
 UNTIL ABS(variaz) < 1e-10;
   tgv:= SIN(anomecc / 2)/COS(anomecc / 2) * SQRT((1 + e) / (1 - e));
   v:= MODULO(2 * ArcoTANG(tgv, 1), 2 * pi);
   End;

   Procedure VettorePosVel;

 Const  c1=102.937348*rad; c2=1.7195269*rad; c3=4.5962e-4*rad; c33=499e-9*rad;
        c4=357.529101*rad; c5=35999.0502962*rad; c6=-1.5594e-4*rad; c66=-289e-9*rad;
        c7=0.01670862; c8=-4.2037e-5; c9=-1.236e-7; c99=4e-11;
        c16=1.000001018;   {Semiasse maggiore}

 Var Lp,M,e,R,te,p,pp   :Double;
   Begin
         { Vettore Posizione }
    Lp:=MODULO(c1+c2*T+c3*T*T+c33*T*T*T,2*pi);  {Longitudine Perielio}
     M:=MODULO(c4+c5*T+c6*T*T+c66*T*T*T,2*pi);  {Anomalia Media}
     e:=c7+c8*T+c9*T*T+c99*T*T*T;               {Obliquit… Eclittica}

    Kepler(M,e,te);
     p:=c16*(1-e*e); pp:=sqrt(p);
     R:=p/(1+e*COS(te));
   Xdata:=COS(Lp)*R*COS(te)-SIN(Lp)*R*SIN(te);
   Ydata:=SIN(Lp)*R*COS(te)+COS(Lp)*R*SIN(te);

         { Vettore Velocit… }

   VXdata:=COS(Lp)*Kgauss/pp*(-SIN(te))-SIN(Lp)*Kgauss/pp*(e+COS(te));
   VYdata:=SIN(Lp)*Kgauss/pp*(-SIN(te))+COS(Lp)*Kgauss/pp*(e+COS(te));

   End;

   Procedure VideoRisultati;
   Begin
  writeln('         Vettore POSIZIONE e VELOCITA'' TERRA [UA; UA/giorno]');
  writeln(ft,'         Vettore POSIZIONE e VELOCITA'' TERRA [UA; UA/giorno]');
  writeln(Xdata:13:8,Ydata:13:8,Zdata:13:8,'³',VXdata:12:9,VYdata:13:9,VZdata:13:9);
  writeln(ft,Xdata:13:8,Ydata:13:8,Zdata:13:8,'³',VXdata:12:9,VYdata:13:9,VZdata:13:9);

   VXdata:=VXdata/Kgauss;
   VYdata:=VYdata/Kgauss;

  writeln;
  writeln(ft);
  writeln('         Vettore POSIZIONE e VELOCITA'' TERRA [Unit… Canoniche]');
  writeln(ft,'         Vettore POSIZIONE e VELOCITA'' TERRA [Unit… Canoniche]');

  writeln(Xdata:13:8,Ydata:13:8,Zdata:13:8,'³',VXdata:12:9,VYdata:13:9,VZdata:13:9);
  writeln(ft,Xdata:13:8,Ydata:13:8,Zdata:13:8,'³',VXdata:12:9,VYdata:13:9,VZdata:13:9);
   End;

{---------------- INIZIO PROGRAMMA -------------------}
  Begin
   ASSIGN (ft,'PV1TERRA.RIS');
   Rewrite(ft);

     CancellaSchermo; GiornoGiuliano;
     CancellaSchermo;
     FormatData;
     VettorePosVel;
     VideoRisultati;

  Close (ft);

  End.
