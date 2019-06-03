 Program PV_TERRA; { POSIZIONE e VELOCITA' ECLITTICHE Terra - Equin. della data
                    - Low Accuracy - Cap.24 Meeus }
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

   Procedure VettorePosVel;

 Const  c1=100.46645*rad; c2=36000.76983*rad; c3=3.032e-4*rad;
        c4=357.52910*rad; c5=35999.05030*rad; c6=-1.559e-4*rad;
        c7=0.016708617; c8=-4.2037e-5; c9=-1.236e-7;
        c10=1.9146*rad; c11=-0.004817*rad; c12=-1.4e-5*rad;
        c13=0.019993*rad; c14=-1.01e-4*rad; c15=0.00029*rad;
        c16=1.000001018;

 Var L0,M,C,e,R,te,v, Xpri,Ypri, Rpri,tepri, epri,vpri, Cpri,
     num,den                     :Double;
   Begin
         { Vettore Posizione }
    L0:=MODULO(c1+c2*T+c3*T*T,2*pi);
     M:=MODULO(c4+c5*T+c6*T*T,2*pi);
     e:=c7+c8*T+c9*T*T;
  C:=(c10+c11*T+c12*T*T)*SIN(M)+(c13+c14*T)*SIN(2*M)+c15*SIN(3*M);
     te:=L0+C; v:=M+C;
     R:=c16*(1-e*e)/(1+e*COS(v));
   Xdata:=R*COS(te);
   Ydata:=R*SIN(te);

         { Vettore Velocit… }
   Cpri:=(c11+2*c12*T)*SIN(M)+COS(M)*(c5+2*c6*T)*(c10+c11*T+c12*T*T);
   Cpri:=Cpri+c14*SIN(2*M)+2*COS(2*M)*(c5+2*c6*T)*(c13+c14*T);
   Cpri:=Cpri+c15*3*COS(3*M)*(c5+2*c6*T);
    tepri:=(c2+2*c3*T)+Cpri;
     vpri:=(c5+2*c6*T)+Cpri;
     epri:=c8+2*c9*T;
   num:=c16*(-2*e*epri)*(1+e*COS(v))-(epri*COS(v)-e*SIN(v)*vpri)*c16*(1-e*e);
   den:=SQR(1+e*COS(v));
 Rpri:=num/den;
   VXdata:=(Rpri*COS(te)-R*tepri*SIN(te))/36525;
   VYdata:=(Rpri*SIN(te)+R*tepri*COS(te))/36525;
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
   ASSIGN (ft,'PV_TERRA.RIS');
   Rewrite(ft);

     CancellaSchermo; GiornoGiuliano;
     CancellaSchermo;
     FormatData;
     VettorePosVel;
     VideoRisultati;

  Close (ft);

  End.
