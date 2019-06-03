  Program SERP_RT;   { Determinazione del vettore geocentrico di un punto
                      di una METEORA, note, tramite lastre fotografiche,
                      le coordinate sferiche (à1,ë1), (à2,ë2) da 2 stazioni
                      terrestri di rilevamento di date coordinate geografiche
                      (Metodo SERPILLI). Si calcolano anche le coordinate
                      del punto: (Long, Lat, h)
   ------------------------------------------------------------------------
                      Metodo SEMPLIFICATO con RT = raggio sfera locale
   ------------------------------------------------------------------------
 }
{$N+}
  Uses  Crt;

  Const
         rad=pi/180;
         DU=6378.14;  {Unit… canonica di lunghezza per campo gravitaz.Terra}
         ameri=6378.388; emeri=0.081991889; {Ellisse di Hayford, Colombo:p.91
                                             f=1/297 -> eý-1=(1-f)ý -> e=emeri}

  Var
      fl,ft                 : Text;
         as: String;
         wo: String[35];
         ann,mes,gio,ore,min,sec,band: Integer;
         x1,y1,z1, r1,  C1,S1,h1, C2,S2,h2,
         fi1,lo1,fi2,lo2,
         al1A,de1A,al1B,de1B, teta0, AA,CC,DD,
         ro1,ro2,
         x11, y11, z11, Xg1, Yg1, Zg1,   x21, y21, z21, Xg2, Yg2, Zg2,
         r1x,r1y,r1z,  r1xver,r1yver,r1zver, r1xdef,r1ydef,r1zdef,
         RT, xkm,ykm,zkm,
         rr,EE,FF,PP,QQ,nu,GG,tt,Long,Lat,h, aw, ha,hb, htila,htilb,
         da,db,epsa,epsb,h0a,h0b
                         : Double;

   Function ArcoTANG(num,den:Double) : Double;
     Var atan : Double;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
   End;

      Function ArcoSIN(valore:Real): Real;
   Begin
     ArcoSIN:=arctan(valore/SQRT(ABS(1-valore*valore)));
  End;

  Function S(x   : Double) : Integer;
  Begin
    if x < 0 then S:= -1 else S:= 1
  End;

 Procedure Stringa;
  Var
      lw             : Integer;
      kw,bw,jw,hw    : Double;
      ksw,jsw,hsw    : String;
  Begin
    lw:=S(aw);
    aw:=abs(aw);
    kw:=trunc(aw);
    bw:=frac(aw)*60;
    jw:=trunc(bw);
    str(kw:2:0,ksw);
    str(jw:2:0,jsw);
      if band=1 then
         begin
        hw:=trunc(frac(bw)*600+0.5)/10;
        str(hw:4:1,hsw);
        as:=ksw+'ø '+jsw+''' '+hsw+'"';
        end
        else begin
        hw:=trunc(frac(bw)*600+0.5)/10;
        str(hw:4:1,hsw);
        as:=ksw+'ø '+jsw+''' '+hsw+'"';
      end;
    if lw<0 then as:='-'+as;

  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
   End;

  Procedure LogoDati;
  Label Skip77;
  Var
     u7,jd,jq,t,wn1,mat1,mat2,tsidd,
     je,ore1,gio1,ore2,torig, min1,min2
      :Double;
   Begin
   ASSIGN (fl,'SERP_RT.DAT');
   Reset(fl);
   ASSIGN (ft,'SERP_RT.RIS');
   Rewrite(ft);
    readln(fl,ann,mes,gio,ore,min,sec);
    readln(fl,fi1,lo1,fi2,lo2);
    readln(fl,al1A,de1A,al1B,de1B,wo);
   Close(fl);

     textcolor(10);
     gotoxy(5,2); writeln('Vettore geocentrico di una METEORA fotografata da 2 stazioni terrestri,');
     gotoxy(5,3); writeln('note le coordinate (à,ë) di 1 punto omologo  -  Metodo Serpilli');
     gotoxy(5,4); writeln('-----------------------------------------------------------------------');
           writeln(ft,'    Vettore geocentrico di una METEORA fotografata da 2 stazioni terrestri,');
           writeln(ft,'    note le coordinate (à,ë) di 1 punto omologo  -  Metodo Serpilli');
           writeln(ft,'    -----------------------------------------------------------------------');

 { Raggio della sfera terrestre locale -> RT, dall'ellissoide di Bessel}
    RT:=ameri*sqrt((1-emeri*emeri)/(1-emeri*emeri*SQR(SIN((fi1+fi2)/2*rad))));
 writeln(ft,' RT= ',RT:8:3);

  textcolor(12);
 gotoxy(53,4);writeln(' RT=',RT:8:3,' km ');
  textcolor(14);

   {Calcolo di JD}
      u7:= 2-INT(ann/100)+INT(INT(ann/100)/4);
      IF mes < 3 THEN begin
          ann:= ann - 1; mes:= mes + 12;
          end;
      jd:= INT(365.25 * ann) + INT(30.6001 * (mes + 1)) + gio + 1720994.5;
      IF jd >= 2299160.5 THEN jd:= jd + u7;

      je:= jd - 2415020;
      jq:= je + ore/24 + (min + sec/60)/1440; t:= jq / 36525;
      ore1:= ore; min1:= min;
      gio1:= ore1 / 24 + min1 / 1440;
      ore2:= INT(gio1 * 24); min2:= (gio1 * 24 - ore2) * 60;

      IF INT(min2 + 0.5) = 60 THEN begin
              ore2:= ore2 + 1; min2:= 0;
                end;
      IF ore2 < 24 THEN goto Skip77;
      IF( ore2 >= 24) AND (ore2 < 48) THEN begin
         ore2:= ore2 - 24; end
           ELSE ore2:= ore2 - 48;
Skip77:
  { Calcolo del T-sid a Greenwich all'ora fissata }
    wn1:= ore + min/60 + sec/3600; torig:= wn1/876600;
   mat1:= 0.276919398 + 100.0021359 * (t-torig) + 1.075E-06 * SQR(t-torig);
   mat2:= mat1 - INT(mat1);
  tsidd:= 1.002737908 * wn1 + mat2 * 24;
  tsidd:= tsidd - INT(tsidd / 24) * 24;  {in ORE}
  teta0:= tsidd*15;     {in GRADI}

 {Calcolo componenti del raggio vettore r1 }
         x11:=COS(de1A*rad)*COS(al1A*rad);
         y11:=COS(de1A*rad)*SIN(al1A*rad);
         z11:=SIN(de1A*rad);
         Xg1:=RT*COS(fi1*rad)*COS((lo1+teta0)*rad);
         Yg1:=RT*COS(fi1*rad)*SIN((lo1+teta0)*rad);
         Zg1:=RT*SIN(fi1*rad);
     {------------------------------------------------------}
         x21:=COS(de1B*rad)*COS(al1B*rad);
         y21:=COS(de1B*rad)*SIN(al1B*rad);
         z21:=SIN(de1B*rad);
         Xg2:=RT*COS(fi2*rad)*COS((lo2+teta0)*rad);
         Yg2:=RT*COS(fi2*rad)*SIN((lo2+teta0)*rad);
         Zg2:=RT*SIN(fi2*rad);

{------------------- PUNTO ------------------------------------------------}
    {Parametri AA,CC,DD per il calcolo di (ro1,ro2), che sono le
     distanze topoc. del Punto della meteora dalle staz. A e B }

   AA:= -(x11*x21+y11*y21+z11*z21);
   CC:= (Xg2-Xg1)*x11+(Yg2-Yg1)*y11+(Zg2-Zg1)*z11;
   DD:= (Xg2-Xg1)*(-x21)+(Yg2-Yg1)*(-y21)+(Zg2-Zg1)*(-z21);
      ro1:=(CC-DD*AA)/(1-AA*AA);
      ro2:=(DD-CC*AA)/(1-AA*AA);
   r1x:=ro1*x11+Xg1;  r1xver:=ro2*x21+Xg2;
   r1y:=ro1*y11+Yg1;  r1yver:=ro2*y21+Yg2;
   r1z:=ro1*z11+Zg1;  r1zver:=ro2*z21+Zg2;
       r1xdef:=(r1x/ro1+r1xver/ro2)/(1/ro1+1/ro2);
       r1ydef:=(r1y/ro1+r1yver/ro2)/(1/ro1+1/ro2);
       r1zdef:=(r1z/ro1+r1zver/ro2)/(1/ro1+1/ro2);

 { Raggio geocentrico (r1) in UNITA' CANONICHE: 1ùDU= 6378.14 km }
     x1:=r1xdef/DU;   y1:=r1ydef/DU;   z1:=r1zdef/DU;
     r1:=SQRT(SQR(x1)+SQR(y1)+SQR(z1));

 { Calcolo delle Coordinate GEOGRAFICHE }
     Long:=ArcoTANG(y1,x1)-teta0*rad;
           if Long<0 then Long:=Long+2*pi;
      Lat:=ArcoSIN(z1/r1);
        h:=r1*DU-RT;

   End;


  Procedure StampaRisultati;
  Var num,den,Hor,
      AzSWa,hanga,AzNEa,
      AzSWb,hangb,AzNEb, hkma,hkmb,deltakm
                          : Double;
   Begin
     writeln('                      Sistema EQUATORIALE Geocentrico');
  textcolor(15); 
     writeln('     Vett_Geoc.   r=',r1:13:10,' DU =',r1*DU:6:0,' km');
  textcolor(14);
     writeln('     x= ',x1:13:10,'     y= ',y1:13:10,'     z= ',z1:13:10,' DU');
     writeln('        ',x1*DU:5:0,' km             ',y1*DU:5:0,' km             ',z1*DU:5:0,' km');
     writeln;
     textcolor(12);
     writeln('                      Coordinate GEOGRAFICHE'+wo);
     textcolor(14);
     writeln('     Long= ',Long/rad:7:4,'ø        Lat= ',Lat/rad:7:4,'ø      h= ',h:7:3,' km');

     writeln(ft,'                      Sistema EQUATORIALE Geocentrico');
     writeln(ft,'     Vett_Geoc.   r=',r1:13:10,' DU =',r1*DU:6:0,' km');
     writeln(ft,'     x= ',x1:13:10,'     y= ',y1:13:10,'     z= ',z1:13:10,' DU');
     writeln(ft,'        ',x1*DU:5:0,' km             ',y1*DU:5:0,' km             ',z1*DU:5:0,' km');
     writeln(ft);
     writeln(ft,'                      Coordinate GEOGRAFICHE'+wo);
     writeln(ft,'     Long= ',Long/rad:7:4,'ø        Lat= ',Lat/rad:7:4,'ø      h= ',h:7:3,' km');

 textcolor(15); gotoxy(6,12);
band:=1; aw:=Long/rad; Stringa;  write('      ',as); write(ft,'           ',as);
band:=2; aw:=Lat/rad;  Stringa;  writeln('        ',as); writeln(ft,'        ',as);
      ha:=ArcoSIN(h/ro1);  da:=ro1*COS(ha);
      hb:=ArcoSIN(h/ro2);  db:=ro2*COS(hb);
    epsa:=arctan(da/RT); epsb:=arctan(db/RT);
     h0a:=RT*(1/COS(epsa)-1);    h0b:=RT*(1/COS(epsb)-1);
     htila:=(h-h0A)*COS(epsa); htilb:=(h-h0B)*COS(epsb);


   Close(ft);
    End;
{-----------------------------------------------------------------------------------------------
   INIZIO Programma}

   Begin
       ColoriSchermo;
       LogoDati;
       StampaRisultati;
       End.
   {FINE Programma}
