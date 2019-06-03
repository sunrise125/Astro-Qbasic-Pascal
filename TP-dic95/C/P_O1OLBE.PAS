  Program O1OLBE;
                   { Determinazione dei Parametri dell'Orbita PARABOLICA di un
                     Corpo Celeste con il metodo di OLBERS, note 3 Osservazioni }
  Uses  Graph,Crt;

    Const
         kgauss=0.01720209895;
         rad=pi/180;
   Var
      nme   : String[16];
      w     : String;
      ft  : Text;
      eps,si,co,                          
      TG1,TG2,TG3,dt1,dt2,dt3,            
      Alfa1,Alfa2,Alfa3,Delt1,Delt2,Delt3,
      XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3,
      L1x,L1y,L1z,L2x,L2y,L2z,L3x,L3y,L3z,
      B11,B12,B21,B22,B31,B32,            
      lam1,lam2,lam3,mu1,mu2,mu3,         
      nu1,nu2,nu3,O2,
      l1,l2,l3,b1,b2,b3, QQ, Costeta, ff,gg,V1,
      Incl,Nodo,Om,Tp,
      V1x,V1y,V1z, sigma1, dd, mm,
      hi,hj,hk,h, ni,nj,n, ei,ej,ek,e,
      x1,y1,z1,x2,y2,z2,x3,y3,z3,ro1,ro2,ro3,
      r1,r2,r3,q,p      : Real;

      Function ArcoTANG(num,den:Real) : Real;
   Var atan : Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

      Function ArcoSIN(valore:Real) : Real;
   Begin
     ArcoSIN:=arctan(valore/SQRT(ABS(1-valore*valore)));
  End;

      Function ArcoCOS(valore:Real) : Real;
   Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(ABS(1-valore*valore)));
  End;

   Procedure Titolo;
   Begin
     textcolor(10);
     gotoxy(13,2); writeln(ft,'CALCOLO DI UN''ORBITA PARABOLICA CON IL METODO DI OLBERS');
     gotoxy(13,3); writeln(ft,'-------------------------------------------------------');
     textcolor(14);
   End;
  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
  End;

  Procedure Pausa;
   Begin
     textcolor(12);
     writeln;
     writeln('                               Premere un tasto');
     repeat until keypressed;
     textcolor(14);
  End;

  Procedure InputDatiTastiera;
   Begin
{
 textcolor(12);
 gotoxy(32,4); writeln('Introduzione Dati'); textcolor(14);
 gotoxy( 8,5); write('1^ Osservazione: Giorno = '); readln(TG1);
 gotoxy(25,6); write('Ascensione Retta [Gradi] = '); readln(ALFA1);
 gotoxy(25,7); write('Declinazione     [Gradi] = '); readln(DELT1);

 gotoxy( 8,8);  write('2^ Osservazione: Giorno = '); readln(TG2);
 gotoxy(25,9);  write('Ascensione Retta [Gradi] = '); readln(ALFA2);
 gotoxy(25,10); write('Declinazione     [Gradi] = '); readln(DELT2);

 gotoxy( 8,11);  write('3^ Osservazione: Giorno = '); readln(TG3);
 gotoxy(25,12);  write('Ascensione Retta [Gradi] = '); readln(ALFA3);
 gotoxy(25,13);  write('Declinazione     [Gradi] = '); readln(DELT3);

 gotoxy(8,15);  writeln('Coordinate Equatoriali del SOLE alle Date delle 3 Osservazioni');
 gotoxy(13,16); writeln('1^ Data                 2^ Data               3^ Data');
 gotoxy(8 ,17); write('XS1 = '); readln(XS1);
 gotoxy(8 ,18); write('YS1 = '); readln(YS1);
 gotoxy(8 ,19); write('ZS1 = '); readln(ZS1);
 gotoxy(32,17); write('XS2 = '); readln(XS2);
 gotoxy(32,18); write('YS2 = '); readln(YS2);
 gotoxy(32,19); write('ZS2 = '); readln(ZS2);
 gotoxy(54,17); write('XS3 = '); readln(XS3);
 gotoxy(54,18); write('YS3 = '); readln(YS3);
 gotoxy(54,19); write('ZS3 = '); readln(ZS3);
 gotoxy(28,20); write('Obliquit… Eclittica î = '); readln(EPS);
                 C1.dat }
   TG1:= 1; TG2:= 11; TG3:= 21;
   ALFA1:= 251.435; ALFA2:= 260.5675; ALFA3:= 269.2625;
   DELT1:= 39.513333; DELT2:= 39.86; DELT3:= 40.046666;
   XS1:=0.157632; YS1:= -0.890519; ZS1:= -0.386111;
   XS2:=0.326970; YS2:= -0.850944; ZS2:= -0.368949;
   XS3:=0.486062; YS3:= -0.784961; ZS3:= -0.340347;
   EPS:= 23.4457889;
{                C6.dat
   TG1:= 1; TG2:= 11; TG3:= 21;
   ALFA1:= 199.3025; ALFA2:= 191.905; ALFA3:= 182.1875;
   DELT1:= -40.201666; DELT2:= -39.24; DELT3:= -36.986666;
   XS1:=0.162009; YS1:= -0.889818; ZS1:= -0.385807;
   XS2:=0.331089; YS2:= -0.849624; ZS2:= -0.368385;
   XS3:=0.489959; YS3:= -0.782997; ZS3:= -0.339493;
   EPS:= 23.4457889;
                C4.dat
   TG1:= 5; TG2:= 9.5; TG3:= 14;
   ALFA1:= 336.661916; ALFA2:= 337.798444; ALFA3:= 339.054527;
   DELT1:= 16.50225; DELT2:= 21.291888; DELT3:= 26.407638;
   XS1:=0.967846; YS1:= 0.232994; ZS1:= 0.101065;
   XS2:=0.946596; YS2:= 0.301141; ZS2:= 0.130626;
   XS3:=0.919738; YS3:= 0.367508; ZS3:= 0.159412;
   EPS:= 23.4490407;

                C8.dat
   TG1:= 1; TG2:= 6; TG3:= 11;
   ALFA1:= 323.4475001; ALFA2:= 324.945; ALFA3:= 326.4475001;
   DELT1:= -16.07; DELT2:= -15.4383333; DELT3:= -14.785;
   XS1:=0.9805379; YS1:= 0.1767759; ZS1:= 0.0766467;
   XS2:=0.9617465; YS2:= 0.2537562; ZS2:= 0.1100241;
   XS3:=0.9358601; YS3:= 0.328865; ZS3:= 0.1425898;
   EPS:= 23.4457889;
}
  End;

  Procedure Video1;
   Begin
     gotoxy(10, 4); writeln(ft,'Vettore L: L1x = ',L1x:9:6,'   L1y = ',L1y:9:6,'   L1z = ',L1z:9:6);
     gotoxy(10, 5); writeln(ft,'           L2x = ',L2x:9:6,'   L2y = ',L2y:9:6,'   L2z = ',L2z:9:6);
     gotoxy(10, 6); writeln(ft,'           L3x = ',L3x:9:6,'   L3y = ',L3y:9:6,'   L3z = ',L3z:9:6);

     gotoxy(10, 8); writeln(ft,' Coseni   lam1 = ',lam1:9:6,'   mu1 = ',mu1:9:6,'   nu1 = ',nu1:9:6);
     gotoxy(10, 9); writeln(ft,'Direttori:lam2 = ',lam2:9:6,'   mu2 = ',mu2:9:6,'   nu2 = ',nu2:9:6);
     gotoxy(10,10); writeln(ft,'          lam3 = ',lam3:9:6,'   mu3 = ',mu3:9:6,'   nu3 = ',nu3:9:6);

     gotoxy(7,12); writeln(ft,'                 Longitudine SOLE (t=t2) = ',O2/rad:10:6,'ø');
     gotoxy(7,13); writeln(ft,'LON_Eclit.: l1 = ',  l1/rad:10:6,'ø    l2 = ', l2/rad:10:6,'ø    l3 = ', l3/rad:10:6,'ø');
     gotoxy(7,14); writeln(ft,'LAT_Eclit.: b1 = ',  b1/rad:10:6,'ø    b2 = ', b2/rad:10:6,'ø    b3 = ', b3/rad:10:6,'ø');

     gotoxy(7,16); writeln(ft,'T.Canonici: ç1 = ', dt1:10:8,'     ç2 = ', dt2:10:8,'     ç3 = ', dt3:10:8);
     gotoxy(7,17); writeln(ft,'            Parametro adimensionale Q = ',QQ:10:8);
     gotoxy(9,19); writeln(ft,'Calcolo, per successive approssimazioni, della Corda c = (P1-P3)');
     gotoxy(9,20); writeln(ft,'----------------------------------------------------------------');

   readln;
   End;

  Procedure Video2;
   Begin
     gotoxy(6, 4); writeln(ft,'Vettore r1 : x1 = ',x1:11:8,'   y1 = ',y1:11:8,'   z1 = ',z1:11:8);
     gotoxy(6, 5); writeln(ft,'Vettore r3 : x3 = ',x3:11:8,'   y3 = ',y3:11:8,'   z3 = ',z3:11:8);
     gotoxy(6, 6); writeln(ft,'             r1 = ',r1:11:8,'   r3 = ',r3:11:8,'    é = ',ArcoCOS(costeta)/rad:9:6,'ø');
     gotoxy(6, 7); writeln(ft,'Distanze Geocentr: rho-1 = ',ro1:11:8,'        rho-3 = ',ro3:11:8);
    gotoxy(24, 9); writeln(ft,'Parametro focale p = ',p:11:8);
     gotoxy(6,11); writeln(ft,'Param.(f,g):     f = ',ff:11:8,'     g = ',gg:11:8);
     gotoxy(6,12); writeln(ft,'Vettore  V1:   Vx1 = ',V1x:11:8,'  Vy1 = ',V1y:11:8,'  Vz1 = ',V1z:11:8);
     gotoxy(6,13); writeln(ft,'                V1 = ',V1:11:8);
     gotoxy(6,14); writeln(ft,'Prod.Scalare:   å1 = ',sigma1:11:8);
     gotoxy(6,16); writeln(ft,'Vettore   h:    hi = ',hi:11:8,'   hj = ',hj:11:8,'   hk = ',hk:11:8);
     gotoxy(6,17); writeln(ft,'                 h = ',h:11:8);
     gotoxy(6,18); writeln(ft,'Vettore   n:    ni = ',ni:11:8,'   nj = ',nj:11:8);
     gotoxy(6,19); writeln(ft,'                 n = ',n:11:8);
     gotoxy(6,20); writeln(ft,'Vettore   e:    ei = ',ei:11:8,'   ej = ',ej:11:8,'   ek = ',ek:11:8);
     gotoxy(6,21); writeln(ft,'                 e = ',e:11:8);

 readln;
   End;

  Procedure TempiCanoniciOLBERS;
   Begin
      dt1:=(TG2-TG1)*kgauss;dt2:=(TG3-TG2)*kgauss; dt3:=(TG3-TG1)*kgauss;
   End;


  Procedure VettoriL123;
   Begin
L1x:= COS(Delt1*rad)*COS(Alfa1*rad); L1y:= COS(Delt1*rad)*SIN(Alfa1*rad); L1z:= SIN(Delt1*rad);
L2x:= COS(Delt2*rad)*COS(Alfa2*rad); L2y:= COS(Delt2*rad)*SIN(Alfa2*rad); L2z:= SIN(Delt2*rad);
L3x:= COS(Delt3*rad)*COS(Alfa3*rad); L3y:= COS(Delt3*rad)*SIN(Alfa3*rad); L3z:= SIN(Delt3*rad);
   End;

  Procedure CosDirettoriRo123;
   Begin
   si:=Sin(eps*rad); co:=Cos(eps*rad);{ Sin(î) & Cos(î) }
   lam1:=L1x; lam2:=L2x; lam3:=L3x;
    mu1:=L1z*si+L1y*co; mu2:=L2z*si+L2y*co; mu3:=L3z*si+L3y*co;
    nu1:=L1z*co-L1y*si; nu2:=L2z*co-L2y*si; nu3:=L3z*co-L3y*si;
   O2:=ArcoTANG(YS2,co*XS2);
   End;

  Procedure LongLatECLITTICA;
  Var
      aux1,aux2,aux3:  Real;
   Begin
    l1:=ArcoTANG(mu1,lam1);
    l2:=ArcoTANG(mu2,lam2);
    l3:=ArcoTANG(mu3,lam3);
    b1:=ArcoSIN(nu1);
    b2:=ArcoSIN(nu2);
    b3:=ArcoSIN(nu3);
    aux1:=Sin(l2-O2)/Sin(b2)*Cos(b2);
    aux2:=Sin(l1-O2)/Sin(b1)*Cos(b1);
    aux3:=Sin(l3-O2)/Sin(b3)*Cos(b3);
  QQ:= dt2/dt1*nu1/nu3*(aux1-aux2)/(aux3-aux1);
   End;

  Procedure IterazioneRo2OLBERS;
  Label Skip;
  Var
      passo,
      c1,c2,cab,fi,csi,eta,fteta  : Real;
   Begin

  ro1:=8.1; passo:=0.1;
     {Inizio CICLO Iterativo}
 repeat
  ro1:=ro1-passo; ro3:= QQ*ro1;
     if ro1<=0.1 then goto Skip;
  x1:=ro1*L1x-XS1; y1:=ro1*L1y-YS1; z1:=ro1*L1z-ZS1;
  x3:=ro3*L3x-XS3; y3:=ro3*L3y-YS3; z3:=ro3*L3z-ZS3;
     r1:=SQRT(x1*x1+y1*y1+z1*z1);   r3:=SQRT(x3*x3+y3*y3+z3*z3);

  c1:=SQRT(SQR(x3-x1)+SQR(y3-y1)+SQR(z3-z1)); {Corda ®geometrica¯ P1-P3}

  fi:=2*dt3/((r1+r3)*SQRT(r1+r3));
  csi:=ArcoSIN(0.75*SQRT(2)*fi)/3;
  eta:=ArcoSIN(SQRT(2)*Sin(csi));
  fteta:=Sin(2*eta)/fi;

  c2:=2*dt3/SQRT(r1+r3)*fteta;              {Corda ®dinamica¯ P1-P3}
          cab:= ABS(c2-c1);
   
   textcolor(10);
   gotoxy(13,21);writeln(ft,'  rho-1 = ',ro1:8:6,'    ëc = ',c2-c1:9:7,'    PASSO=',passo:8:6); textcolor(14);
  
     if c2-c1 > 0 then
         begin
         ro1:=ro1+passo;
          passo:=passo/10;
         end;

 until cab < 0.0000001;
     {Fine CICLO Iterativo}

Skip:
 textcolor(12);writeln(ft);
  gotoxy(30,22);writeln(ft,'Corda c =',c1:11:8); textcolor(14);

  readln;
   End;

  Procedure ParametriOrbitaPARABOLICA;
  Label Salta;
  Var
      xe1,ye1,ze1, Ve1x,Ve1y,Ve1z,
      xx1,yy1,zz1, V1xx,V1yy,V1zz  : Real;

     Begin
{ Problema di GAUSS per il calcolo del Vettore velocit… V1 }
     costeta:=(x1*x3+y1*y3+z1*z3)/(r1*r3);

        {Parametro focale p e distanza perielia q}
     p:=(r1*r3*(1-costeta))/(r1+r3-2*COS(ArcoCOS(costeta)/2)*SQRT(r1*r3));
     q:=p/2;

        {Parametri f,g}
    ff:=1-r3*(1-costeta)/p;  gg:=r1*r3*SIN(ArcoCOS(costeta))/SQRT(p);

   V1x:=(x3-ff*x1)/gg;
   V1y:=(y3-ff*y1)/gg;
   V1z:=(z3-ff*z1)/gg;
        V1:= SQRT(V1x*V1x+V1y*V1y+V1z*V1z);

    { å1 = prodotto scalare vettori r1,V1 }
   sigma1:=x1*V1x+y1*V1y+z1*V1z;

      { Tp = Epoca passaggio al perielio }
     dd:=sigma1;
     mm:=0.5*(p*dd+dd*dd*dd/3);
  Tp:= -(mm+dt1)/kgauss;

  {Passaggio da coord.Equatoriali a Eclittiche}
        xe1:=x1; ye1:=y1*co+z1*si; ze1:=-y1*si+z1*co;
       Ve1x:=V1x; Ve1y:=V1y*co+V1z*si; Ve1z:=-V1y*si+V1z*co;
 xx1:=xe1; yy1:=ye1; zz1:=ze1;
 V1xx:=Ve1x; V1yy:=Ve1y; V1zz:=Ve1z;

 { Componenti Vettori h,n,e }
  hi:=yy1*V1zz-zz1*V1yy;   hj:=zz1*V1xx-xx1*V1zz;   hk:=xx1*V1yy-yy1*V1xx;
            h:=SQRT(hi*hi+hj*hj+hk*hk);
  ni:=-hj; nj:=hi; n:=SQRT(ni*ni+nj*nj);

ei:= (V1*V1-1/r1)*xx1-sigma1*V1xx;
ej:= (V1*V1-1/r1)*yy1-sigma1*V1yy;
ek:= (V1*V1-1/r1)*zz1-sigma1*V1zz;
  e:= SQRT(ei*ei+ej*ej+ ek*ek);

     {  Parametri Orbitali i,ê,ì }
   Incl:=ArcoCOS(hk/h);
   Nodo:=ArcoCOS(ni/n);
         if hi<0 then Nodo:=2*pi-Nodo;
   Om:=ArcoCOS((ni*ei+nj*ej)/(n*e));
        if ek<0 then Om:=2*pi-Om;

    ColoriSchermo;Titolo; Video2;

       {Matrice di Trasformazione --> Serve per il GRAFICO dell'Orbita}
 B11:= COS(nodo) * COS(om) - SIN(nodo) * COS(incl) * SIN(om);
 B12:= -COS(nodo) * SIN(om) - SIN(nodo) * COS(incl) * COS(om);
 B21:= SIN(nodo) * COS(om) + COS(nodo) * COS(incl) * SIN(om);
 B22:= -SIN(nodo) * SIN(om) + COS(nodo) * COS(incl) * COS(om);
 B31:= SIN(incl) * SIN(om); B32:= SIN(incl) * COS(om);
 textcolor(10);clrscr;
 gotoxy(18,3);writeln(ft,'RISULTATI FINALI - Parametri Orbita PARABOLICA'); textcolor(1);
 gotoxy(34,4);writeln(nme);textcolor(14);
 gotoxy(23,6);writeln(ft,'Pass.Perielio  Tp = ',Tp:8:6,' giorni da t2');
 gotoxy(23,7);writeln(ft,'Dist.Perielio   q = ',q:10:8,' UA');
 writeln(ft);
 gotoxy(23,9);writeln(ft,'Inclinazione    i = ',Incl/rad:12:8,'ø');
 gotoxy(23,10);writeln(ft,'Long-Nodo-Asc   ê = ',Nodo/rad:12:8,'ø');
 gotoxy(23,11);writeln(ft,'Argom-Perielio  w = ',Om/rad:12:8,'ø');

 gotoxy(22,13);writeln(ft,'[Eccentricit…    e = ',e:12:8,'ø]');
  
   readln;
    End;

  Procedure GRAFICOparabola;
 Label vai,SKIP;
  Var
      ch                      : String;
      gd,gm,k                 : Integer;
      xc,yc,x,y,colore        : Word;
      ka,xorbp,yorbp,ro,
      xeclp,yeclp,zeclp,rsole : Real;
   Begin
    gd:=detect;initgraph(gd,gm,'');
    if graphresult <> grOk then halt(1);
       ka:=30;
vai:
    SetFillStyle(1,6);      { Campitura [0..11] e Colore SFONDO }
    FloodFill(0,0,12);
    SetColor(14);

    SetTextStyle(4,0,5);
    OutText('            ORBITA');
        SetFillStyle(1,12);
     xc:=GetmaxX div 2; yc:=GetmaxY div 2;
     SetTextStyle(4,0,3); {Font,Direz.,Size}
     OutTextXY(xc+250,yc-30,'Punto');OutTextXY(xc+240,yc-8,'Gamma');
     SetTextStyle(0,0,1); {Font,Direz.,Size}
     OutTextXY(xc-310,yc+220,'LINEA VERDE=>Orbita SOPRA l`Eclittica;  LINEA ROSSA=>Orbita SOTTO l`Eclittica');
     OutTextXY(xc+176,yc+200,'SOLE Fuori Scala');
    SetLineStyle(0,0,1);
    Rectangle(xc-315,yc+160,xc-75,yc+208);
      SetColor(10);
      OutTextXY(xc-310,yc+165,' + => Ingrandisce     Grafico');
      OutTextXY(xc-310,yc+180,' - => Rimpicciolisce  Grafico');
      OutTextXY(xc-310,yc+195,'USCITA  con Qualsiasi Freccia');
      SetColor(14);
    rsole:=0.3*ka;
    SetLineStyle(3,0,1);
    line(xc,yc,2*xc,yc);
    FillEllipse(xc,yc,round(rsole),round(rsole));
       {PLOT Orbita}
    colore:=14;
      for k:=0 to 360 do
      Begin
          if (k>175) and (k<185) then goto SKIP;
               ro:=ka*2*q/(1+cos(k*rad));
      XORBP:= ro * COS(k*rad); YORBP:= ro * SIN(k*rad);{Coordinate ORBITALI

      Coordinate ECLITTICHE di P }
 XECLP:= B11 * XORBP + B12 * YORBP;
 YECLP:= B21 * XORBP + B22 * YORBP;
 ZECLP:= B31 * XORBP + B32 * YORBP;
            x:=xc+round(xeclp);  y:=yc-round(yeclp);
              if zeclp>0 then PutPixel(x,y,10) ELSE PutPixel(x,y,12);
SKIP:
         end;
  repeat
     ch:=readkey;
       if ch='+' then
         begin
         ClearDevice;
          ka:=ka*2;
          goto vai;
          end;

       if ch='-' then
         begin
         ClearDevice;
          ka:=ka/2;
          goto vai;
         end;
  until keypressed;

    closegraph;
    restoreCrtMode;
  End;
{-----------------------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
     Assign(ft,'CC1.RIS');
     Rewrite(ft);
       ColoriSchermo; Titolo;
       InputDatiTastiera; Pausa; ColoriSchermo;
       TempiCanoniciOLBERS;
       VettoriL123;
       CosDirettoriRo123;
       LongLatECLITTICA;
    Titolo; Video1;
       IterazioneRo2OLBERS;
       ParametriOrbitaPARABOLICA;
       GRAFICOparabola;
       ColoriSchermo;
     Close(ft);
     End.
   {FINE Programma}
