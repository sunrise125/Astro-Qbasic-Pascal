  Program O_OLBE;
                   { Determinazione dei Parametri dell'Orbita PARABOLICA di un
                     Corpo Celeste con il metodo di OLBERS, note 3 Osservazioni}
  Uses  Graph,Crt;

    Const
         np=3;
         kgauss=0.01720209895;
         rad=pi/180;
   Type
        MatrINTER = Array [1..np] of Integer;
         MatrREAL = Array [1..np] of Real;
   Var
      nme   : String[16];
      w     : String;
      eps,si,co,                          
      TG1,TG2,TG3,dt1,dt2,dt3,            
      Alfa1,Alfa2,Alfa3,Delt1,Delt2,Delt3,
      XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3,
      L1x,L1y,L1z,L2x,L2y,L2z,L3x,L3y,L3z,
      B11,B12,B21,B22,B31,B32,            
      lam1,lam2,lam3,mu1,mu2,mu3,         
      nu1,nu2,nu3,O2,                     
      l1,l2,l3,b1,b2,b3,                  
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
     ArcoSIN:=arctan(valore/SQRT(1-valore*valore));
  End;

      Function ArcoCOS(valore:Real) : Real;
   Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(ABS(1-valore*valore)));
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
  End;

  Procedure LogoDatiOLBERS;
   Begin
     textcolor(10);
     gotoxy(13,2); writeln('CALCOLO DI UN''ORBITA PARABOLICA CON IL METODO DI OLBERS');
     gotoxy(13,3); writeln('-------------------------------------------------------');
     textcolor(14);
     gotoxy(28,4); writeln('Nome del File-Dati = ');gotoxy(49,4);
     read (w);
     textcolor(10);
     gotoxy(23,5);writeln('RIEPILOGO Dati delle 3 OSSERVAZIONI');textcolor(14);
  End;

  Procedure SLICE;
   Var
      fl                       : Text;
      XS,YS,ZS,ARs,DEs,giorno  : MatrREAL;
      ARh,ARm,DEg,DEp          : MatrINTER;
      segno,i                  : Integer;

    Begin
     ASSIGN(fl,w+'.dat');
          Reset(fl);
     for
     i:=1 to np do
      Begin
       readln(fl,giorno[i]);
       readln(fl,XS[i],YS[i],ZS[i]);
       readln(fl,ARh[i],ARm[i],ARs[i],DEg[i],DEp[i],DEs[i]);
  
      write(giorno[i]:11:5,XS[i]:14:10,YS[i]:14:10,ZS[i]:14:10,'³');
      write(ARh[i]:2,ARm[i]:3,ARs[i]:6:2,' ',DEg[i]:3,DEp[i]:3,DEs[i]:5:1);
      writeln;
     End;
      readln(fl,nme,eps);
      write('   '+nme,' î=',eps:9:6);

      TG1:=giorno[1];TG2:=giorno[2];TG3:=giorno[3];

       XS1:=XS[1];YS1:=YS[1];ZS1:=ZS[1];
       XS2:=XS[2];YS2:=YS[2];ZS2:=ZS[2];
       XS3:=XS[3];YS3:=YS[3];ZS3:=ZS[3];

    Alfa1:=(ARs[1]/3600+ARm[1]/60+ARh[1])*15;
    Alfa2:=(ARs[2]/3600+ARm[2]/60+ARh[2])*15;
    Alfa3:=(ARs[3]/3600+ARm[3]/60+ARh[3])*15;
       if  DEg[1]  < 0 then segno:=-1 else segno:=1;
    Delt1:=segno*(DEs[1]/3600+DEp[1]/60+ABS(DEg[1]));
       if  DEg[2]  < 0 then segno:=-1 else segno:=1;
    Delt2:=segno*(DEs[2]/3600+DEp[2]/60+ABS(DEg[2]));
       if  DEg[3]  < 0 then segno:=-1 else segno:=1;
    Delt3:=segno*(DEs[3]/3600+DEp[3]/60+ABS(DEg[3]));
      readln;readln;
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
   Begin
    l1:=ArcoTANG(mu1,lam1);
    l2:=ArcoTANG(mu2,lam2);
    l3:=ArcoTANG(mu3,lam3);
    b1:=ArcoSIN(nu1);
    b2:=ArcoSIN(nu2);
    b3:=ArcoSIN(nu3);
   End;

  Procedure IterazioneRo2OLBERS;
  Label Skip;
  Const srad=206265;
  Var
      aux1,aux2,aux3,QQ,          
      Alfa2c,Delt2c,ResAl2,ResDe2,
      passo,
      c1,c2,cab,fi,csi,eta,fteta  : Real;
   Begin
    aux1:=Sin(l2-O2)/Sin(b2)*Cos(b2);
    aux2:=Sin(l1-O2)/Sin(b1)*Cos(b1);
    aux3:=Sin(l3-O2)/Sin(b3)*Cos(b3);

  QQ:= dt2/dt1*nu1/nu3*(aux1-aux2)/(aux3-aux1);

  ro1:=8.1; passo:=0.1;
     {Inizio CICLO Iterativo}
 repeat
  ro1:=ro1-passo; ro3:= QQ*ro1;
     if ro1<=0.1 then goto Skip;
  x1:=ro1*L1x-XS1; y1:=ro1*L1y-YS1; z1:=ro1*L1z-ZS1;
  x3:=ro3*L3x-XS3; y3:=ro3*L3y-YS3; z3:=ro3*L3z-ZS3;
     r1:=SQRT(x1*x1+y1*y1+z1*z1);   r3:=SQRT(x3*x3+y3*y3+z3*z3);

  c1:=SQRT(SQR(x3-x1)+SQR(y3-y1)+SQR(z3-z1)); {1ø valore della corda c13}

  fi:=2*dt3/((r1+r3)*SQRT(r1+r3));
  csi:=ArcoSIN(0.75*SQRT(2)*fi)/3;
  eta:=ArcoSIN(SQRT(2)*Sin(csi));
  fteta:=Sin(2*eta)/fi;

  c2:=2*dt3/SQRT(r1+r3)*fteta;              {2ø valore della corda c13}
          cab:= ABS(c2-c1);
   
   textcolor(10);
   gotoxy(18,11);writeln('  ro1 = ',ro1:8:6,'    ëc = ',c2-c1:9:7); textcolor(14);
   gotoxy(58,11);writeln('PASSO=',passo:8:6);
  
     if c2-c1 > 0 then
         begin
         ro1:=ro1+passo;
          passo:=passo/10;
         end;

 until cab < 0.0000001;
     {Fine CICLO Iterativo}

Skip:
 textcolor(12);writeln;
  gotoxy(32,13);writeln('c-finale=',c1:11:8); textcolor(14);

  repeat until keypressed;
   End;

  Procedure ParametriOrbitaPARABOLICA;
  Label Salta;
  Var
       Costeta,ff,gg,V1,
      xe1,ye1,ze1, Ve1x,Ve1y,Ve1z,
      Incl,Nodo,Om,Tp,
      V1x,V1y,V1z, sigma1, dd, mm,
      hi,hj,hk,h,
      ni,nj,n,
      ei,ej,ek,e                   : Real;

     Begin
{ Problema di GAUSS per il calcolo del Vettore velocit… V1 }
     costeta:=(x1*x3+y1*y3+z1*z3)/(r1*r3);

        {Parametro focale p e distanza perielia q}
     p:=(r1*r3*(1-costeta))/(r1+r3-2*COS(ArcoCOS(costeta)/2)*SQRT(r1*r3));
     q:=p/2;

        {Parametri f,g}
    ff:=1-r3*(1-costeta)/p;  gg:=r1*r3*SIN(ArcoCOS(costeta))/SQRT(p);

{writeln(' r1=',r1:10:8,'    r3=',r3:10:8,'  é=',ArcoCOS(costeta)/rad:10:6,'   ç3=',dt3/kgauss:6:0);
writeln(' ff=',ff:10:8,'    gg=',gg:10:8);
halt;}
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
 x1:=xe1; y1:=ye1; z1:=ze1;
 V1x:=Ve1x; V1y:=Ve1y; V1z:=Ve1z;

 { Componenti Vettori h,n,e }
  hi:=y1*V1z-z1*V1y;   hj:=z1*V1x-x1*V1z;   hk:=x1*V1y-y1*V1x;
            h:=SQRT(hi*hi+hj*hj+hk*hk);
  ni:=-hj; nj:=hi; n:=SQRT(ni*ni+nj*nj);

ei:= (V1*V1-1/r1)*x1-sigma1*V1x;
ej:= (V1*V1-1/r1)*y1-sigma1*V1y;
ek:= (V1*V1-1/r1)*z1-sigma1*V1z;
  e:= SQRT(ei*ei+ej*ej+ ek*ek);

     {  Parametri Orbitali i,ê,ì }
   Incl:=ArcoCOS(hk/h);
   Nodo:=ArcoCOS(ni/n);
         if hi<0 then Nodo:=2*pi-Nodo;
   Om:=ArcoCOS((ni*ei+nj*ej)/(n*e));
        if ek<0 then Om:=2*pi-Om;

       {Matrice di Trasformazione --> Serve per il GRAFICO dell'Orbita}
 B11:= COS(nodo) * COS(om) - SIN(nodo) * COS(incl) * SIN(om);
 B12:= -COS(nodo) * SIN(om) - SIN(nodo) * COS(incl) * COS(om);
 B21:= SIN(nodo) * COS(om) + COS(nodo) * COS(incl) * SIN(om);
 B22:= -SIN(nodo) * SIN(om) + COS(nodo) * COS(incl) * COS(om);
 B31:= SIN(incl) * SIN(om); B32:= SIN(incl) * COS(om);
 textcolor(10);clrscr;
 gotoxy(18,3);writeln('RISULTATI FINALI - Parametri Orbita PARABOLICA'); textcolor(1);
 gotoxy(34,4);writeln(nme);textcolor(14);
 gotoxy(23,6);writeln('Pass.Perielio  Tp = ',Tp:8:6,' giorni da t2');
 gotoxy(23,7);writeln('Dist.Perielio   q = ',q:10:8,' UA');
 writeln;
 gotoxy(23,9);writeln('Inclinazione    i = ',Incl/rad:12:8,'ø');
 gotoxy(23,10);writeln('Long-Nodo-Asc   ê = ',Nodo/rad:12:8,'ø');
 gotoxy(23,11);writeln('Argom-Perielio  ì = ',Om/rad:12:8,'ø');

 gotoxy(22,13);writeln('[Eccentricit…    e = ',e:12:8,'ø]');
  
   readln; repeat until keypressed;
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
     OutTextXY(xc-310,yc+220,'LINEA VERDE=>Orbita SOPRA l''Eclittica;  LINEA ROSSA=>Orbita SOTTO l''Eclittica');
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
       ColoriSchermo;
       LogoDatiOLBERS;
       SLICE;
       TempiCanoniciOLBERS;
       VettoriL123;
       CosDirettoriRo123;
       LongLatECLITTICA;
       IterazioneRo2OLBERS;
       ParametriOrbitaPARABOLICA;
       GRAFICOparabola;
       ColoriSchermo;
     End.
   {FINE Programma}
