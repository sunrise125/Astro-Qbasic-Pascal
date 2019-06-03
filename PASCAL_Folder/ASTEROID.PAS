Program ASTEROID;   {EFFEMERIDI Asteroidi
                     Dbase = STAMP97, 6678 pianetini catalogati
                     Epoca 1997 12 18.0 }

  Uses  Crt,Graph;
   Label Skip;
   Const LUNG=122; rad=pi/180; kgauss=0.01720209895; uv0=29.78469169;
   Var
      wo: String[24];
      Hmag,Gmag,m,om,lnodo,incl,e,n,a,
      g,mm,yy,
      r11,r21,r31,r12,r22,r32,p,
       Rr11,Rr21,Rr31,Rr12,Rr22,Rr32,pt,e1,xt,yt,zt,
      jd,aw,xeclp,yeclp,zeclp,
      xeclpsave,yeclpsave,xtsave,ytsave,
      Vxeclp,Vyeclp,Vzeclp       : Real;
      fb          : file of Byte;
      xxx,u,band  : Integer;
      k,j,zoz     : LongInt;
      byt         : Byte;
      Us          : String[1];
      ncat,c3,c4  : String[5];
      nome        : String[18];
      c5,c9,c7,c6 : String[9];
      c8          : String[8];
      c10,c11     : String[10];
      as          : String;

      Function MODULO(x,y:Real) : Real;
   Var modu : Real;
   Begin
     modu:=(x/y-INT(x/y))*y;
     MODULO:=modu;
  End;

  Function S(x   : Real) : Integer;
  Begin
    if x < 0 then S := -1 else S := 1
  End;

      Function ArcoTANG(num,den:Real): Real;
    Var atan: Real;
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

      Function ArcoCOS(valore:Real): Real;
   Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(ABS(1-valore*valore)));
  End;

 Procedure Stringa;
  Var
      lw             : Integer;
      kw,bw,jw,hw    : Real;
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
        hw:=trunc(frac(bw)*6000+0.5)/100;
        str(hw:5:2,hsw);
        as:=ksw+'h '+jsw+'m '+hsw+'s';
        end
        else begin
        hw:=trunc(frac(bw)*600+0.5)/10;
        str(hw:4:1,hsw);
        as:=ksw+'ø '+jsw+'` '+hsw+'``';
      end;
    if lw<0 then as:='-'+as;

  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
  End;

  Procedure Titolo;
   Begin;
     textcolor(10);
     gotoxy(17,2); writeln('EFFEMERIDI DI ASTEROIDI CATALOGATI (1..6678)');
     gotoxy(17,3); writeln('---------------- Ast_97.dat ----------------');
     textcolor(14);
  End;

   Procedure Kepler(m,e:Real; VAR v:Real);
   Var anomecc,variaz,tgv: Real;
   Begin
   anomecc:= m;  {Valore iniziale di E assunto pari all'anom.media M}
 repeat
    variaz:= (m + e * SIN(anomecc) - anomecc) / (1 - e * COS(anomecc));
    anomecc:= anomecc + variaz;
 UNTIL ABS(variaz) < 1e-10;
   tgv:= SIN(anomecc / 2)/COS(anomecc / 2) * SQRT((1 + e) / (1 - e));
   v:= MODULO(2 * ArcoTANG(tgv, 1), 2 * pi);
   End;

   Procedure ElementiOrbitali;
   Begin
  textcolor(12); gotoxy(23,9); write('Elementi Orbitali:'+ wo); textcolor(14);
  gotoxy(20,10); write('   a= ',a:10:7,'       e= ',e:9:7);
  gotoxy(11,11); write('   i= ',incl:8:5,'ø     w= ',om:9:5,'ø    ê= ',lnodo:9:5,'ø');
  gotoxy(28,12); write('   M= ',m:9:5,'ø (1997 12 18.0)');
   End;

  Procedure GiornoGiuliano;
 Var m,ha,mn: WORD;
        y: INTEGER;
  a,b,kk,c,d: REAL;

    BEGIN
gotoxy(20,4); write('ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
gotoxy(20,5); write('³                                    ³');
gotoxy(20,6); write('³                                    ³');
gotoxy(20,7); write('³                                    ³');
gotoxy(20,8); write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
repeat gotoxy(22,5);write(' Giorno = ');readln(g);
      until (g<32) and (g>=0);
repeat gotoxy(22,6);write('   Mese = ');readln(m);
      until m in [1..12];
       gotoxy(22,7);write('   Anno = ');readln(y);
repeat gotoxy(43,5);write(' Ore TU = ');readln(ha);
      until ha in [0..23];
repeat gotoxy(43,6);write(' Min TU = ');readln(mn);
      until mn in [0..59];

       yy:=y; mm:=m;
          IF m <= 2 then
           begin
           y:=y-1; m:=m+12;
           end;
              b:=0; kk:=yy+mm/100+g/10000;
          IF kk > 1582.101499 then    {Data >= 15.10.1582--->Rif.Gregoriana}
           begin
          a:=int(y/100); b:=2 - a + int(a/4);
           end;
              c:=int(365.25*y);
              d:=int(30.6001*(m+1));
      jd:=b+c+d+g+1720994.5+ha/24+mn/1440;
     end;

   Procedure CALCOLO;
   Label Skip1;
   Const JDt0=2450800.5; JD2000=2451545.0;
   Var
       card: String[1];
       M1,a1,om1,lnodo1,incl1,
       M2,a2,e2,om2,lnodo2,incl2,
       Mast,Mter,te1,te2,
       x,y,z,Vx,Vy,Vz,dx,dy,dz,deltap,rp,
       Vxt,Vyt,Vzt,
       yeq,zeq,Vyeq,Vzeq,yteq,zteq,Vyteq,Vzteq,CO,SI,
       angfase,tg2,aux1,aux2,gi,me,an,
       magn,magn0,magn1,magn2,
       AR,DE,V,provet,elong,rs
                 : Real;
 Begin
  {Elementi orbitali Asteroide selezionato dal Dbase}
  M2:= m * rad; a2:= a; e2:= e;
 om2:= om * rad; lnodo2:= lnodo*rad; incl2:= incl*rad;

{ Terra (Fonte: Astronomical Algorithms di J.Meeus) }
  a1:= 1.000001018; e1:= 0.01670862;
 om1:= 102.937348 * rad; lnodo1:= 0; incl1:= 0;

{ Elementi della MATRICE di Trasformazione  1)=Terra;  2)=Asteroide }
Rr11:= COS(om1) * COS(lnodo1) - SIN(om1) * COS(incl1) * SIN(lnodo1);
Rr21:= COS(om1) * SIN(lnodo1) + SIN(om1) * COS(incl1) * COS(lnodo1);
Rr31:= SIN(om1) * SIN(incl1);
Rr12:= -SIN(om1) * COS(lnodo1) - COS(om1) * COS(incl1) * SIN(lnodo1);
Rr22:= -SIN(om1) * SIN(lnodo1) + COS(om1) * COS(incl1) * COS(lnodo1);
Rr32:= COS(om1) * SIN(incl1);

r11:= COS(om2) * COS(lnodo2) - SIN(om2) * COS(incl2) * SIN(lnodo2);
r21:= COS(om2) * SIN(lnodo2) + SIN(om2) * COS(incl2) * COS(lnodo2);
r31:= SIN(om2) * SIN(incl2);
r12:= -SIN(om2) * COS(lnodo2) - COS(om2) * COS(incl2) * SIN(lnodo2);
r22:= -SIN(om2) * SIN(lnodo2) + COS(om2) * COS(incl2) * COS(lnodo2);
r32:= COS(om2) * SIN(incl2);

   pt:= a1 * (1 - SQR(e1));  {Parametro orbitale Terra}
    p:= a2 * (1 - SQR(e2));  {   "         "     Asteroide}
 Skip1:
  GiornoGiuliano;
 { Anomalie MEDIE }
    Mast:= M2 + kgauss / (a2 *SQRT(a2)) * (jd-JDt0);
    Mter:= -2.470899 + 0.985609113 * (jd-JD2000);
    Mast:= MODULO(Mast, 2 * pi); Mter:= rad*MODULO(Mter, 360);
    Kepler(Mter,e1,te1);
    Kepler(Mast,e2,te2);
    x:= p / (1 + e2 * COS(te2)) * (r11 * COS(te2) + r12 * SIN(te2));
    y:= p / (1 + e2 * COS(te2)) * (r21 * COS(te2) + r22 * SIN(te2));
    z:= p / (1 + e2 * COS(te2)) * (r31 * COS(te2) + r32 * SIN(te2));
          Vx:= kgauss / SQRT(p) * (-r11 * SIN(te2) + r12 * (COS(te2) + e2));
          Vy:= kgauss / SQRT(p) * (-r21 * SIN(te2) + r22 * (COS(te2) + e2));
          Vz:= kgauss / SQRT(p) * (-r31 * SIN(te2) + r32 * (COS(te2) + e2));

    xeclp:=x;     yeclp:=y;    zeclp:=z;  {Per memorizzazione}
  Vxeclp:=Vx;   Vyeclp:=Vy;  Vzeclp:=Vz;  {Per memorizzazione}

    xt:= pt / (1 + e1 * COS(te1)) * (Rr11 * COS(te1) + Rr12 * SIN(te1));
    yt:= pt / (1 + e1 * COS(te1)) * (Rr21 * COS(te1) + Rr22 * SIN(te1));
    zt:= pt / (1 + e1 * COS(te1)) * (Rr31 * COS(te1) + Rr32 * SIN(te1));
          Vxt:= kgauss / SQRT(pt) * (-Rr11 * SIN(te1) + Rr12 * (COS(te1) + e1));
          Vyt:= kgauss / SQRT(pt) * (-Rr21 * SIN(te1) + Rr22 * (COS(te1) + e1));
          Vzt:= kgauss / SQRT(pt) * (-Rr31 * SIN(te1) + Rr32 * (COS(te1) + e1));

    xeclpsave:=x;  yeclpsave:=y;  {Per il grafico}
    xtsave:=xt;    ytsave:=yt;

{ Cambio di coordinate: da eclittiche ad equatoriali}
  CO:= COS(23.4392911 * rad); SI:= SIN(23.4392911 * rad);
  yeq:= y * CO - z * SI;       yteq:= yt * CO - zt * SI;
  zeq:= y * SI + z * CO;       zteq:= yt * SI + zt * CO;

  Vyeq:= Vy * CO - Vz * SI;    Vyteq:= Vyt * CO - Vzt * SI;
  Vzeq:= Vy * SI + Vz * CO;    Vzteq:= Vyt * SI + Vzt * CO;
        y:=yeq;  z:=zeq;  Vy:=Vyeq;  Vz:=Vzeq;
       yt:=yteq;  zt:=zteq;  Vyt:=Vyteq;  Vzt:=Vzteq;
 {Fine Cambio di coordinate}

    dx:= x - xt; dy:= y - yt; dz:= z - zt;
    deltap:= SQRT(SQR(dx) + SQR(dy) + SQR(dz));
    rp:= SQRT(SQR(x) + SQR(y) + SQR(z));
    rs:= SQRT(SQR(xt) + SQR(yt) + SQR(zt));
     V:= SQRT(SQR(Vx) + SQR(Vy) + SQR(Vz));

     AR:= ArcoTANG(dy, dx)/(15*rad);
     DE:= ArcoSIN(dz/deltap)/rad;

 { Angolo di Fase di P }
   angfase:= ArcoCOS((x*dx+y*dy+z*dz)/(rp*deltap));

 { Elongazione di P }
provet:= -(dx*xt + dy*yt + dz*zt); {Prodotto vettoriale (P-T)ù(S-T)}
 elong:= ArcoCOS(provet/ (rs * deltap));
 IF S(dx*yt - dy*xt) = 1 THEN card:= 'E' ELSE card:= 'W';

{ MAGNITUDINE di P  --> m = Hmag + 5ùlog(rùdelta) - 2.5 *
                          *log[(1-Gmag)*e^(-3.33*tg(á/2)^.63)+
                               +Gmag*e^(-1.87*tg(á/2)^1.22)]
                       --------------------------------------------------
                        dove á=angolo di fase  e log = Logaritmo decimale }
        
         magn0:= Hmag + 5 * LN(rp * deltap) / LN(10);
           tg2:= SIN(angfase/2)/COS(angfase/2);
          aux1:= EXP(0.63*LN(tg2)); aux2:=EXP(1.22*LN(tg2));
         magn1:= (1 - Gmag) * EXP(-3.33 * aux1);
         magn2:= Gmag * EXP(-1.87 * aux2);

   magn:= magn0 - 2.5 * LN(magn1 + magn2) / LN(10);

gotoxy(2,13); write('ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
gotoxy(2,14); write('³ POSIZIONE Eclittica                     VELOCITA'' Eclittica               ³');
gotoxy(2,15); write('³                                                                           ³');
gotoxy(2,16); write('³                                                                           ³');
gotoxy(2,17); write('³                                                                           ³');
gotoxy(2,18); write('³                                                                           ³');
gotoxy(2,19); write('³                                                                           ³');
gotoxy(2,20); write('³                                                                           ³');
gotoxy(2,21); write('³                                                                           ³');
gotoxy(2,22); write('³                                                                           ³');
gotoxy(2,23); write('³                                                                           ³');
gotoxy(2,24); write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');

gotoxy(4,15);
write('x = ',xeclp:10:7,' UA     Vx=',Vxeclp:11:8,' UA/g=',Vxeclp/kgauss:10:7,' UA/UT0=',Vxeclp/kgauss*uv0:8:4,' Km/s');
gotoxy(4,16);
write('y = ',yeclp:10:7,' UA     Vy=',Vyeclp:11:8,' UA/g=',Vyeclp/kgauss:10:7,' UA/UT0=',Vyeclp/kgauss*uv0:8:4,' Km/s');
gotoxy(4,17);
write('z = ',zeclp:10:7,' UA     Vz=',Vzeclp:11:8,' UA/g=',Vzeclp/kgauss:10:7,' UA/UT0=',Vzeclp/kgauss*uv0:8:4,' Km/s');
textcolor(11);
gotoxy(4,18); write('r = ',rp:10:7,' UA     V =',V:11:8,' UA/g=',V/kgauss:10:7,' UA/UT0=',V/kgauss*uv0:8:4,' Km/s');
textcolor(14);
gotoxy(13,19); write('AngFase á =',angfase/rad:8:4,'ø  Elongaz.=',abs(elong)/rad:8:4,' '+card,'   Magnitud.= ',magn:4:1);
gotoxy(33,20); textcolor(6);writeln('Equinozio J2000');
 textcolor(15); gotoxy(23,21);
band:=1; aw:=AR; Stringa;  write('à= ',as);
band:=2; aw:=DE;      Stringa;  writeln('   ë= ',as);
 textcolor(14);
       gotoxy(27,22);
       textcolor(12); writeln('r=',rp:8:5,'   Delta=',deltap:8:5); textcolor(14);
 gotoxy(25,23); textcolor(10); writeln('Grafico,Altra Data,Uscita: (g/d/u)');
 textcolor(14);

  repeat
     Us:=readkey;
   if (Us='d') then goto Skip1;
   if (Us='u') then halt;
   until (Us='g');
   End;

  Procedure GRAFICOconica;
 Label vai;
   Var
      ch,dt1,dt2,dt3    : String;
      gd,gm                : Integer;
      xc,yc,x,y,colore     : Word;
      ka,rsole,k,
      xm,ym,zm,am,em,pm,lnodom,inclm,omm,
      rm11,rm12,rm21,rm22,rm31,rm32,
      xg,yg,zg,ag,eg,pg,lnodog,inclg,omg,
      rg11,rg12,rg21,rg22,rg31,rg32

      : Real;
   Begin
    gd:=detect;initgraph(gd,gm,'');
    if graphresult <> grOk then halt(1);
     ka:=30;

         {Orbita MARTE}
     omm:=286.43045*rad;
     lnodom:=49.57444*rad;
     inclm :=1.85005*rad;
     em:=0.0933765;
     am:=1.5237606;
              pm:=am*(1-em*em);
rm11:= COS(omm) * COS(lnodom) - SIN(omm) * COS(inclm) * SIN(lnodom);
rm21:= COS(omm) * SIN(lnodom) + SIN(omm) * COS(inclm) * COS(lnodom);
rm31:= SIN(omm) * SIN(inclm);
rm12:= -SIN(omm) * COS(lnodom) - COS(omm) * COS(inclm) * SIN(lnodom);
rm22:= -SIN(omm) * SIN(lnodom) + COS(omm) * COS(inclm) * COS(lnodom);
rm32:= COS(omm) * SIN(inclm);
         {Fine Orbita MARTE}

         {Orbita GIOVE}
     omg:=275.24335*rad;
     lnodog:=100.47058*rad;
     inclg :=1.30463*rad;
     eg:=0.0484029;
     ag:=5.2024651;
              pg:=ag*(1-eg*eg);
rg11:= COS(omg) * COS(lnodog) - SIN(omg) * COS(inclg) * SIN(lnodog);
rg21:= COS(omg) * SIN(lnodog) + SIN(omg) * COS(inclg) * COS(lnodog);
rg31:= SIN(omg) * SIN(inclg);
rg12:= -SIN(omg) * COS(lnodog) - COS(omg) * COS(inclg) * SIN(lnodog);
rg22:= -SIN(omg) * SIN(lnodog) + COS(omg) * COS(inclg) * COS(lnodog);
rg32:= COS(omg) * SIN(inclg);
         {Fine Orbita GIOVE}

vai:
    SetFillStyle(1,6);      { Campitura [0..11] e Colore SFONDO }
    FloodFill(0,0,12);
    SetColor(14);
    SetTextStyle(1,0,5);
    OutText(ncat+'-'+nome);
    SetTextStyle(1,0,3);
    xc:=GetmaxX div 2; yc:=GetmaxY div 2;
    str(g:2:0,dt1); str(mm:2:0,dt2); str(yy:4:0,dt3);
    OutTextXY(xc+172,yc-222,dt1+'/'+dt2+'/'+dt3);
        SetFillStyle(1,12);
     SetTextStyle(1,0,3); {Font,Direz.,Size}
     OutTextXY(xc+250,yc-30,'Punto');OutTextXY(xc+232,yc-8,'Gamma');
     SetTextStyle(0,0,1); {Font,Direz.,Size}
     OutTextXY(xc-310,yc+220,'LINEA VERDE=>Orbita SOPRA l''Eclittica;  LINEA ROSSA=>Orbita SOTTO l''Eclittica');
     OutTextXY(xc+58,yc+200,'Altre ORBITE: Terra,Marte,Giove');
    SetLineStyle(0,0,1);
    Rectangle(xc-315,yc+160,xc-75,yc+208);
      SetColor(10);
      OutTextXY(xc-310,yc+165,' + => Ingrandisce     Grafico');
      OutTextXY(xc-310,yc+180,' - => Rimpicciolisce  Grafico');
      OutTextXY(xc-310,yc+195,'USCITA  con Qualsiasi Freccia');
      SetColor(14);
    rsole:=3;
    SetLineStyle(3,0,1);
    line(xc,yc,2*xc,yc);
    FillEllipse(xc,yc,round(rsole),round(rsole));

    FillEllipse(xc+round(ka*xeclpsave),yc-round(ka*yeclpsave),2,2); {Spot asteroide}
    FillEllipse(xc+round(ka*xtsave),yc-round(ka*ytsave),2,2);       {Spot Terra}


       {PLOT Orbita}
    colore:=14;
    k:=0.0;
repeat
    xeclp:=ka * p / (1 + e * COS(k)) * (r11 * COS(k) + r12 * SIN(k));
    yeclp:=ka * p / (1 + e * COS(k)) * (r21 * COS(k) + r22 * SIN(k));
    zeclp:=ka * p / (1 + e * COS(k)) * (r31 * COS(k) + r32 * SIN(k));
            x:=xc+round(xeclp);  y:=yc-round(yeclp);
              if zeclp>0 then PutPixel(x,y,10) ELSE PutPixel(x,y,12);
    xt:=ka * pt / (1 + e1 * COS(k)) * (Rr11 * COS(k) + Rr12 * SIN(k));
    yt:=ka * pt / (1 + e1 * COS(k)) * (Rr21 * COS(k) + Rr22 * SIN(k));
    zt:=ka * pt / (1 + e1 * COS(k)) * (Rr31 * COS(k) + Rr32 * SIN(k));
            x:=xc+round(xt);  y:=yc-round(yt);
              PutPixel(x,y,11);
    xm:=ka * pm / (1 + em * COS(k)) * (rm11 * COS(k) + rm12 * SIN(k));
    ym:=ka * pm / (1 + em * COS(k)) * (rm21 * COS(k) + rm22 * SIN(k));
    zm:=ka * pm / (1 + em * COS(k)) * (rm31 * COS(k) + rm32 * SIN(k));
            x:=xc+round(xm);  y:=yc-round(ym);
              if zm>0 then PutPixel(x,y,14) ELSE PutPixel(x,y,15);
    xg:=ka * pg / (1 + eg * COS(k)) * (rg11 * COS(k) + rg12 * SIN(k));
    yg:=ka * pg / (1 + eg * COS(k)) * (rg21 * COS(k) + rg22 * SIN(k));
    zg:=ka * pg / (1 + eg * COS(k)) * (rg31 * COS(k) + rg32 * SIN(k));
            x:=xc+round(xg);  y:=yc-round(yg);
              if zg>0 then PutPixel(x,y,13) ELSE PutPixel(x,y,11);
    k:=k+rad/1;  {......per infittire i punti inserire 2 o 3 al denominatore}
 until k>2*pi;

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
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
Skip:
     ColoriSchermo;
     Titolo;
   ASSIGN (fb,'AST_97.DAT');      {File Dati}
   Reset(fb);
   ncat:=''; nome:=''; c3:=''; c4:=''; c5:=''; c6:='';
   c7:=''; c8:='';  c9:=''; c10:=''; c11:='';
repeat gotoxy(28,4); write('N_cat [0 per uscire]= '); readln(zoz);
      until (zoz<=6678) and (zoz>=0);
   if zoz=0 then halt;
    For j:=1 to LUNG do begin
     seek(fb,LUNG*(zoz-1)+j-1);
     read(fb,byt);
        if (j>=1) and (j<=5) then ncat:=ncat+Chr(byt);
        if (j>=7) and (j<=24) then nome:=nome+Chr(byt);
        if (j>=25) and (j<=29) then c3:=c3+Chr(byt); VAL(c3,Hmag,xxx);
        if (j>=32) and (j<=36) then c4:=c4+Chr(byt); VAL(c4,Gmag,xxx);
        if (j>=40) and (j<=48) then c5:=c5+Chr(byt); VAL(c5,m,xxx);
        if (j>=52) and (j<=60) then c6:=c6+Chr(byt); VAL(c6,om,xxx);
        if (j>=64) and (j<=72) then c7:=c7+Chr(byt); VAL(c7,lnodo,xxx);
        if (j>=76) and (j<=83) then c8:=c8+Chr(byt); VAL(c8,incl,xxx);
        if (j>=87) and (j<=95) then c9:=c9+Chr(byt); VAL(c9,e,xxx);
        if (j>=99) and (j<=108) then c10:=c10+Chr(byt);  VAL(c10,n,xxx);
        if (j>=111) and (j<=120) then c11:=c11+Chr(byt); VAL(c11,a,xxx);
    end;
        wo:=ncat+' '+nome;
    Close(fb);
     ElementiOrbitali;
     CALCOLO;
     GraficoConica;
   goto Skip;
  End.

   {FINE Programma}
