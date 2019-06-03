Program FILTROAS;   {EFFEMERIDI Asteroidi per screening à,ë
                     Catalogo.dat - Epoca 1997 12 18.0 }

  Uses  Crt;
   Const LUNG=122; rad=pi/180; kgauss=0.01720209895;
   Var
      wo: String[24];
      nomecat: String[12];
      Hmag,Gmag,m,om,lnodo,incl,e,n,a,
      g,mm,yy,
      r11,r21,r31,r12,r22,r32,p,
       Rr11,Rr21,Rr31,Rr12,Rr22,Rr32,pt,e1,xt,yt,zt,
      jd,aw       : Real;
      fb          : file of Byte;
      ft          : Text;
      xxx,u,band  : Integer;
      k,j,zoz     : LongInt;
      byt         : Byte;
      ncat,c3,c4  : String[5];
      nome        : String[18];
      c5,c9,c7,c6 : String[9];
      c8          : String[8];
      c10,c11     : String[10];
      as          : String;
      AR,DE, arinf,AR0,arsup,deinf,DE0, desup,
      ARore,ARmin,DEgra,DEpri,oredif,mindif,orestep,grastep,
      dedif,maglim: Real;
      cont,kount: Integer;
      aaa: String;

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
     gotoxy(19,2); writeln('FILTRO DI ASTEROIDI CATALOGATI (1..6678)');
     gotoxy(19,3); writeln('----------------------------------------');
     textcolor(14);
  End;

   Procedure KEPLER(m,e:Real; VAR v:Real);
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

   Procedure SAL;    { Stato Avanzamento Lavori }
   Begin
  textcolor(12); gotoxy(16,9); write(wo); textcolor(14);
                 gotoxy(16,10); write('Completato: ',zoz/kount*100:4:1,' %');
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
repeat gotoxy(27,5);write(' Giorno.dec = ');readln(g);
      until (g<32) and (g>=0);
repeat gotoxy(27,6);write('       Mese = ');readln(m);
      until m in [1..12];
       gotoxy(27,7);write('       Anno = ');readln(y);

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
      jd:=b+c+d+g+1720994.5;
     end;

   Procedure CALCOLO;
   Const JDt0=2450800.5; JD2000=2451545.0;
   Var
       M1,a1,om1,lnodo1,incl1,
       M2,a2,e2,om2,lnodo2,incl2,
       Mast,Mter,te1,te2,
       x,y,z,dx,dy,dz,deltap,rp,
       yeq,zeq,yteq,zteq,CO,SI,
       angfase,tg2,aux1,aux2,gi,me,an,
       magn,magn0,magn1,magn2,
       rs
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
 { Anomalie MEDIE }
    Mast:= M2 + kgauss / (a2 *SQRT(a2)) * (jd-JDt0);
    Mter:= -2.470899 + 0.985609113 * (jd-JD2000);
    Mast:= MODULO(Mast, 2 * pi); Mter:= rad*MODULO(Mter, 360);
    Kepler(Mter,e1,te1);
    Kepler(Mast,e2,te2);
    x:= p / (1 + e2 * COS(te2)) * (r11 * COS(te2) + r12 * SIN(te2));
    y:= p / (1 + e2 * COS(te2)) * (r21 * COS(te2) + r22 * SIN(te2));
    z:= p / (1 + e2 * COS(te2)) * (r31 * COS(te2) + r32 * SIN(te2));

    xt:= pt / (1 + e1 * COS(te1)) * (Rr11 * COS(te1) + Rr12 * SIN(te1));
    yt:= pt / (1 + e1 * COS(te1)) * (Rr21 * COS(te1) + Rr22 * SIN(te1));
    zt:= pt / (1 + e1 * COS(te1)) * (Rr31 * COS(te1) + Rr32 * SIN(te1));

{ Cambio di coordinate: da eclittiche ad equatoriali}
  CO:= COS(23.4392911 * rad); SI:= SIN(23.4392911 * rad);
  yeq:= y * CO - z * SI;       yteq:= yt * CO - zt * SI;
  zeq:= y * SI + z * CO;       zteq:= yt * SI + zt * CO;
        y:=yeq;  z:=zeq;
       yt:=yteq;  zt:=zteq;
 {Fine Cambio di coordinate}

    dx:= x - xt; dy:= y - yt; dz:= z - zt;
    deltap:= SQRT(SQR(dx) + SQR(dy) + SQR(dz));
    rp:= SQRT(SQR(x) + SQR(y) + SQR(z));
    rs:= SQRT(SQR(xt) + SQR(yt) + SQR(zt));

     AR:= ArcoTANG(dy, dx)/(15*rad);
     DE:= ArcoSIN(dz/deltap)/rad;

 { Angolo di Fase di P }
   angfase:= ArcoCOS((x*dx+y*dy+z*dz)/(rp*deltap));

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

  if (AR>=arinf) and (AR<=arsup) and (DE>=deinf) and (DE<=desup)
     and (magn<=maglim) then begin
  cont:=cont+1;
textcolor(14);
textcolor(12); gotoxy(21,19); write(wo); write(ft,' '+ wo);textcolor(14);
gotoxy(46,19); write('m= ',magn:4:1); write(ft,'   m= ',magn:4:1);
gotoxy(15,11); write('Corpi Selez:',cont:3);
gotoxy(32,21); textcolor(6);writeln('Equinozio J2000');
 textcolor(15); gotoxy(23,20);
band:=1; aw:=AR; Stringa;  write('à= ',as); write(ft,'   à= ',as);
band:=2; aw:=DE;      Stringa;  writeln('   ë= ',as); writeln(ft,'   ë= ',as);
 textcolor(14);
gotoxy(31,22); writeln('Magn.limite= ',maglim:4:1);
  end;

   End;

  Procedure Setaccio;
   Begin;
 gotoxy(32,13);write('Centro Lastra');
 gotoxy(25,14);write('AR:   Ore = ');gotoxy(37,14); readln(ARore);
 gotoxy(42,14);write('Min = ');gotoxy(48,14); readln(ARmin);
  AR0:=ARore+ARmin/60;
 gotoxy(25,15);write('DE: Gradi.dec = ');gotoxy(41,15); readln(DE0);

 gotoxy(30,16);write('Larghezza Campo');
 gotoxy(25,17);write('in AR:   Ore = ');gotoxy(40,17); readln(oredif);
 gotoxy(42,17);write('Min = ');gotoxy(48,17); readln(mindif);
   orestep:=oredif+mindif/60;
   arinf:=AR0-orestep/2; arsup:=AR0+orestep/2;
 gotoxy(25,18);write('in DE: Gradi.dec = ');gotoxy(44,18); readln(dedif);
   grastep:=dedif;
   deinf:=DE0-grastep/2; desup:=DE0+grastep/2;

 gotoxy(25,19);write('Magnitudine LIMITE = ');gotoxy(46,19); readln(maglim);
 gotoxy(28,21);write('Nome Catalogo= ');gotoxy(43,21); readln(nomecat);

 clrscr;
  gotoxy(42,9); write('AR0=',AR0:7:3,'  DE0=',DE0:7:3);
  gotoxy(42,10); write('Limiti in AR->',arinf:7:3,' ö',arsup:7:3);
  gotoxy(42,11); write('Limiti in DE->',deinf:7:3,' ö',desup:7:3);
  End;
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
    ColoriSchermo;
    Titolo;
    GiornoGiuliano;
    Setaccio;

    kount:=0; cont:=0;
   ASSIGN (ft,nomecat);      {Conta i Records del Catalogo}
   Reset(ft);
 while not eof(ft) do begin
    kount:=kount+1;
    readln(ft,aaa);
 end;
   Close (ft);

   ASSIGN (fb,nomecat);      {File Dati}
   Reset(fb);

   ASSIGN (ft,'FILTROAS.RIS');      {File Risultati}
   Rewrite(ft);

 for zoz:=1 to kount do begin
   ncat:=''; nome:=''; c3:=''; c4:=''; c5:=''; c6:='';
   c7:=''; c8:='';  c9:=''; c10:=''; c11:='';

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

     SAL;
     CALCOLO;
  end;
    Close(fb);
    Close(ft);
  End.

   {FINE Programma}
